#!/usr/bin/env python3
"""Sequence mapping utilities supporting minimap2 and BLAST with a simulator fallback.

Features:
- map_sequences_to_genomes(queries: dict, genomes: dict, min_identity, min_coverage)
  -> in-memory exact/subsequence mapping used for unit tests and fallback.
- CLI for file-based mapping: build indices (optional), run minimap2/blastn if available, or fall back to simulator.
- Output JSON lines / TSV with fields: query, mapped, chromosome, position, identity, coverage, tool, score, notes

This module is intentionally conservative (no external dependency hard fail) and provides
usable mapping for small datasets and unit tests.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any


# Simple, dependency-free progress indicator
class SimpleProgress:
    def __init__(self, total: int, label: str = '', unit: str = 'items'):
        self.total = total
        self.label = label
        self.unit = unit
        self.start_time = time.time()
        self.last = 0

    def update(self, current: int):
        now = time.time()
        elapsed = now - self.start_time
        if current <= 0:
            rate = 0
            eta = 0
        else:
            rate = current / elapsed if elapsed > 0 else 0
            eta = (self.total - current) / rate if rate > 0 else 0
        pct = (current / self.total * 100) if self.total > 0 else 100
        bar_len = 30
        filled = int(bar_len * current / self.total) if self.total > 0 else bar_len
        bar = '[' + '=' * filled + ' ' * (bar_len - filled) + ']'
        sys.stderr.write(f"\r{self.label} {bar} {current}/{self.total} {self.unit} ({pct:.1f}%) Elapsed: {elapsed:.1f}s ETA: {eta:.1f}s   ")
        sys.stderr.flush()

    def finish(self):
        self.update(self.total)
        sys.stderr.write("\n")
        sys.stderr.flush()


# Very small spinner for long-running subprocesses
class Spinner:
    def __init__(self, label: str = 'Working'):
        self.chars = ['|', '/', '-', '\\']
        self.i = 0
        self.label = label
        self.last_write = 0

    def spin(self):
        c = self.chars[self.i % len(self.chars)]
        self.i += 1
        sys.stderr.write(f"\r{self.label} {c}   ")
        sys.stderr.flush()

    def done(self):
        sys.stderr.write('\r')
        sys.stderr.flush()


# Progress file helper
def _write_progress_file(progress_path: Optional[Path], stage: str, current: int, total: int):
    if not progress_path:
        return
    payload = {
        'stage': stage,
        'current': current,
        'total': total,
        'percent': float(current) / total if total else 100.0,
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    }
    try:
        # atomic write if available
        try:
            import atomic_io
            atomic_io.write_json_atomic(progress_path, payload, indent=2)
        except Exception:
            tmp = progress_path.with_suffix('.progress.tmp')
            tmp.write_text(json.dumps(payload))
            tmp.replace(progress_path)
    except Exception:
        pass


def map_sequences_to_genomes(queries: Dict[str, str], genomes: Dict[str, str], min_identity: float = 0.9, min_coverage: float = 0.5, fast_mode: bool = True, k: int = 15, max_contig_len: int = 5_000_000, max_positions_per_kmer: int = 20, progress_path: Optional[Path] = None) -> Dict[str, Dict[str, Any]]:
    """Map query sequences to genomes.

    If fast_mode=True, use a hash-based k-mer index (k-mer size `k`) for fast approximate matching.
    Very large contigs (> max_contig_len) are trimmed (head+tail) to prevent extremely slow scans.

    Returns mapping: {query_id: {mapped: bool, chromosome, position, identity, coverage, tool, score, notes}}
    """
    results: Dict[str, Dict[str, Any]] = {}

    # Preprocess genomes: uppercase and trim large contigs
    processed: Dict[str, str] = {}
    for chrom, seq in genomes.items():
        s = (seq or '').upper()
        if max_contig_len and len(s) > max_contig_len:
            s = s[:max_contig_len//2] + s[-max_contig_len//2:]
        processed[chrom] = s

    # Fast k-mer index mode
    if fast_mode:
        if k <= 0:
            k = 15
        kmer_index: Dict[str, List[Tuple[str, int]]] = {}
        # Build k-mer index (show progress)
        chroms = list(processed.items())
        prog = SimpleProgress(total=len(chroms), label='Indexing')
        for idx_chrom, (chrom, seq) in enumerate(chroms, start=1):
            L = len(seq)
            if L < k:
                prog.update(idx_chrom)
                _write_progress_file(progress_path, 'indexing', idx_chrom, len(chroms))
                continue
            # sample positions with stride for very large sequences to reduce memory
            stride = 1 if L < 1_000_000 else max(1, L // 500_000)
            for i in range(0, L - k + 1, stride):
                kmer = seq[i:i+k]
                lst = kmer_index.setdefault(kmer, [])
                if len(lst) < max_positions_per_kmer:
                    lst.append((chrom, i))
            prog.update(idx_chrom)
            _write_progress_file(progress_path, 'indexing', idx_chrom, len(chroms))
        prog.finish()
        _write_progress_file(progress_path, 'indexing', len(chroms), len(chroms))

        # Map queries (show progress)
        q_items = list(queries.items())
        prog_q = SimpleProgress(total=len(q_items), label='Mapping')
        for idx_q, (qid, qseq) in enumerate(q_items, start=1):
            qseq = (qseq or '').upper()
            qlen = len(qseq)
            if qlen < k:
                results[qid] = {
                    'mapped': False,
                    'chromosome': None,
                    'position': None,
                    'identity': 0.0,
                    'coverage': 0.0,
                    'tool': 'simulator-fast',
                    'score': 0,
                    'notes': 'query shorter than k'
                }
                prog_q.update(idx_q)
                _write_progress_file(progress_path, 'mapping', idx_q, len(q_items))
                continue

            counts: Dict[Tuple[str, int], int] = {}
            # use stride over query kmers to reduce runtime
            qstep = max(1, k // 2)
            for i in range(0, qlen - k + 1, qstep):
                kmer = qseq[i:i+k]
                if kmer not in kmer_index:
                    continue
                for chrom, pos in kmer_index[kmer]:
                    key = (chrom, pos - i)
                    counts[key] = counts.get(key, 0) + 1

            if not counts:
                results[qid] = {
                    'mapped': False,
                    'chromosome': None,
                    'position': None,
                    'identity': 0.0,
                    'coverage': 0.0,
                    'tool': 'simulator-fast',
                    'score': 0,
                    'notes': 'no k-mer matches'
                }
                prog_q.update(idx_q)
                _write_progress_file(progress_path, 'mapping', idx_q, len(q_items))
                continue

            (best_chrom, best_offset), best_count = max(counts.items(), key=lambda x: x[1])
            coverage = min(1.0, (best_count * k) / max(1, qlen))
            identity = 1.0 if coverage >= min_coverage else coverage
            position = max(1, best_offset + 1)
            results[qid] = {
                'mapped': True,
                'chromosome': best_chrom,
                'position': position,
                'identity': identity,
                'coverage': coverage,
                'tool': 'simulator-fast',
                'score': best_count * k,
                'notes': f'kmer_count={best_count} k={k} max_positions={max_positions_per_kmer}'
            }
            prog_q.update(idx_q)
            _write_progress_file(progress_path, 'mapping', idx_q, len(q_items))
        prog_q.finish()
        _write_progress_file(progress_path, 'mapping', len(q_items), len(q_items))

        return results

    # Slow exact + partial match mode (legacy)
    for qid, qseq in queries.items():
        qseq = qseq.strip().upper()
        qlen = len(qseq)
        best = None
        # Simple exact substring search
        for chrom, cseq in processed.items():
            idx = cseq.find(qseq)
            if idx != -1:
                res = {
                    'mapped': True,
                    'chromosome': chrom,
                    'position': idx + 1,
                    'identity': 1.0,
                    'coverage': 1.0,
                    'tool': 'simulator-exact',
                    'score': qlen,
                    'notes': 'exact substring match'
                }
                best = res
                break

        if not best:
            for chrom, cseq in processed.items():
                for window in range(len(qseq), max(int(len(qseq)*min_coverage) - 1, 0), -10):
                    for start in range(0, qlen - window + 1):
                        subseq = qseq[start:start+window]
                        if subseq in cseq:
                            coverage = window / qlen
                            identity = 1.0
                            if coverage >= min_coverage and identity >= min_identity:
                                res = {
                                    'mapped': True,
                                    'chromosome': chrom,
                                    'position': cseq.find(subseq) + 1,
                                    'identity': identity,
                                    'coverage': coverage,
                                    'tool': 'simulator-partial',
                                    'score': window,
                                    'notes': f'partial match window {window}'
                                }
                                if not best or res['score'] > best['score']:
                                    best = res
                    if best:
                        break
                if best:
                    break

        if not best:
            results[qid] = {
                'mapped': False,
                'chromosome': None,
                'position': None,
                'identity': 0.0,
                'coverage': 0.0,
                'tool': 'simulator-none',
                'score': 0,
                'notes': 'no match found'
            }
        else:
            results[qid] = best

    return results


# -------------------
# CLI / file helpers
# -------------------


def is_executable_in_path(name: str) -> bool:
    return shutil.which(name) is not None


def run_minimap2(genome_fasta: Path, query_fasta: Path, out_paf: Path, threads: int = 1, preset: str = 'map-ont', show_progress: bool = False) -> int:
    """Run minimap2 and write a PAF file. Requires minimap2 in PATH.

    If show_progress is True, display a small spinner while the process runs.
    """
    if not is_executable_in_path('minimap2'):
        raise RuntimeError('minimap2 is not available in PATH')
    cmd = ['minimap2', '-t', str(threads), '-x', preset, str(genome_fasta), str(query_fasta)]
    with open(out_paf, 'w', encoding='utf-8') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        spinner = Spinner('minimap2') if show_progress else None
        try:
            while proc.poll() is None:
                if spinner:
                    spinner.spin()
                time.sleep(0.2)
            ret = proc.returncode
            if spinner:
                spinner.done()
            return ret
        finally:
            try:
                # drain stderr if any
                stderr = proc.stderr.read() if proc.stderr else ''
                if stderr:
                    sys.stderr.write('\n'.join(stderr.splitlines()[-5:]) + '\n')
            except Exception:
                pass


def parse_paf_line(line: str) -> Optional[Dict[str, Any]]:
    # PAF columns: qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch?, alnlen?, mapq?
    try:
        cols = line.strip().split('\t')
        qname = cols[0]
        qlen = int(cols[1])
        qstart = int(cols[2])
        qend = int(cols[3])
        strand = cols[4]
        tname = cols[5]
        tlen = int(cols[6])
        tstart = int(cols[7])
        tend = int(cols[8])
        # Parse optional fields for NM (edit distance)
        nm = None
        for c in cols[12:]:
            if c.startswith('NM:i:'):
                nm = int(c.split(':')[-1])
                break
        aln_len = qend - qstart
        identity = ((aln_len - nm) / aln_len) if (nm is not None and aln_len > 0) else None
        coverage = aln_len / qlen if qlen > 0 else 0.0
        return {
            'qname': qname,
            'qlen': qlen,
            'qstart': qstart,
            'qend': qend,
            'strand': strand,
            'tname': tname,
            'tlen': tlen,
            'tstart': tstart,
            'tend': tend,
            'aln_len': aln_len,
            'nm': nm,
            'identity': identity,
            'coverage': coverage
        }
    except Exception:
        return None


def parse_paf(paf_path: Path, show_progress: bool = False) -> Dict[str, List[Dict[str, Any]]]:
    results: Dict[str, List[Dict[str, Any]]] = {}
    total = paf_path.stat().st_size if paf_path.exists() else 0
    prog = SimpleProgress(total=total, label=f'Parsing {paf_path.name}', unit='bytes') if show_progress else None
    with open(paf_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            rec = parse_paf_line(line)
            if not rec:
                continue
            q = rec['qname']
            results.setdefault(q, []).append(rec)
            if prog:
                prog.update(f.tell())
    if prog:
        prog.finish()
    return results


def run_blastn_and_parse(genome_fasta: Path, query_fasta: Path, out_tab: Path, threads: int = 1, show_progress: bool = False) -> int:
    # Prefer blastn from NCBI BLAST+ (blastn) if available; requires makeblastdb prior
    if not is_executable_in_path('blastn'):
        raise RuntimeError('blastn not in PATH')
    # Ensure db exists (use fasta as -dbtype nucl requires makeblastdb)
    db_prefix = genome_fasta.with_suffix('').with_suffix('')
    # Run blastn: outfmt 6 tabular: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    cmd = ['blastn', '-query', str(query_fasta), '-db', str(db_prefix), '-num_threads', str(threads), '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore']
    with open(out_tab, 'w', encoding='utf-8') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        spinner = Spinner('blastn') if show_progress else None
        try:
            while proc.poll() is None:
                if spinner:
                    spinner.spin()
                time.sleep(0.2)
            ret = proc.returncode
            if spinner:
                spinner.done()
            return ret
        finally:
            try:
                stderr = proc.stderr.read() if proc.stderr else ''
                if stderr:
                    sys.stderr.write('\n'.join(stderr.splitlines()[-5:]) + '\n')
            except Exception:
                pass


def cli_main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--genomes', nargs='+', help='Genome FASTA files to map against', required=True)
    p.add_argument('--queries', help='Query FASTA file with sequences to map', required=True)
    p.add_argument('--out', help='Output file (jsonl)', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--method', choices=['minimap2', 'blastn', 'simulator', 'auto'], default='auto')
    p.add_argument('--min-identity', type=float, default=0.9)
    p.add_argument('--min-coverage', type=float, default=0.5)
    p.add_argument('--fast', action='store_true', help='Use fast simulator k-mer mode when using simulator')
    p.add_argument('--k', type=int, default=15, help='k-mer size for fast simulator')
    p.add_argument('--max-contig-bases', type=int, default=5000000, help='Max bases per contig to consider in simulator (head+tail sampling)')
    p.add_argument('--max-positions-per-kmer', type=int, default=20, help='Cap positions stored per k-mer to limit memory')
    p.add_argument('--no-progress', action='store_true', help='Disable progress output (for CI or quiet runs)')
    args = p.parse_args(argv)
    show_progress = not args.no_progress

    genomes = args.genomes
    queries = args.queries
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    # Progress file path (persistent heartbeat for UI/CI)
    progress_path: Optional[Path] = outp.with_suffix('.progress.json') if show_progress else None

    # Read queries into memory (FASTA simple parser)
    queries_dict: Dict[str, str] = {}
    with open(queries, 'r', encoding='utf-8') as f:
        cur = None
        seqs = []
        for line in f:
            if line.startswith('>'):
                if cur:
                    queries_dict[cur] = ''.join(seqs)
                cur = line[1:].strip().split()[0]
                seqs = []
            else:
                seqs.append(line.strip())
        if cur:
            queries_dict[cur] = ''.join(seqs)

    # Read genomes into memory too (for simulator fallback)
    genomes_dict: Dict[str, str] = {}
    for g in genomes:
        name = Path(g).stem
        seq = []
        with open(g, 'r', encoding='utf-8') as f:
            cur = None
            seqs = []
            for line in f:
                if line.startswith('>'):
                    if cur:
                        genomes_dict[cur] = ''.join(seqs)
                    cur = line[1:].strip().split()[0]
                    seqs = []
                else:
                    seqs.append(line.strip())
            if cur:
                genomes_dict[cur] = ''.join(seqs)

    # Decide method
    method = args.method
    if method == 'auto':
        if is_executable_in_path('minimap2'):
            method = 'minimap2'
        elif is_executable_in_path('blastn'):
            method = 'blastn'
        else:
            method = 'simulator'

    results = {}

    if method == 'simulator':
        results = map_sequences_to_genomes(
            queries_dict,
            genomes_dict,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            fast_mode=args.fast,
            k=args.k,
            max_contig_len=args.max_contig_bases,
            max_positions_per_kmer=args.max_positions_per_kmer,
            progress_path=progress_path,
        )
    else:
        # For now, if binaries available, we still run an external tool per genome and parse outputs.
        # Implementing full PAF/BLAST parsing for multi-genome is handled by running each genome as a separate reference.
        for g in genomes:
            gpath = Path(g)
            if method == 'minimap2' and is_executable_in_path('minimap2'):
                paf = outp.with_suffix('.paf')
                try:
                    _write_progress_file(progress_path, 'minimap2', 0, 1)
                    rc = run_minimap2(gpath, Path(queries), paf, threads=args.threads, show_progress=show_progress)
                    _write_progress_file(progress_path, 'minimap2', 1, 1)
                    if rc == 0 and paf.exists():
                        paf_results = parse_paf(paf, show_progress=show_progress)
                        # Convert paf_results to our output schema for each query
                        for qid, recs in paf_results.items():
                            # select best by aln_len
                            best = max(recs, key=lambda r: r.get('aln_len', 0))
                            results[qid] = {
                                'mapped': True,
                                'chromosome': best['tname'],
                                'position': best['tstart'] + 1,
                                'identity': best['identity'] or 0.0,
                                'coverage': best['coverage'],
                                'tool': 'minimap2',
                                'score': best.get('aln_len', 0),
                                'notes': ''
                            }
                except Exception as e:
                    _write_progress_file(progress_path, 'minimap2', 0, 1)
                    print('minimap2 run failed, falling back to simulator:', e, file=sys.stderr)
            elif method == 'blastn' and is_executable_in_path('blastn'):
                tab = outp.with_suffix('.blast.tsv')
                try:
                    _write_progress_file(progress_path, 'blastn', 0, 1)
                    rc = run_blastn_and_parse(Path(g), Path(queries), tab, threads=args.threads, show_progress=show_progress)
                    _write_progress_file(progress_path, 'blastn', 1, 1)
                    # Very simple parser with progress
                    if rc == 0 and tab.exists():
                        total = tab.stat().st_size if tab.exists() else 0
                        prog = SimpleProgress(total=total, label=f'Parsing {tab.name}', unit='bytes') if show_progress else None
                        with open(tab, 'r', encoding='utf-8') as f:
                            for line in f:
                                parts = line.strip().split('\t')
                                if len(parts) < 10:
                                    if prog:
                                        prog.update(f.tell())
                                    continue
                                qid = parts[0]
                                sseqid = parts[1]
                                pident = float(parts[2]) / 100.0
                                length = int(parts[3])
                                qstart = int(parts[4]); qend = int(parts[5])
                                sstart = int(parts[6]); send = int(parts[7])
                                evalue = parts[8]
                                bitscore = float(parts[9])
                                coverage = (qend - qstart + 1) / len(queries_dict.get(qid, '')) if qid in queries_dict else 0.0
                                if qid not in results or results[qid].get('score', 0) < bitscore:
                                    results[qid] = {
                                        'mapped': True,
                                        'chromosome': sseqid,
                                        'position': sstart,
                                        'identity': pident,
                                        'coverage': coverage,
                                        'tool': 'blastn',
                                        'score': bitscore,
                                        'notes': f'evalue={evalue}'
                                    }
                                if prog:
                                    prog.update(f.tell())
                        if prog:
                            prog.finish()
                except Exception as e:
                    _write_progress_file(progress_path, 'blastn', 0, 1)
                    print('blastn run failed, falling back to simulator:', e, file=sys.stderr)
        # Fill unmapped using simulator as fallback
        sim_results = map_sequences_to_genomes(
            queries_dict,
            genomes_dict,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            fast_mode=args.fast,
            k=args.k,
            max_contig_len=args.max_contig_bases,
            max_positions_per_kmer=args.max_positions_per_kmer,
            progress_path=progress_path,
        )
        for qid, r in sim_results.items():
            if qid not in results:
                results[qid] = r

    # Write outputs as JSON lines
    with open(outp, 'w', encoding='utf-8') as out:
        for qid, rec in results.items():
            line = {'query': qid}
            line.update(rec)
            out.write(json.dumps(line) + '\n')

    return 0


if __name__ == '__main__':
    raise SystemExit(cli_main())
