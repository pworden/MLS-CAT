
#!/usr/bin/env python3
"""
mls_khajiit.py — Unified MLST/MLSA phylogenetic pipeline

Combines BLAST extraction, alignment, supermatrix construction, RAxML tree inference,
MinVar rooting, and tseemann/mlst typing into a single script. Input query sequences
can be whole genomes or sequence amplicons.

ARGUMENTS:
    --input-query-sequences       Folder of input query sequences (genomes or amplicons)
    --target-typing-genes         Folder of target gene files
    --output                      Output parent directory
    -e, --query-seq-extension     Input query sequence file extensions (default: fna fasta)
    --type-genes-extension        Allowed target gene file extensions (default: fasta fa fas fna tfa)
    --multi-allele                Target gene files are a multi-allele MLST scheme (best-hit BLAST only)
    --min-genes INT               Minimum genes per sequence to include in tree (default: 2)
    --run-mlst                    Run tseemann/mlst typing on input query sequences
    --custom-mlst-db STR          Path to a custom mlst database directory (see README: Custom MLST schemes)
    --mlst-scheme-name STR        Force a specific mlst scheme name (e.g. yourspecies)
    --mlst-minid FLOAT            MLST minimum allele identity; if unset, uses --min-identity value
    --mlst-mincov FLOAT           MLST minimum coverage for partial alleles (default: mlst default 50%)
    --no-trim-amplicon-input      Skip trimal and use raw MAFFT alignments (recommended for amplicon-only input)
    --bootstraps INT              Number of RAxML bootstrap replicates (default: 100)
    --threads INT                 Number of threads for BLAST, MLST, and RAxML (default: 20)
    --skip-extraction             Skip BLAST extraction stage
    --skip-consolidation          Skip consolidation stage
    --skip-alignment              Skip MAFFT/trimal/supermatrix stages
    --skip-tree                   Skip RAxML/MinVar stages
    -f, --force                   Overwrite existing directories
    -v, --verbose                 Verbose output

NOTE: For a species not in the standard pubmlst database, add the scheme to tseemann/mlst
      first — see the README's "Custom MLST schemes" section.

USAGE:
conda activate MLS_Khajiit

# 1. MLSA tree only — reference genomes, single representative per gene:
mls-khajiit \
    --input-query-sequences /path/to/genomes \
    --target-typing-genes /path/to/custom_mlst_db/yourspecies \
    --output /path/to/output \
    --query-seq-extension fasta \
    --type-genes-extension tfa \
    --multi-allele \
    --min-genes 3 \
    --verbose --force

# 2. Amplicons — MLSA + MLST with custom scheme; use --no-trim-amplicon-input for amplicon-only input:
mls-khajiit \
    --input-query-sequences /path/to/amplicon_consensus_sequences \
    --target-typing-genes /path/to/custom_mlst_db/yourspecies \
    --output /path/to/output \
    --query-seq-extension fasta \
    --type-genes-extension tfa \
    --multi-allele \
    --min-genes 2 \
    --run-mlst \
    --custom-mlst-db /path/to/custom_mlst_db \
    --mlst-scheme-name yourspecies \
    --no-trim-amplicon-input \
    --verbose --force

# 3. Combined reference genomes + amplicons — MLSA tree (trimal appropriate for mixed input):
mls-khajiit \
    --input-query-sequences /path/to/combined_input \
    --target-typing-genes /path/to/custom_mlst_db/yourspecies \
    --output /path/to/output \
    --query-seq-extension fasta \
    --type-genes-extension tfa \
    --multi-allele \
    --min-genes 2 \
    --verbose --force

# 4. Full MLSA + MLST with matched identity thresholds and outgroup rooting:
mls-khajiit \
    --input-query-sequences /path/to/combined_input \
    --target-typing-genes /path/to/custom_mlst_db/yourspecies \
    --output /path/to/output \
    --query-seq-extension fasta \
    --type-genes-extension tfa \
    --multi-allele \
    --min-genes 2 \
    --min-identity 85.0 \
    --mlst-minid 85.0 \
    --mlst-mincov 50.0 \
    --run-mlst \
    --custom-mlst-db /path/to/custom_mlst_db \
    --mlst-scheme-name yourspecies \
    --no-trim-amplicon-input \
    --outgroup /path/to/input/outgroup_sample.fasta \
    --bootstraps 1000 \
    --threads 32 \
    --force \
    --verbose
"""

import os
import sys
import shutil
import subprocess
import glob
import argparse
import logging
import re
import csv

from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus

log = logging.getLogger("mage_pipeline")
pipeline_warnings = []


def setup_logging(log_file, verbose):
    """Configure file and console logging."""
    level = logging.DEBUG if verbose else logging.INFO
    log.setLevel(logging.DEBUG)

    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(logging.Formatter("%(message)s"))

    log.addHandler(fh)
    log.addHandler(ch)
    log.info(f"Log file: {log_file}")


# =============================================================================
# Argument parsing and entry point
# =============================================================================

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Unified MLST/MLSA phylogenetic pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument("--input-query-sequences", required=True,
                          help="Folder of input query sequences (genomes or amplicons)")
    required.add_argument("--target-typing-genes", required=False, default=None,
                          help="Folder of target gene files (required unless --mlst-only)")
    required.add_argument("--output", required=True,
                          help="Output parent directory")

    # Optional arguments
    parser.add_argument("-e", "--query-seq-extension", nargs="+",
                        default=["fna", "fasta"],
                        help="Input query sequence file extensions (default: fna fasta)")
    parser.add_argument("--type-genes-extension", nargs="+",
                        default=["fasta", "fa", "fas", "fna", "tfa"],
                        help="Allowed target gene file extensions (default: fasta fa fas fna tfa)")
    parser.add_argument("--multi-allele", action="store_true",
                        help="Flag that target gene files are a multi-allele MLST scheme. "
                             "BLAST will search all alleles and keep only the single best hit per gene per genome.")
    parser.add_argument("--min-genes", type=int, default=2, metavar="INT",
                        help="Minimum genes per genome to include in tree (default: 0)")
    parser.add_argument("--run-mlst", action="store_true",
                        help="Run tseemann/mlst typing on input query sequences")
    parser.add_argument("--custom-mlst-db", default=None,
                        help="Path to custom mlst database directory (see README: Custom MLST schemes). "
                             "Use for species not in the tseemann/mlst database.")
    parser.add_argument("--mlst-scheme-name", default=None,
                        help="MLST scheme name (e.g. yourspecies). "
                             "Required when --run-mlst or --mlst-only is set.")
    parser.add_argument("--skip-extraction", action="store_true",
                        help="Skip BLAST extraction stage")
    parser.add_argument("--skip-consolidation", action="store_true",
                        help="Skip consolidation stage")
    parser.add_argument("--skip-alignment", action="store_true",
                        help="Skip MAFFT/trimal/supermatrix stages")
    parser.add_argument("--no-trim-amplicon-input", action="store_true",
                        help="Skip trimal trimming and use raw MAFFT alignments for supermatrix. "
                             "Recommended when input sequences are amplicons (already clean, uniform length).")
    parser.add_argument("--min-identity", type=float, default=85.0, metavar="FLOAT",
                        help="Minimum percent identity for BLAST hits (default: 85)")
    parser.add_argument("--mlst-minid", type=float, default=None, metavar="FLOAT",
                        help="MLST minimum allele identity threshold; if None, uses --min-identity value (default: None)")
    parser.add_argument("--mlst-mincov", type=float, default=None, metavar="FLOAT",
                        help="MLST minimum coverage for partial alleles (default: None, uses mlst default 50)")
    parser.add_argument("--skip-tree", action="store_true",
                        help="Skip RAxML/MinVar stages")
    parser.add_argument("--no-minvar-root", action="store_true",
                        help="Skip FastRoot MinVar rooting of the final tree.")
    parser.add_argument("--outgroup", default=None, metavar="FASTA_PATH",
                        help="Path to the fasta file (as named in --input-query-sequences) of the "
                             "sample to use as outgroup. Roots the final tree on that taxon instead "
                             "of MinVar rooting (implies --no-minvar-root; MinVar will not run). "
                             "Writes MAGE_tree_output_OutgroupRooted.final. Errors (and skips) if "
                             "the derived taxon name matches zero or more than one tip.")
    parser.add_argument("--mlst-only", action="store_true",
                        help="Run MLST typing only — skips all MLSA stages. "
                             "--target-typing-genes is not required. Implies --run-mlst.")
    parser.add_argument("--bootstraps", type=int, default=100, metavar="INT",
                        help="Number of RAxML bootstrap replicates (default: 100)")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Overwrite existing directories")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    parser.add_argument("--threads", type=int, default=20, metavar="INT",
                        help="Number of threads for BLAST, MLST, and RAxML (default: 20)")

    return parser.parse_args()


# =============================================================================
# Stage 1: BLAST extraction
# =============================================================================

def blast_extraction(query_seqs_folder, targets_folder, extracted_dir, gene_extensions, min_identity=85.0, num_threads=32):
    """BLAST each target gene against all query sequences."""
    query_seqs = [f for f in os.listdir(query_seqs_folder)
                  if os.path.isfile(os.path.join(query_seqs_folder, f))]
    targets = [f for f in os.listdir(targets_folder)
               if os.path.isfile(os.path.join(targets_folder, f))
               and os.path.splitext(f)[1].lstrip('.') in gene_extensions]

    for target_file in targets:
        shortname = os.path.splitext(target_file)[0]
        for query_seq in query_seqs:
            out_path = os.path.join(extracted_dir, f"{shortname}_{query_seq}")
            result = subprocess.run(
                    ["blastn",
                     "-task", "blastn",
                     "-query", os.path.join(targets_folder, target_file),
                     "-subject", os.path.join(query_seqs_folder, query_seq),
                     "-max_target_seqs", "1",
                     "-evalue", "1e-10",
                     "-perc_identity", str(min_identity),
                     "-num_threads", str(num_threads),
                     "-outfmt", "6 bitscore sseq"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
            # Keep only the single best hit (highest bitscore) across all alleles
            best_seq = None
            best_score = -1
            for line in result.stdout.decode('utf-8').splitlines():
                parts = line.split('\t')
                if len(parts) == 2:
                    try:
                        score = float(parts[0])
                        if score > best_score:
                            best_score = score
                            best_seq = parts[1].replace('-', '')
                    except ValueError:
                        pass
            stderr_text = result.stderr.decode('utf-8').strip()
            if result.returncode != 0:
                log.warning(f"  BLAST non-zero exit ({result.returncode}) for {target_file} vs {query_seq}")
            if stderr_text:
                log.debug(f"  BLAST stderr [{target_file} vs {query_seq}]: {stderr_text}")
            with open(out_path, "w") as out_handle:
                out_handle.write(f">{target_file}_{query_seq}\n")
                if best_seq:
                    out_handle.write(best_seq + '\n')

    log.info(f"  Extracted: {len(targets)} genes x {len(query_seqs)} query sequences")



# =============================================================================
# Stage 2: Consolidation
# =============================================================================

def consolidate_targets(extracted_dir, consolidated_dir, targets_dir, gene_extensions):
    """Consolidate extracted targets into per-gene FASTA files."""
    target_files = [f for f in os.listdir(targets_dir)
                    if os.path.isfile(os.path.join(targets_dir, f))
                    and os.path.splitext(f)[1].lstrip('.') in gene_extensions]
    gene_names = sorted([os.path.splitext(f)[0] for f in target_files])

    total_seqs = 0
    for gene in gene_names:
        input_files = sorted(glob.glob(os.path.join(extracted_dir, f"{gene}_*")))
        if not input_files:
            log.warning(f"  No extracted files for {gene}")
            continue

        output_file = os.path.join(consolidated_dir, f"{gene}_references.fasta")
        seq_count = 0
        empty_count = 0
        with open(output_file, 'w') as outf:
            for input_file in input_files:
                with open(input_file, 'r') as inf:
                    content = inf.read().strip()
                    if not content:
                        empty_count += 1
                        continue
                    lines = content.split('\n')
                    if len(lines) > 1 and lines[0].startswith('>'):
                        outf.write(content + '\n')
                        seq_count += 1
                    elif content.startswith('>'):
                        empty_count += 1
        total_seqs += seq_count
        log.debug(f"  {gene}: {seq_count} sequences, {empty_count} empty")

    log.info(f"  Consolidation complete: {total_seqs} total sequences")


# =============================================================================
# Stage 3: MAFFT alignment
# =============================================================================

def align_genes(consolidated_dir, alignments_dir):
    """Align each gene using MAFFT."""
    ref_files = sorted(glob.glob(os.path.join(consolidated_dir, "*_references.fasta")))
    if not ref_files:
        log.warning("  No consolidated reference files found")
        return

    for ref_file in ref_files:
        gene = os.path.basename(ref_file).replace("_references.fasta", "")
        output_file = os.path.join(alignments_dir, f"{gene}.mafft.aln")
        log.info(f"  Aligning {gene}...")
        with open(output_file, "w") as out:
            result = subprocess.run(
                ["mafft", "--adjustdirection", ref_file],
                stdout=out, stderr=subprocess.PIPE, text=True,
                env=os.environ.copy()
            )
            if result.returncode != 0:
                log.error(f"  MAFFT failed for {gene} (exit {result.returncode})")
            if result.stderr.strip():
                log.debug(f"  MAFFT stderr [{gene}]: {result.stderr.strip()}")


# =============================================================================
# Stage 4: trimal
# =============================================================================

def trim_alignments(alignments_dir):
    """Trim alignments using trimal (both -automated1 and -nogaps versions).
    Generates two sets of trimmed alignments for comparison:
      - *.trimmed.fasta (keeps gaps, uses -automated1)
      - *.trimmed_nogaps.fasta (removes all gaps, uses -nogaps)
    Skips gracefully if trimal is not installed.
    """
    if not shutil.which("trimal"):
        log.warning("  trimal not found in PATH — skipping trimming, using raw MAFFT alignments")
        return

    aln_files = sorted(glob.glob(os.path.join(alignments_dir, "*.mafft.aln")))
    for aln_file in aln_files:
        gene = os.path.basename(aln_file).split(".")[0]
        log.info(f"  Trimming {gene}...")

        # Version 1: -automated1 (preserves gaps)
        output_file = os.path.join(alignments_dir, f"{gene}.trimmed.fasta")
        result = subprocess.run(
            ["trimal", "-in", aln_file, "-out", output_file, "-fasta", "-automated1"],
            stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True,
            env=os.environ.copy()
        )
        if result.returncode != 0:
            log.error(f"  trimal -automated1 failed for {gene} (exit {result.returncode})")
        if result.stderr.strip():
            log.debug(f"  trimal -automated1 stderr [{gene}]: {result.stderr.strip()}")

        # Version 2: -nogaps (removes all gaps, MAGe approach)
        output_file_nogaps = os.path.join(alignments_dir, f"{gene}.trimmed_nogaps.fasta")
        result = subprocess.run(
            ["trimal", "-in", aln_file, "-out", output_file_nogaps, "-fasta", "-nogaps"],
            stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True,
            env=os.environ.copy()
        )
        if result.returncode != 0:
            log.error(f"  trimal -nogaps failed for {gene} (exit {result.returncode})")
        if result.stderr.strip():
            log.debug(f"  trimal -nogaps stderr [{gene}]: {result.stderr.strip()}")


# =============================================================================
# Stage 5: Supermatrix
# =============================================================================

def create_supermatrix(alignments_dir, supermatrix_dir, min_genes, nogaps=False):
    """Create concatenated supermatrix from per-gene alignments.

    Args:
        nogaps: If True, uses *.trimmed_nogaps.fasta (gap-removed version).
               If False, uses *.trimmed.fasta or raw *.mafft.aln (gap-preserving).

    Returns: supermatrix_fasta path
    """
    suffix = "_nogaps" if nogaps else ""

    # Try to find trimmed files (with or without nogaps variant)
    if nogaps:
        trimmed_files = sorted(glob.glob(os.path.join(alignments_dir, "*.trimmed_nogaps.fasta")))
        if trimmed_files:
            aln_files = trimmed_files
            log.info("  Using trimmed_nogaps alignments for supermatrix")
    else:
        trimmed_files = sorted(glob.glob(os.path.join(alignments_dir, "*.trimmed.fasta")))
        if trimmed_files:
            aln_files = trimmed_files
            log.info("  Using trimmed alignments for supermatrix")

    # Fallback to raw MAFFT if no trimmed files
    if not trimmed_files:
        aln_files = sorted(glob.glob(os.path.join(alignments_dir, "*.mafft.aln")))
        log.info("  Using raw MAFFT alignments for supermatrix (trimal not run)")

    if not aln_files:
        sys.exit("No alignments found")

    all_taxa = set()
    gene_seqs = {}
    gene_lengths = {}

    for aln_file in aln_files:
        gene_name = os.path.basename(aln_file).split(".")[0]
        log.info(f"  Processing {gene_name}...")
        aln_lengths = {}
        gene_seqs[gene_name] = {}

        for record in SeqIO.parse(aln_file, "fasta"):
            taxon_id = record.id
            # Strip gene file prefix (e.g. "arfB.fna_" or "arfB.tfa_") to get genome filename
            # Header format from BLAST extraction: {gene_name}.{ext}_{genome_filename}
            if taxon_id.startswith(gene_name + ".") and "_" in taxon_id:
                taxon_id = taxon_id[len(gene_name):].split("_", 1)[1]
            seq_str = str(record.seq)
            gene_seqs[gene_name][taxon_id] = seq_str
            aln_lengths[taxon_id] = len(seq_str)
            all_taxa.add(taxon_id)

        if aln_lengths:
            lengths = list(aln_lengths.values())
            if len(set(lengths)) > 1:
                log.warning(f"  Inconsistent alignment lengths in {gene_name}")
            else:
                gene_lengths[gene_name] = lengths[0]
                log.debug(f"    {gene_name}: {lengths[0]} bp, {len(lengths)} taxa")

    concatenated_seqs = {}
    taxon_gene_counts = {}
    total_genes = len(gene_seqs)

    for taxon_id in sorted(all_taxa):
        seq_parts = []
        genes_found = 0
        for gene_name in sorted(gene_seqs.keys()):
            if taxon_id in gene_seqs[gene_name]:
                seq_parts.append(gene_seqs[gene_name][taxon_id])
                genes_found += 1
            else:
                gap_length = gene_lengths.get(gene_name, 0)
                seq_parts.append("-" * gap_length)
        concatenated_seqs[taxon_id] = "".join(seq_parts)
        taxon_gene_counts[taxon_id] = (genes_found, total_genes)

    concatenated_seqs_with_ratio = {}
    excluded_taxa = []
    for taxon_id, seq in concatenated_seqs.items():
        genes_found, total = taxon_gene_counts[taxon_id]
        if genes_found < min_genes:
            excluded_taxa.append((taxon_id, genes_found, total))
            continue
        taxon_base = taxon_id.rsplit('.', 1)[0] if '.' in taxon_id else taxon_id
        taxon_with_ratio = f"{taxon_base}_{genes_found}genesFrom{total}"
        concatenated_seqs_with_ratio[taxon_with_ratio] = seq

    if excluded_taxa:
        log.info(f"\n  Excluded {len(excluded_taxa)} taxa with <{min_genes} genes:")
        for taxon, found, total in excluded_taxa:
            log.info(f"    {taxon}: {found}/{total} genes")

    supermatrix_fasta = os.path.join(supermatrix_dir, f"supermatrix{suffix}.fasta")
    with open(supermatrix_fasta, "w") as out:
        for taxon, seq in concatenated_seqs_with_ratio.items():
            out.write(f">{taxon}\n{seq}\n")

    log.info("  Final supermatrix sequence lengths:")
    for taxon, seq in sorted(concatenated_seqs_with_ratio.items()):
        log.info(f"    {taxon}: {len(seq)} bp")

    records = [SeqRecord(seq, id=taxon, description="")
               for taxon, seq in concatenated_seqs_with_ratio.items()]
    nex = Nexus.Nexus()
    nex.sequences = records
    output_nex = os.path.join(supermatrix_dir, f"supermatrix{suffix}.nex")
    nex.write_nexus_data(output_nex)
    log.info(f"  Supermatrix created: {output_nex}")

    # Gene presence/absence table (only write once for the non-nogaps version)
    if not nogaps:
        in_tree_taxa = set(concatenated_seqs_with_ratio.keys())
        sorted_genes = sorted(gene_seqs.keys())
        presence_csv = os.path.join(supermatrix_dir, "gene_presence_table.csv")
        with open(presence_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["sample"] + sorted_genes + ["genes_found", "total_genes", "in_tree"])
            for taxon_id in sorted(all_taxa):
                taxon_base = taxon_id.rsplit(".", 1)[0] if "." in taxon_id else taxon_id
                genes_found, total = taxon_gene_counts[taxon_id]
                taxon_with_ratio = f"{taxon_base}_{genes_found}genesFrom{total}"
                in_tree = "yes" if taxon_with_ratio in in_tree_taxa else "no"
                row = [taxon_base]
                for gene in sorted_genes:
                    row.append("found" if taxon_id in gene_seqs[gene] else "missing")
                row += [genes_found, total, in_tree]
                writer.writerow(row)
        log.info(f"  Gene presence table: {presence_csv}")

    return supermatrix_fasta


# =============================================================================
# Stage 6: RAxML + MinVar rooting
# =============================================================================

def minvar_root_tree(tree_file, output_unrooted, output_minvar):
    """Root a tree using FastRoot while preserving bootstrap values.

    FastRoot finds the optimal root position that minimizes variance in
    branch lengths. This version ensures bootstrap values are preserved
    through the rooting process using BioPython tree manipulation.
    """
    try:
        # Read tree with BioPython to preserve bootstrap values
        trees = Phylo.read(tree_file, "newick")

        # Write unrooted version (with bootstrap values intact)
        Phylo.write(trees, output_unrooted, "newick")
        log.info(f"  Unrooted tree: {output_unrooted}")

        # Use FastRoot to find optimal root (MV = minimum variance, default)
        # FastRoot reads the tree, finds the best root position, and outputs rooted tree
        result = subprocess.run(
            ["FastRoot.py", "-i", tree_file, "-o", output_minvar, "-s", "newick"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env=os.environ.copy()
        )

        if result.returncode != 0:
            log.error(f"  FastRoot failed (exit {result.returncode})")
            if result.stderr.strip():
                log.debug(f"  FastRoot stderr: {result.stderr.strip()}")
        else:
            # Read the rooted tree from FastRoot
            rooted_tree = Phylo.read(output_minvar, "newick")

            # Read the original unrooted tree to extract bootstrap values
            original_tree = Phylo.read(tree_file, "newick")

            # Transfer bootstrap values from original tree to rooted tree
            # (internal node labels in BioPython are stored as node.name)
            def transfer_bootstraps(orig_clade, root_clade, bootstrap_map=None):
                """Recursively match clades and transfer bootstrap values."""
                if bootstrap_map is None:
                    bootstrap_map = {}

                # Create a signature for each clade (sorted leaf names under it)
                def get_clade_signature(clade):
                    if clade.is_terminal():
                        return frozenset([clade.name])
                    leaves = set()
                    for terminal in clade.get_terminals():
                        leaves.add(terminal.name)
                    return frozenset(leaves)

                orig_sig = get_clade_signature(orig_clade)
                if orig_sig not in bootstrap_map and orig_clade.name:
                    bootstrap_map[orig_sig] = orig_clade.name

                # Recursively process subtrees
                if not orig_clade.is_terminal():
                    for child in orig_clade.clades:
                        transfer_bootstraps(child, root_clade, bootstrap_map)

                return bootstrap_map

            # Build bootstrap value mapping from original tree
            bootstrap_map = transfer_bootstraps(original_tree.root, rooted_tree.root)

            # Apply bootstrap values to rooted tree
            def apply_bootstraps(clade, bootstrap_map):
                if not clade.is_terminal():
                    sig = frozenset([t.name for t in clade.get_terminals()])
                    if sig in bootstrap_map and bootstrap_map[sig]:
                        clade.name = bootstrap_map[sig]
                    for child in clade.clades:
                        apply_bootstraps(child, bootstrap_map)

            apply_bootstraps(rooted_tree.root, bootstrap_map)

            # Write rooted tree back, preserving bootstrap values
            Phylo.write(rooted_tree, output_minvar, "newick")
            log.info(f"  MinVar rooted tree (FastRoot): {output_minvar}")
            if result.stdout.strip():
                log.debug(f"  FastRoot output: {result.stdout.strip()}")

    except Exception as e:
        log.error(f"  Error in MinVar rooting: {e}", exc_info=True)


def outgroup_root_tree(tree_file, output_file, outgroup_fasta_path):
    """Root a tree on the taxon corresponding to outgroup_fasta_path.

    outgroup_fasta_path is the fasta file (as it sits in --input-query-sequences)
    of the sample to use as outgroup, e.g.
    ".../X.arboricola_CFBP_2528_GCF_001013475.1.fasta". The taxon name in the
    final tree is that filename with its extension stripped, plus a
    "_<N>genesFrom<total>" suffix appended by create_supermatrix() — so we
    match any tip whose name starts with "<basename>_".

    This function re-roots the same in-memory
    BioPython tree object directly, so bootstrap/confidence values are
    preserved automatically.
    """
    try:
        basename = os.path.splitext(os.path.basename(outgroup_fasta_path))[0]
        tree = Phylo.read(tree_file, "newick")
        matches = [t for t in tree.get_terminals()
                   if (t.name or "").startswith(basename + "_") or t.name == basename]
        if not matches:
            log.error(f"  No taxon matches outgroup sample '{basename}' — "
                      f"skipping outgroup rooting")
            return
        if len(matches) > 1:
            log.error(f"  Outgroup sample '{basename}' matches {len(matches)} taxa "
                      f"(ambiguous): {[m.name for m in matches]} — skipping outgroup rooting")
            return
        tree.root_with_outgroup(matches[0])
        Phylo.write(tree, output_file, "newick")
        log.info(f"  Outgroup-rooted tree (outgroup: {matches[0].name}): {output_file}")
    except Exception as e:
        log.error(f"  Error in outgroup rooting: {e}", exc_info=True)


def run_raxml(supermatrix_fasta, supermatrix_dir, tree_dir, num_threads=32, bootstraps=1000, nogaps=False,
              no_minvar_root=False, outgroup=None):
    """Run RAxML best tree search, bootstrap, and consensus; then apply MinVar rooting.

    Args:
        nogaps: If True, uses nogaps version of tree directory and output files.
    """
    suffix = "_nogaps" if nogaps else ""
    taxon_mapping = {}
    try:
        records = list(SeqIO.parse(supermatrix_fasta, "fasta"))
        if len(records) < 4:
            msg = (f"Too few sequences in supermatrix ({len(records)}) — RAxML requires "
                   f"at least 4 taxa to infer a meaningful tree. Tree stage skipped.")
            log.error(f"  {msg}")
            pipeline_warnings.append(msg)
            return
        for i, record in enumerate(records, 1):
            short_name = f"taxon_{i}"
            taxon_mapping[short_name] = record.id

        alignment_file = os.path.join(supermatrix_dir, "supermatrix_short_names.phy")
        for i, record in enumerate(records, 1):
            record.id = f"taxon_{i}"
            record.description = ""
        SeqIO.write(records, alignment_file, "phylip")

        mapping_file = os.path.join(supermatrix_dir, "taxon_mapping.txt")
        with open(mapping_file, "w") as mf:
            mf.write("Short_Name\tOriginal_Name\n")
            for short, orig in taxon_mapping.items():
                mf.write(f"{short}\t{orig}\n")
        log.debug(f"  PHYLIP: {alignment_file}")
        log.debug(f"  Mapping: {mapping_file}")
    except Exception as e:
        log.error(f"  Error creating PHYLIP: {e}", exc_info=True)
        sys.exit(1)

    # Rapid bootstrap analysis (-f a = ML search + bootstrap in one pass)
    log.info(f"  Running RAxML rapid bootstrap ({bootstraps} replicates)...")
    alignment_relpath = os.path.relpath(alignment_file, tree_dir)
    result = subprocess.run(
        ["raxmlHPC", "-T", str(num_threads), "-f", "a", "-p", "12345", "-x", "12345",
         "-m", "GTRGAMMA", "-s", alignment_relpath, "-N", str(bootstraps), "-n", "bootstrap"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        cwd=tree_dir, env=os.environ.copy()
    )
    if result.returncode != 0:
        log.error(f"  RAxML rapid bootstrap failed (exit {result.returncode})")
    if result.stderr.strip():
        log.debug(f"  RAxML rapid bootstrap stderr: {result.stderr.strip()}")
    if result.returncode == 0:
        log.info("  Rapid bootstrap analysis completed")

    log.debug("  RAxML output files:")
    for f in sorted(os.listdir(tree_dir)):
        log.debug(f"    {f}")

    # Find best available tree (from rapid bootstrap analysis)
    bootstrap_tree = os.path.join(tree_dir, "RAxML_bipartitions.bootstrap")
    result_tree_file = os.path.join(tree_dir, "RAxML_result.bootstrap")

    tree_to_use = None
    tree_source = None
    if os.path.exists(bootstrap_tree):
        tree_to_use = bootstrap_tree
        tree_source = "rapid bootstrap"
    elif os.path.exists(result_tree_file):
        tree_to_use = result_tree_file
        tree_source = "result tree"

    final_tree = os.path.join(tree_dir, f"MAGE_tree_output{suffix}.final")

    if tree_to_use:
        # Read tree with BioPython to preserve bootstrap values
        tree = Phylo.read(tree_to_use, "newick")

        # Replace short taxon names with original names (preserving tree structure)
        def replace_taxon_names(clade, taxon_mapping):
            if clade.is_terminal():
                if clade.name in taxon_mapping:
                    clade.name = taxon_mapping[clade.name]
            else:
                for child in clade.clades:
                    replace_taxon_names(child, taxon_mapping)

        replace_taxon_names(tree.root, taxon_mapping)

        # Write final tree preserving bootstrap values
        Phylo.write(tree, final_tree, "newick")
        log.info(f"  Final tree ({tree_source}): {final_tree}")

        # --outgroup implies --no-minvar-root: MinVar never runs if an outgroup was given
        if outgroup:
            if no_minvar_root:
                log.info("  Skipping MinVar rooting (--outgroup given; --no-minvar-root also set)")
            else:
                log.info("  Skipping MinVar rooting (--outgroup given takes precedence)")
        elif no_minvar_root:
            log.info("  Skipping MinVar rooting (--no-minvar-root)")
        else:
            log.info("  Applying MinVar rooting...")
            unrooted = os.path.join(tree_dir, f"MAGE_tree_output_unrooted{suffix}.final")
            minvar_rooted = os.path.join(tree_dir, f"MAGE_tree_output_MinVarRooted{suffix}.final")
            minvar_root_tree(final_tree, unrooted, minvar_rooted)

        if outgroup:
            log.info(f"  Applying outgroup rooting (sample: {os.path.basename(outgroup)})...")
            outgroup_rooted = os.path.join(tree_dir, f"MAGE_tree_output_OutgroupRooted{suffix}.final")
            outgroup_root_tree(final_tree, outgroup_rooted, outgroup)
    else:
        log.warning("  No RAxML tree files found")


# =============================================================================
# Stage 7: MLST typing (tseemann/mlst)
# =============================================================================

def build_per_sample_fastas_from_extracted(extracted_dir, query_seq_files, mlst_input_dir):
    """Build per-sample FASTA files for MLST input from Stage 1 extracted targets.

    Each output FASTA contains one sequence per gene (the best BLAST hit already
    selected by Stage 1), giving tseemann/mlst clean, unambiguous single-allele input.
    Returns list of paths to the created FASTA files.
    """
    os.makedirs(mlst_input_dir, exist_ok=True)
    output_files = []

    for query_seq_path in sorted(query_seq_files):
        query_seq_file = os.path.basename(query_seq_path)
        extracted_files = sorted(glob.glob(os.path.join(extracted_dir, f"*_{query_seq_file}")))

        if not extracted_files:
            log.debug(f"  No extracted targets found for {query_seq_file} — skipping")
            continue

        sample_fasta = os.path.join(mlst_input_dir, query_seq_file)
        seqs_written = 0
        with open(sample_fasta, 'w') as out:
            for extracted_file in extracted_files:
                with open(extracted_file, 'r') as f:
                    content = f.read().strip()
                if not content:
                    continue
                lines = content.split('\n')
                if len(lines) >= 2 and lines[0].startswith('>') and lines[1].strip():
                    out.write(content + '\n')
                    seqs_written += 1

        if seqs_written > 0:
            output_files.append(sample_fasta)
            log.debug(f"  {query_seq_file}: {seqs_written} gene sequences written")
        else:
            os.remove(sample_fasta)
            log.debug(f"  {query_seq_file}: no sequences — skipped")

    log.info(f"  Built {len(output_files)} per-sample FASTA files for MLST input")
    return output_files


def run_mlst_typing(input_query_sequences, output, query_seq_extension,
                    custom_mlst_db, mlst_scheme_name, input_files=None,
                    mlst_minid=None, mlst_mincov=None, num_threads=32):
    """Run tseemann/mlst on input query sequences.

    If input_files is provided, those files are used directly (e.g. per-sample FASTAs
    built from Stage 1 extracted targets). Otherwise, files are collected from
    input_query_sequences using query_seq_extension.

    Args:
        mlst_minid: MLST minimum allele identity (default None = mlst default 95)
        mlst_mincov: MLST minimum coverage for partial alleles (default None = mlst default 50)
    """
    if input_files is not None:
        seq_files = sorted(input_files)
    else:
        seq_files = []
        for ext in query_seq_extension:
            seq_files.extend(glob.glob(os.path.join(input_query_sequences, f"*.{ext}")))
        seq_files = sorted(set(seq_files))

    if not seq_files:
        log.warning(f"  No sequence files found in {input_query_sequences}")
        return

    output_tsv = os.path.join(output, "mlst_results.tsv")

    # Build mlst command — prefer mlst directly on PATH (e.g. installed into the
    # same env as this script); fall back to a separate conda env named "mlst"
    # if it isn't.
    if shutil.which("mlst"):
        cmd = ["mlst"]
        log.info("  Using mlst found on PATH")
    else:
        cmd = ["conda", "run", "-n", "mlst", "mlst"]
        log.info("  mlst not found on PATH — using conda env 'mlst'")

    if custom_mlst_db:
        blast_db = os.path.join(custom_mlst_db, "blast", "mlst")
        cmd += ["--datadir", custom_mlst_db, "--blastdb", blast_db]
        log.info(f"  Using custom MLST database: {custom_mlst_db}")

    if mlst_scheme_name:
        # Normalize to tseemann/mlst convention: lowercase, no underscores/hyphens/spaces
        # e.g. Your_Species -> yourspecies
        mlst_scheme_name_normalized = re.sub(r'[^a-z0-9]', '', mlst_scheme_name.lower())
        cmd += ["--scheme", mlst_scheme_name_normalized]
        log.info(f"  Using scheme: {mlst_scheme_name_normalized}")

    if mlst_minid is not None:
        cmd += ["--minid", str(mlst_minid)]
        log.info(f"  MLST minimum identity: {mlst_minid}%")
    else:
        log.info(f"  MLST minimum identity: default (95%)")

    if mlst_mincov is not None:
        cmd += ["--mincov", str(mlst_mincov)]
        log.info(f"  MLST minimum coverage: {mlst_mincov}%")
    else:
        log.info(f"  MLST minimum coverage: default (50%)")

    cmd += ["--threads", str(num_threads)]
    log.info(f"  MLST threads: {num_threads}")

    cmd += seq_files

    log.info(f"  Running mlst on {len(seq_files)} sequences...")
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env=os.environ.copy()
    )
    if result.returncode != 0:
        log.error(f"  mlst failed (exit {result.returncode})")
    if result.stderr.strip():
        log.debug(f"  mlst stderr: {result.stderr.strip()}")
    if result.returncode != 0:
        return

    log.info(result.stdout.strip())

    rows = [line.split('\t') for line in result.stdout.splitlines() if line.strip()]
    if rows:
        # Parse gene names from fields like "adk(1)", "recA(-)", "tpi(~3)"
        gene_names = [f.split('(')[0] for f in rows[0][3:]]
        # Preferred column order for the X. citri pv. citri scheme; falls back to
        # mlst's own gene order for any other scheme (no genes match the preference list).
        gene_order_preference = ["nusG", "rpmA", "arfB", "rplQ", "queD", "fdx", "exbD"]
        ordered_genes = ([g for g in gene_order_preference if g in gene_names]
                         + [g for g in gene_names if g not in gene_order_preference])

        # Write TSV with gene columns reordered (path/scheme/ST columns left as mlst gave them)
        with open(output_tsv, 'w') as f:
            for row in rows:
                fields_by_gene = {field.split('(')[0]: field for field in row[3:]}
                f.write('\t'.join(row[:3] + [fields_by_gene[g] for g in ordered_genes]) + '\n')

        # Write CSV with header row and basename-only filenames
        output_csv = os.path.join(output, "mlst_summary_results.csv")
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['genome', 'scheme', 'ST'] + ordered_genes)
            for row in rows:
                genome = os.path.basename(row[0])
                scheme = row[1]
                st = row[2]
                alleles_by_gene = dict(zip(gene_names, [f.split('(')[1].rstrip(')') for f in row[3:]]))
                writer.writerow([genome, scheme, st] + [alleles_by_gene[g] for g in ordered_genes])
        log.info(f"  CSV saved: {output_csv}")
    else:
        with open(output_tsv, 'w') as f:
            f.write(result.stdout)

    log.info(f"  Results saved: {output_tsv}")


# =============================================================================
# Orchestration
# =============================================================================

def pipeline(args):
    """Orchestrate all pipeline stages."""
    output = os.path.abspath(args.output)
    os.makedirs(output, exist_ok=True)

    extracted_dir = os.path.join(output, "1_extracted_targets")
    consolidated_dir = os.path.join(output, "2_consolidated_targets")
    alignments_dir = os.path.join(output, "3_alignments")
    supermatrix_dir = os.path.join(output, "4_supermatrix")
    tree_dir = os.path.join(output, "5_Final_MLSA_Tree")
    tree_dir_nogaps = os.path.join(output, "5_Final_MLSA_Tree_nogaps")

    mlst_dir = os.path.join(output, "6_MLST")

    if not args.mlst_only:
        dirs_to_create = [extracted_dir, consolidated_dir, alignments_dir, supermatrix_dir, tree_dir, tree_dir_nogaps]
        for d in dirs_to_create:
            if os.path.exists(d) and not args.force:
                sys.exit(f"Directory exists: {d}. Use --force to overwrite.")
            if os.path.exists(d) and args.force:
                shutil.rmtree(d)
            os.makedirs(d)

    targets_for_blast = os.path.abspath(args.target_typing_genes) if args.target_typing_genes else None

    if not args.mlst_only and not args.skip_extraction:
        log.info("\n[STAGE 1] Extracting target genes via BLAST...")
        blast_extraction(os.path.abspath(args.input_query_sequences), targets_for_blast,
                         extracted_dir, args.type_genes_extension, args.min_identity, args.threads)

    if not args.mlst_only and not args.skip_consolidation:
        log.info("\n[STAGE 2] Consolidating extracted targets...")
        consolidate_targets(extracted_dir, consolidated_dir, targets_for_blast,
                            args.type_genes_extension)

    if not args.mlst_only and not args.skip_alignment:
        log.info("\n[STAGE 3] Aligning genes with MAFFT...")
        align_genes(consolidated_dir, alignments_dir)
        if args.no_trim_amplicon_input:
            log.info("\n[STAGE 4] Skipping trimal (--no-trim-amplicon-input set) — using raw MAFFT alignments")
        else:
            log.info("\n[STAGE 4] Trimming alignments with trimal...")
            trim_alignments(alignments_dir)
        log.info("\n[STAGE 5] Creating supermatrices (with gaps and without gaps)...")
        supermatrix_fasta = create_supermatrix(alignments_dir, supermatrix_dir, args.min_genes, nogaps=False)
        supermatrix_fasta_nogaps = create_supermatrix(alignments_dir, supermatrix_dir, args.min_genes, nogaps=True)

    if not args.mlst_only and not args.skip_tree:
        supermatrix_fasta = os.path.join(supermatrix_dir, "supermatrix.fasta")
        supermatrix_fasta_nogaps = os.path.join(supermatrix_dir, "supermatrix_nogaps.fasta")
        if not os.path.exists(supermatrix_fasta):
            sys.exit(f"Supermatrix not found: {supermatrix_fasta}. "
                     f"Run without --skip-tree or run alignment first.")
        log.info("\n[STAGE 6a] Running RAxML on gap-preserved alignment and MinVar rooting...")
        run_raxml(supermatrix_fasta, supermatrix_dir, tree_dir, num_threads=args.threads, bootstraps=args.bootstraps, nogaps=False,
                  no_minvar_root=args.no_minvar_root, outgroup=args.outgroup)
        log.info("\n[STAGE 6b] Running RAxML on gap-removed alignment and MinVar rooting...")
        if os.path.exists(supermatrix_fasta_nogaps):
            run_raxml(supermatrix_fasta_nogaps, supermatrix_dir, tree_dir_nogaps, num_threads=args.threads, bootstraps=args.bootstraps, nogaps=True,
                      no_minvar_root=args.no_minvar_root, outgroup=args.outgroup)
        else:
            log.warning("  Gap-removed supermatrix not found, skipping RAxML for nogaps version")

    if args.run_mlst or args.mlst_only:
        os.makedirs(mlst_dir, exist_ok=True)
        log.info("\n[STAGE 7] Running MLST typing (tseemann/mlst)...")

        mlst_input_files = None
        if not args.mlst_only and os.path.exists(extracted_dir):
            # Build clean per-sample FASTAs from Stage 1 extracted targets.
            # This ensures tseemann/mlst receives one sequence per gene per sample
            # (best BLAST hit already selected), avoiding multi-hit confusion from
            # amplicon inputs.
            log.info("  Building per-sample FASTAs from Stage 1 extracted targets...")
            query_seq_files = []
            for ext in args.query_seq_extension:
                query_seq_files.extend(
                    glob.glob(os.path.join(os.path.abspath(args.input_query_sequences), f"*.{ext}"))
                )
            mlst_input_dir = os.path.join(mlst_dir, "mlst_input_fastas")
            mlst_input_files = build_per_sample_fastas_from_extracted(
                extracted_dir, sorted(set(query_seq_files)), mlst_input_dir
            )

        run_mlst_typing(os.path.abspath(args.input_query_sequences), mlst_dir,
                        args.query_seq_extension, args.custom_mlst_db,
                        args.mlst_scheme_name, input_files=mlst_input_files,
                        mlst_minid=args.mlst_minid, mlst_mincov=args.mlst_mincov, num_threads=args.threads)


def main():
    """Main entry point."""
    args = parse_arguments()
    os.makedirs(args.output, exist_ok=True)
    setup_logging(os.path.join(args.output, "pipeline.log"), args.verbose)

    log.info("=" * 80)
    log.info("mls_khajiit.py — Unified MLST/MLSA Phylogenetic Pipeline")
    log.info("=" * 80)

    # Validate inputs
    if not os.path.isdir(args.input_query_sequences):
        sys.exit(f"Error: --input-query-sequences directory not found: {args.input_query_sequences}")

    if not args.mlst_only:
        if not args.target_typing_genes:
            sys.exit("Error: --target-typing-genes is required unless --mlst-only is set")
        if not os.path.isdir(args.target_typing_genes):
            sys.exit(f"Error: --target-typing-genes directory not found: {args.target_typing_genes}")

    if (args.run_mlst or args.mlst_only) and not args.mlst_scheme_name:
        sys.exit("Error: --mlst-scheme-name is required when --run-mlst or --mlst-only is set")

    if args.custom_mlst_db:
        if not os.path.isdir(args.custom_mlst_db):
            sys.exit(f"Error: --custom-mlst-db directory not found: {args.custom_mlst_db}")
        if args.mlst_scheme_name:
            scheme_normalized = re.sub(r'[^a-z0-9]', '', args.mlst_scheme_name.lower())
            scheme_dir = os.path.join(args.custom_mlst_db, scheme_normalized)
            if not os.path.isdir(scheme_dir):
                sys.exit(f"Error: scheme directory not found in custom database: {scheme_dir}")
        blast_db = os.path.join(args.custom_mlst_db, "blast", "mlst.njs")
        if not os.path.exists(blast_db):
            sys.exit(f"Error: BLAST database not found in custom database: "
                     f"{os.path.join(args.custom_mlst_db, 'blast', 'mlst')}")

    log.debug(f"Input query sequences: {os.path.abspath(args.input_query_sequences)}")
    log.debug(f"Targets:               {os.path.abspath(args.target_typing_genes) if args.target_typing_genes else 'N/A (mlst-only)'}")
    log.debug(f"Output:                {os.path.abspath(args.output)}")
    log.debug(f"Query seq extension:   {args.query_seq_extension}")
    log.debug(f"Gene extensions:       {args.type_genes_extension}")
    log.debug(f"Min genes:             {args.min_genes}")
    log.debug(f"Multi-allele input:    {args.multi_allele}")
    log.debug(f"BLAST min-identity:    {args.min_identity}%")
    log.debug(f"MLST minid:            {args.mlst_minid if args.mlst_minid is not None else f'{args.min_identity}% (from --min-identity)'}")
    log.debug(f"MLST mincov:           {args.mlst_mincov if args.mlst_mincov is not None else 'default (50%)'}")
    log.debug(f"Run mlst:              {args.run_mlst}")

    pipeline(args)

    log.info("\n" + "=" * 80)
    log.info("Pipeline complete.")
    if pipeline_warnings:
        log.info("\nWARNINGS:")
        for w in pipeline_warnings:
            log.info(f"  ! {w}")
    log.info("=" * 80)


if __name__ == "__main__":
    main()
