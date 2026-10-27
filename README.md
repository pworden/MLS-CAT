# MLS-Khajiit

Unified MLSA (Multilocus Sequence Analysis) and MLST (Multilocus Sequence Typing) phylogenetic pipeline. Can be run as a combined MLSA + MLST workflow, or as MLST-only.

MLSA pipeline based on MAGe by Jwebster89 (https://github.com/Jwebster89/mage).
MLST typing powered by tseemann/mlst (https://github.com/tseemann/mlst).

---

## What it does

### MLSA pipeline (Stages 1–6)

| Stage | Description |
|-------|-------------|
| 1 | **BLAST extraction** — identifies and extracts target typing gene sequences from each input genome/amplicon |
| 2 | **Consolidation** — merges per-genome BLAST hits into per-gene FASTA files |
| 3 | **Alignment** — aligns each gene with MAFFT |
| 4 | **Trimming** — trims alignments with trimAl |
| 5 | **Supermatrix** — concatenates trimmed alignments into a single supermatrix FASTA |
| 6 | **Tree** — infers a maximum-likelihood tree with RAxML, then roots it (MinVar by default, or on a specified outgroup, or left unrooted) |

### MLST typing (Stage 7)

Runs tseemann/mlst on the input genomes/amplicons to assign Sequence Types (STs). Supports both the standard pubmlst database and custom schemes added directly to tseemann/mlst (see "Custom MLST schemes" below).

---

## Installation

```bash
git clone https://github.com/pworden/MLS-Khajiit.git
cd MLS-Khajiit
conda env create -f environment.yml
conda activate MLS_Khajiit
```

This sets up everything needed — `blastn`, `mafft`, `trimal`, `raxml`, `mlst`, and the `mls-khajiit` command itself.

Already have an environment with those tools installed? Just activate it and run `pip install -e .` instead.

---

## Requirements

- Conda environment `MLS_Khajiit` (or any environment, any name) containing: `blastn`, `mafft`, `trimal`, `raxml`, `biopython`, `mlst` (all included by `environment.yml`)

---

## Input

| Input | Description |
|-------|-------------|
| Query sequences | Folder of genome FASTA files or amplicon consensus sequences (one file per sample) |
| Target typing genes | Folder of reference gene FASTA files used as BLAST queries (one file per gene). For multi-allele MLST schemes, each file contains all known alleles for that gene. Not required for `--mlst-only`. |

---

## Output

All output is written to the directory specified by `--output`, organised into numbered subdirectories:

```
output/
├── 1_extracted_targets/      BLAST output files (one per gene × genome combination)
├── 2_consolidated_targets/   Per-gene FASTA files with sequences from all genomes
├── 3_alignments/             MAFFT alignments and trimAl-trimmed versions
├── 4_supermatrix/            Concatenated supermatrix FASTA
├── 5_Final_MLSA_Tree/        RAxML tree files and rooted final tree (MinVar- or outgroup-rooted, depending on flags)
├── 6_MLST/
│   ├── mlst_results.tsv      Raw tseemann/mlst output (full file paths)
│   └── mlst_summary_results.csv  Formatted summary with header row and filenames only
└── pipeline.log              Full run log
```

### mlst_summary_results.csv columns

| Column | Description |
|--------|-------------|
| `genome` | Input filename (no path) |
| `scheme` | MLST scheme used |
| `ST` | Assigned Sequence Type (`-` if unassigned) |
| `adk`, `atpA`, ... | Allele number for each scheme gene (`-` absent, `~N` novel-ish, `N?` partial) |

---

## Usage

### Full MLSA + MLST pipeline (multi-allele scheme with custom database)

```bash
conda activate MLS_Khajiit

mls-khajiit \
    --input-query-sequences /path/to/genomes \
    --target-typing-genes /path/to/mlst_genes \
    --output /path/to/output \
    --query-seq-extension fna \
    --type-genes-extension fas \
    --min-genes 3 \
    --multi-allele \
    --run-mlst \
    --custom-mlst-db /path/to/custom_mlst_db \
    --mlst-scheme-name yourspecies \
    --verbose --force
```

### MLSA only (single representative sequence per gene, no MLST)

```bash
mls-khajiit \
    --input-query-sequences /path/to/amplicons \
    --target-typing-genes /path/to/mlst_genes \
    --output /path/to/output \
    --query-seq-extension fasta \
    --type-genes-extension fna \
    --min-genes 3 \
    --verbose --force
```

### MLST only (no tree building)

```bash
mls-khajiit \
    --input-query-sequences /path/to/genomes \
    --output /path/to/output \
    --query-seq-extension fna \
    --mlst-only \
    --custom-mlst-db /path/to/custom_mlst_db \
    --mlst-scheme-name yourspecies \
    --verbose
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `--input-query-sequences` | Folder of input genome or amplicon FASTA files |
| `--target-typing-genes` | Folder of target gene reference FASTA files (required unless `--mlst-only`) |
| `--output` | Output directory |

### Optional — General

| Argument | Default | Description |
|----------|---------|-------------|
| `-e`, `--query-seq-extension` | `fna fasta` | File extension(s) of input query sequences |
| `--type-genes-extension` | `fasta fa fas fna tfa` | File extension(s) of target gene files |
| `--multi-allele` | off | Target gene files contain multiple alleles per gene (MLST scheme format). BLAST keeps only the single best hit per gene per genome. |
| `--min-genes` | `2` | Minimum number of genes a genome must have to be included in the tree |
| `--min-identity` | `85.0` | Minimum percent identity for BLAST hits |
| `--no-trim-amplicon-input` | off | Skip trimAl trimming and use raw MAFFT alignments for the supermatrix. Recommended when input sequences are amplicons (already clean, uniform length). |
| `--bootstraps` | `100` | Number of RAxML bootstrap replicates |
| `--threads` | `20` | Number of threads for BLAST, MLST, and RAxML |
| `-f`, `--force` | off | Overwrite existing output directories |
| `-v`, `--verbose` | off | Verbose logging output |

### Optional — MLST

| Argument | Default | Description |
|----------|---------|-------------|
| `--run-mlst` | off | Run MLST typing (Stage 7) in addition to MLSA |
| `--mlst-only` | off | Run MLST typing only — skips all MLSA stages. `--target-typing-genes` not required. Implies `--run-mlst`. |
| `--custom-mlst-db` | None | Path to a custom MLST database directory (see "Custom MLST schemes" below). Use for species not in the standard pubmlst database. |
| `--mlst-scheme-name` | None | Scheme name to pass to tseemann/mlst (e.g. `yourspecies`). Automatically normalised to tseemann/mlst convention (lowercase, no underscores). Required when using `--custom-mlst-db`. |
| `--mlst-minid` | None | MLST minimum allele identity threshold. If unset, uses the `--min-identity` value. |
| `--mlst-mincov` | None | MLST minimum coverage for partial alleles. If unset, uses tseemann/mlst's own default (50). |

### Optional — Tree rooting

| Argument | Default | Description |
|----------|---------|-------------|
| `--no-minvar-root` | off | Skip FastRoot MinVar rooting of the final tree (produces only the plain, unrooted `.final` tree). |
| `--outgroup` | None | Path to the FASTA file (as named in `--input-query-sequences`) of the sample to root on, instead of MinVar rooting. Implies `--no-minvar-root` — MinVar never runs if `--outgroup` is set. Writes `MAGE_tree_output_OutgroupRooted.final`. Errors and skips (no file written) if the derived taxon name matches zero or more than one tip in the tree. |

### Optional — Skip stages

| Argument | Description |
|----------|-------------|
| `--skip-extraction` | Skip Stage 1 (BLAST extraction) |
| `--skip-consolidation` | Skip Stage 2 (consolidation) |
| `--skip-alignment` | Skip Stages 3–5 (alignment, trimming, supermatrix) |
| `--skip-tree` | Skip Stage 6 (RAxML tree) |

---

## Custom MLST schemes

For species not in the standard pubmlst database, add the scheme directly to tseemann/mlst following its own instructions (https://github.com/tseemann/mlst#adding-a-scheme).

### Adding a scheme

> **Note:** A custom MLST scheme can be created in a format compatible with the MLST app by following the instructions at this GitHub repo: https://github.com/tseemann/mlst

1. Create a folder named for the scheme under `mlst/db/pubmlst/`, e.g. `mlst/db/pubmlst/yourspecies/`.
2. Put two kinds of file in it:
   - `yourspecies.txt` — the profiles file: tab-separated, one ST per row, with a header row containing an `ST` column plus one column per locus. (`mlst_clade`, `clonal_complex`, `species`, `CC`, and `Lineage` columns are ignored if present.)
   - `LOCUS.tfa` for each locus named in the profiles file — nucleotide FASTA, one file per locus, with sequence IDs as `>allele_number` (avoid ambiguous IUPAC bases where possible).
3. Run `mlst/scripts/mlst-make_blast_db` to rebuild the BLAST index.
4. Confirm it's registered: `mlst --info | grep yourspecies`.
5. Test it: `mlst --scheme yourspecies some_genome.fasta`.

Once added this way, use it with `mls-khajiit` via `--mlst-scheme-name yourspecies` — no `--custom-mlst-db` needed, since the scheme is now part of mlst's own database.

`--custom-mlst-db` (see Arguments below) is for pointing at a separate, self-built database elsewhere instead of touching the bundled mlst install — it expects `DIR/yourspecies/yourspecies.txt` + `.tfa` files, and `DIR/blast/mlst*` (a `makeblastdb` nucleotide database built from those alleles, with headers prefixed `yourspecies.`). Building that database isn't covered by this tool; the steps above are the simpler path for most cases.

---

## Notes

- By default the final tree is MinVar-rooted and saved as `MAGE_tree_output_MinVarRooted.final`, uploadable directly to iTOL for visualisation. Pass `--outgroup FASTA_PATH` to root on a specific taxon instead (writes `MAGE_tree_output_OutgroupRooted.final`), or `--no-minvar-root` to skip rooting entirely and root manually (e.g. in iTOL).
- Genomes with fewer than `--min-genes` extracted genes are excluded from the supermatrix and tree but are still typed by MLST if `--run-mlst` or `--mlst-only` is set.
- All pipeline progress and tool output is recorded in `pipeline.log`.
