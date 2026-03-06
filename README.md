# Synthetic Visium HD FASTQ Simulator

A Nextflow DSL2 pipeline for generating synthetic paired-end FASTQ files that mimic Visium HD spatial transcriptomics data. The simulator uses real Visium HD count data as input and produces realistic reads with 36 nt spatial barcodes, a two-stage barcode error model (array synthesis + sequencing), PCR duplication, and full ground-truth output for pipeline validation.

---

## Input data

The simulation requires the Visium HD mouse brain `feature_slice.h5` file from 10x Genomics. Download it with:

```bash
mkdir -p data
wget -P data/ https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain/Visium_HD_3prime_Mouse_Brain_feature_slice.h5
```

You also need a transcript FASTA for cDNA read generation (GENCODE vM33 mouse):

```bash
mkdir -p data/reference
wget -P data/reference/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz
gunzip data/reference/gencode.vM33.transcripts.fa.gz
```

Set `feature_h5` and `reference_fasta` in `nextflow.config` (or via `--feature_h5` / `--reference_fasta` on the command line) to point to these files.

---

## Quick start

```bash
# 1. Create and activate the Conda environment
conda env create -f environment.yml
conda activate stdata_sim

# 2. Run the simulation
nextflow run main.nf -resume
```

Default settings simulate **50 million reads** from **100,000 spots** and **1,000 genes**. Edit `nextflow.config` to change parameters.

---

## Pipeline overview

The simulation proceeds in four stages that mirror the real experimental workflow:

**Stage 1 - Spatial UMI counts, spot selection, and barcode synthesis errors** (`scripts/run_sccigen_sim.R`)

Subsamples G genes (stratified by expression level) and S spots from the Visium HD H5. Assigns each spot a random 36 nt barcode, then applies synthesis errors once per spot to model array printing defects. Both the perfect barcode (`true_barcode`) and the synthesis-mutated barcode are written to `spot_coordinates.tsv`.

**Stage 2 - Molecule sampling with PCR family sizes**

Converts the count matrix into a pool of molecules and samples until the target read depth is reached. Each molecule receives a family size drawn from a geometric distribution (PCR duplication model) and a unique UMI shared by all its PCR copies.

**Stage 3 - cDNA read generation** (`bin/generate_cdna_reads.py`)

Resolves each molecule to a transcript in the reference FASTA and extracts a deterministic read pair, ensuring PCR duplicates originate from the same fragment.

**Stage 4 - PCR expansion and sequencing errors** (`bin/expand_pcr_duplicates.py`)

Expands each molecule into f_m reads. Applies per-read sequencing errors to the barcode region of Read 1 only; Read 2 is emitted unchanged. Reads are assigned globally unique IDs of the form `read_s{shard}_{n}` across all parallel shards.

### Barcode error model

| Stage | Applied | Rates (default) | Represents |
|---|---|---|---|
| Synthesis | Once per spot | sub=0.03, ins=0.04, del=0.14 | Array printing defects |
| Sequencing | Per read (R1 barcode only) | sub=0.001, ins=0.0001, del=0.0001 | Illumina base-calling errors |

---

## Configuration

All parameters are in `nextflow.config`. Key settings:

| Parameter | Default | Description |
|---|---|---|
| `total_reads` | 50,000,000 | Target number of reads |
| `n_spots` | 100,000 | Number of simulated spots |
| `n_genes` | 1,000 | Genes to sample from H5 |
| `required_genes` | `"Bcl11b,Tbr1,..."` | Genes force-included regardless of sampling |
| `barcode_length` | 36 | Barcode length (nt) |
| `umi_length` | 10 | UMI length (nt) |
| `synth_sub_rate` | 0.03 | Synthesis substitution rate |
| `synth_ins_rate` | 0.04 | Synthesis insertion rate |
| `synth_del_rate` | 0.14 | Synthesis deletion rate |
| `seq_sub_rate` | 0.001 | Sequencing substitution rate |
| `seq_ins_rate` | 0.0001 | Sequencing insertion rate |
| `seq_del_rate` | 0.0001 | Sequencing deletion rate |
| `pcr_geom_p` | 0.5 | PCR geometric dist. parameter (mean = 1/p) |
| `pcr_max_family` | 20 | PCR family size cap |
| `sim_shards` | 8 | Parallelisation shards |
| `feature_h5` | null | Path to Visium HD `feature_slice.h5` |
| `reference_fasta` | -- | Path to transcript FASTA |
| `out_dir` | `results/` | Output directory |

---

## Outputs

| File | Description |
|---|---|
| `simulated_R1.fastq.gz` | Read 1 (barcode + UMI + cDNA suffix) |
| `simulated_R2.fastq.gz` | Read 2 (cDNA) |
| `spot_coordinates.tsv` | spot_id, barcode (synthesis-mutated), x, y, true_barcode |
| `read_to_spot.tsv` | read_id, spot_id, transcript_id, umi, molecule_id |
| `barcode_mutations.tsv` | Per-read: original_barcode, synthesis_barcode, mutated_barcode, mutation_events |
| `molecules.tsv` | One row per unique molecule with family size and fragment start |
| `expression_matrix.tsv` | Spot x gene UMI count matrix |
| `summary_metrics.json` | Read/spot counts summary |

The three barcode columns in `barcode_mutations.tsv` reflect the two-stage error model:

- `original_barcode` -- perfect designed sequence (`true_barcode` from `spot_coordinates.tsv`)
- `synthesis_barcode` -- after array synthesis errors (shared by all reads from a spot)
- `mutated_barcode` -- after per-read sequencing errors (unique per read)

---

## Utility scripts

**`scripts/read_to_spot_csv.py`** -- Convert `read_to_spot.tsv` to a deduplicated CSV with `barcode,x,y` (one row per unique spot/barcode):

```bash
python3 scripts/read_to_spot_csv.py \
    --read_to_spot  results/<run>/read_to_spot.tsv \
    --spot_coords   results/<run>/spot_coordinates.tsv \
    --out           results/<run>/read_to_spot.csv
```

---

## Requirements

- Nextflow 21.04+
- Python 3.8+
- R 4.2+ with `rhdf5`, `Matrix`
- Conda (recommended) -- see `environment.yml`

### R dependencies

```r
install.packages("Matrix")
BiocManager::install("rhdf5")
```

---

## Repository structure

```
.
+-- main.nf                   # Main Nextflow workflow
+-- nextflow.config           # Parameters and profiles
+-- environment.yml           # Conda environment definition
+-- Dockerfile                # Container definition
+-- modules/                  # Nextflow process modules
|   +-- run_srt_sim.nf
|   +-- expand_pcr_duplicates.nf
|   +-- ...
+-- bin/                      # Python helper scripts
|   +-- expand_pcr_duplicates.py
|   +-- attach_spatial_barcodes.py
|   +-- generate_cdna_reads.py
|   +-- ...
+-- scripts/
|   +-- run_sccigen_sim.R     # Stage 1: spatial simulation + synthesis errors
|   +-- read_to_spot_csv.py   # Utility: deduplicate read_to_spot to barcode CSV
+-- params_example.yaml       # Example parameter file
+-- manuscript_methods_simulation.tex  # Methods section LaTeX source
```

---

## Citation

If you use this simulator, please cite the associated manuscript (in preparation).
