# QIIME 2 + PICRUSt2 Pipeline for `proj_SRP145069`

This README explains how to use `qiime2_script_proj_SRP145069_relative_with_manifest.sh`, what input files are required, how the folders should be arranged, and what output files are produced.

The workflow processes paired-end 16S reads with QIIME 2, classifies ASVs with a SILVA V4 classifier, exports the feature table and taxonomy, and then runs PICRUSt2 custom trait prediction using `./concatenated_feature_table.tsv` as the PICRUSt2 feature-table input.

---

## 1. Files you must place in the project directory

Put all required files in the same working directory as the shell script.

Expected project directory:

```text
proj_SRP145069/
├── qiime2_script_proj_SRP145069_relative_with_manifest.sh
├── make_manifest.py
├── silva138_AB_V4_classifier.qza
├── feature_table.csv
├── concatenated_feature_table.tsv
└── reads/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
```

The `reads/` folder must be in the same directory as the script. Put the raw paired-end FASTQ files inside `reads/`.

---

## 2. Required input files

### `./reads/`

This folder contains the raw paired-end FASTQ files.

Accepted read-name styles include:

```text
sample1_R1.fastq.gz / sample1_R2.fastq.gz
sample1_1.fastq.gz  / sample1_2.fastq.gz
sample1.forward.fastq.gz / sample1.reverse.fastq.gz
```

The helper script `make_manifest.py` detects forward and reverse reads from the file names and writes a QIIME 2 paired-end manifest.

---

### `./make_manifest.py`

This script is included with the pipeline.

It is called automatically by the shell script:

```bash
python ./make_manifest.py ./reads ./manifest_mmdb1.csv
```

It creates:

```text
manifest_mmdb1.csv
```

The manifest has three columns:

```text
sample-id,absolute-filepath,direction
```

Although QIIME 2 uses the column name `absolute-filepath`, this version writes portable relative paths such as `reads/sample1_R1.fastq.gz`.

---

### `./silva138_AB_V4_classifier.qza`

Download or prepare the SILVA V4 QIIME 2 classifier and place it in the same directory as the script.

The file must be named exactly:

```text
silva138_AB_V4_classifier.qza
```

The taxonomy classification step uses:

```bash
qiime feature-classifier classify-sklearn   --i-classifier ./silva138_AB_V4_classifier.qza   --i-reads rep-seqs.qza   --o-classification taxonomy.qza
```

Use a classifier compatible with your QIIME 2 version and the amplified 16S region. This script expects a V4 classifier.

---

### `./concatenated_feature_table.tsv`

This is the feature table used by PICRUSt2.

The script runs PICRUSt2 using:

```bash
picrust2_pipeline.py   -s exported-rep-seqs/dna-sequences.fasta   -i ./concatenated_feature_table.tsv   -o picrust2_output   --processes 15   --custom_trait_tables ./feature_table.csv   --no_pathways
```

Important: the feature IDs in `concatenated_feature_table.tsv` must match the sequence IDs in:

```text
exported-rep-seqs/dna-sequences.fasta
```

If the IDs do not match, PICRUSt2 may fail or discard features.

---

### `./feature_table.csv`

This is the custom PICRUSt2 trait table.

Place it in the same directory as the script and name it exactly:

```text
feature_table.csv
```

The script passes it to PICRUSt2 using:

```bash
--custom_trait_tables ./feature_table.csv
```

---

## 3. Required software

You need QIIME 2, BIOM, conda, and PICRUSt2.

Before running the shell script, activate your QIIME 2 environment. Example:

```bash
conda activate qiime2-2024.2
```

Use your actual QIIME 2 environment name if different.

The script later switches into the PICRUSt2 environment using:

```bash
conda activate picrust2
```

Therefore, you must already have a conda environment named `picrust2` with PICRUSt2 installed.

---

## 4. How to run the pipeline

From inside the project directory:

```bash
bash qiime2_script_proj_SRP145069_relative_with_manifest.sh
```

You can also run the script from another directory because it begins with:

```bash
cd "$(dirname "$0")"
```

That forces all relative paths to resolve from the directory containing the script.

---

## 5. What the pipeline does

The shell script performs these steps:

1. Changes into the script directory.
2. Generates `manifest_mmdb1.csv` from `./reads/` using `make_manifest.py`.
3. Imports paired-end reads into QIIME 2 as `single-end-demux.qza`.
4. Runs `qiime cutadapt trim-paired`.
5. Creates the read summary visualization `demux.qzv`.
6. Runs DADA2 denoising.
7. Produces the ASV table `table.qza`.
8. Produces representative sequences `rep-seqs.qza`.
9. Classifies representative sequences using the local SILVA classifier.
10. Exports taxonomy, representative sequences, and the feature table.
11. Converts the QIIME 2 BIOM table into `feature-table.tsv`.
12. Runs PICRUSt2 custom trait prediction using `./concatenated_feature_table.tsv`.
13. Copies the main outputs into `../all_output`.

---

## 6. Main intermediate outputs

These files are created inside the project directory.

### `manifest_mmdb1.csv`

QIIME 2 paired-end manifest generated from the `reads/` folder.

---

### `single-end-demux.qza`

Imported QIIME 2 sequence artifact.

Despite the filename, the import type is paired-end data:

```bash
SampleData[PairedEndSequencesWithQuality]
```

---

### `single-end-trim-demux.qza`

Trimmed paired-end reads produced by cutadapt.

---

### `demux.qzv`

QIIME 2 visualization summarizing the trimmed reads.

View it with:

```bash
qiime tools view demux.qzv
```

---

### `table.qza`

DADA2 denoised ASV feature table.

---

### `rep-seqs.qza`

Representative ASV sequences from DADA2.

---

### `denoising-stats.qza`

DADA2 denoising statistics.

---

### `taxonomy.qza`

Taxonomic classification of representative sequences using the SILVA classifier.

---

### `exported-taxonomy/taxonomy.tsv`

Exported taxonomy table.

---

### `exported-feature-table/feature-table.biom`

Exported BIOM feature table from QIIME 2.

---

### `feature-table.tsv`

TSV version of the exported BIOM table.

This is copied to the final output folder as:

```text
../all_output/proj_SRP145069_feature-table.tsv
```

---

### `exported-rep-seqs/dna-sequences.fasta`

Representative sequences exported as FASTA.

PICRUSt2 uses this file as the sequence input.

---

### `picrust2_output/`

PICRUSt2 output directory.

Because the script uses `--no_pathways`, it produces trait prediction outputs but does not infer pathway abundances.

---

## 7. Final outputs copied to `../all_output`

The script creates the output folder:

```text
../all_output
```

Then it copies or decompresses the main output files there.

### `../all_output/proj_SRP145069_pred_metagenome_unstrat.tsv`

Unstratified predicted metagenome table from PICRUSt2.

This is generated from:

```text
./picrust2_output/feature_table_metagenome_out/pred_metagenome_unstrat.tsv.gz
```

This is usually the main custom-trait prediction output.

---

### `../all_output/proj_SRP145069_feature_table_predicted.tsv`

Predicted feature table from PICRUSt2.

This is generated from:

```text
./picrust2_output/feature_table_predicted.tsv.gz
```

---

### `../all_output/proj_SRP145069_taxonomy.tsv`

Exported QIIME 2 taxonomy table.

This is copied from:

```text
./exported-taxonomy/taxonomy.tsv
```

---

### `../all_output/proj_SRP145069_feature-table.tsv`

DADA2 feature table exported from QIIME 2 and converted from BIOM to TSV.

This is copied from:

```text
./feature-table.tsv
```

---

## 8. Common checks before running

Before running the script, confirm that these paths exist:

```text
./qiime2_script_proj_SRP145069_relative_with_manifest.sh
./make_manifest.py
./reads/
./silva138_AB_V4_classifier.qza
./concatenated_feature_table.tsv
./feature_table.csv
```

Also check that your paired-end FASTQ files have recognizable R1/R2 names.

Good examples:

```text
SRR123456_R1.fastq.gz
SRR123456_R2.fastq.gz
```

Problematic examples:

```text
SRR123456.fastq.gz
SRR123456_readA.fastq.gz
SRR123456_readB.fastq.gz
```

If the names do not clearly identify forward and reverse reads, rename them before running `make_manifest.py`.

---

## 9. Troubleshooting

### `ERROR: reads directory not found`

Create a folder named `reads` in the same directory as the script and place the FASTQ files inside it.

---

### `ERROR: no complete paired-end samples were found`

The read files were not recognized as R1/R2 pairs.

Rename them using a clear paired-end pattern such as:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

---

### SILVA classifier not found

Make sure the classifier is in the same directory as the script and named exactly:

```text
silva138_AB_V4_classifier.qza
```

---

### PICRUSt2 drops many features

Check that the IDs in `concatenated_feature_table.tsv` match the sequence IDs in `exported-rep-seqs/dna-sequences.fasta`.

---

### Conda activation fails inside the script

Run this once before executing the pipeline:

```bash
conda init bash
```

Then close and reopen the terminal, activate QIIME 2, and rerun the script.
