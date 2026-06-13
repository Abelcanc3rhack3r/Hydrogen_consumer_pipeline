cd "$(dirname "$0")"

python ./make_manifest.py ./reads ./manifest_mmdb1.csv 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manifest_mmdb1.csv \
  --output-path single-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences single-end-demux.qza \
  --p-cores 15 \
  --o-trimmed-sequences single-end-trim-demux.qza \
  #--p-adapter-f GGACTACHVGGGTWTCTAAT    \
  #--p-front-f GTGCCAGCMGCCGCGGTAA \
  #--p-adapter-r TTACCGCGGCKGCTGGCAC   \
  #--p-front-r ATTAGAWACCCBDGTAGTCC
qiime demux summarize \
  --i-data single-end-trim-demux.qza\
  --o-visualization demux.qzv





  #denoise single reads using DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-trim-demux.qza \
  --p-trim-left 13 \
  --p-trunc-len 240 \
  --p-n-threads 8 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ./silva138_AB_V4_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

  qiime tools export \
     --input-path taxonomy.qza \
     --output-path exported-taxonomy
  qiime tools export \
     --input-path table.qza \
      --output-path exported-feature-table

  biom convert -i exported-feature-table/feature-table.biom -o feature-table.tsv --to-tsv


  qiime tools export \
     --input-path rep-seqs.qza \
     --output-path exported-rep-seqs
  conda activate picrust2

  picrust2_pipeline.py \
    -s exported-rep-seqs/dna-sequences.fasta \
    -i ./concatenated_feature_table.tsv \
    -o picrust2_output \
    --processes 15 \
 --custom_trait_tables ./feature_table.csv \
    --no_pathways

   mkdir -p ../all_output

  gunzip -c ./picrust2_output/feature_table_metagenome_out/pred_metagenome_unstrat.tsv.gz > ../all_output/proj_SRP145069_pred_metagenome_unstrat.tsv
    gunzip -c ./picrust2_output/feature_table_predicted.tsv.gz > ../all_output/proj_SRP145069_feature_table_predicted.tsv
    cp ./exported-taxonomy/taxonomy.tsv ../all_output/proj_SRP145069_taxonomy.tsv

       cp ./feature-table.tsv ../all_output/proj_SRP145069_feature-table.tsv