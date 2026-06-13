cd /oceanstor/home/e1103389/MMDB/proj_SRP145069

qiime vsearch uchime-denovo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --output-dir uchime-dn-out

qiime metadata tabulate \
  --m-input-file uchime-dn-out/stats.qza \
  --o-visualization uchime-dn-out/stats.qzv


qiime feature-table filter-features \
  --i-table table.qza \
  --m-metadata-file uchime-dn-out/nonchimeras.qza \
  --o-filtered-table table-nonchimeric.qza
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --m-metadata-file uchime-dn-out/nonchimeras.qza \
  --o-filtered-data rep-seqs-nonchimeric.qza

qiime feature-classifier classify-sklearn \
  --i-classifier /oceanstor/home/e1103389/MMDB/test/silva138_AB_V4_classifier.qza \
  --i-reads rep-seqs-nonchimeric.qza \
  --o-classification taxonomy.qza

  qiime tools export \
     --input-path taxonomy.qza \
     --output-path exported-taxonomy
  qiime tools export \
     --input-path table-nonchimeric.qza \
      --output-path exported-feature-table
  conda activate picrust2

  picrust2_pipeline.py \
    -s exported-rep-seqs/dna-sequences.fasta \
    -i exported-feature-table/feature-table.biom \
    -o picrust2_output \
    --processes 15 \
 --custom_trait_tables  /oceanstor/home/e1103389/JGI_picrust/feature_table2.tsv  \
    --no_pathways

   mkdir -p ../all_output

  gunzip -c ./picrust2_output/feature_table_metagenome_out/pred_metagenome_unstrat.tsv.gz > ../all_output/proj_SRP145069_pred_metagenome_unstrat.tsv
    gunzip -c ./picrust2_output/feature_table_predicted.tsv.gz > ../all_output/proj_SRP145069_feature_table_predicted.tsv
    cp ./exported-taxonomy/taxonomy.tsv ../all_output/proj_SRP145069_taxonomy.tsv

       cp ./feature-table.tsv ../all_output/proj_SRP145069_feature-table.tsv