FASTA_GZ=$1 
GTF_GZ=$2
kb ref --workflow=standard -i index.idx -g t2g.txt -f1 cdna.fa \
  --include-attribute gene_biotype:protein_coding \
  --include-attribute gene_biotype:lncRNA \
  --include-attribute gene_biotype:lincRNA \
  --include-attribute gene_biotype:antisense \
  --include-attribute gene_biotype:IG_LV_gene \
  --include-attribute gene_biotype:IG_V_gene \
  --include-attribute gene_biotype:IG_V_pseudogene \
  --include-attribute gene_biotype:IG_D_gene \
  --include-attribute gene_biotype:IG_J_gene \
  --include-attribute gene_biotype:IG_J_pseudogene \
  --include-attribute gene_biotype:IG_C_gene \
  --include-attribute gene_biotype:IG_C_pseudogene \
  --include-attribute gene_biotype:TR_V_gene \
  --include-attribute gene_biotype:TR_V_pseudogene \
  --include-attribute gene_biotype:TR_D_gene \
  --include-attribute gene_biotype:TR_J_gene \
  --include-attribute gene_biotype:TR_J_pseudogene \
  --include-attribute gene_biotype:TR_C_gene \
  $FASTA_GZ $GTF_GZ
