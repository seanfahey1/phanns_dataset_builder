#!/usr/bin/env bash

if [ -n "$1" ]
then
  p3-get-feature-sequence -i "$1" > OTH.fasta
else
  echo "pass in a file with a BRC-ID on each line"
#  https://www.bv-brc.org/view/GenomeList/?eq(genome_name,phage)#view_tab=proteins&filter=and(or(eq(feature_type,CDS),eq(feature_type,mat_peptide)),eq(annotation,PATRIC),between(aa_length,50,1000000))
fi
