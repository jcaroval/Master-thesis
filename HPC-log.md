####Content:
• 45 Illumina Short-Read (PE150) Samples (Acrobeles emmatus)
• 1 Reference Genome (Acrobeloides tricornis)
• 1 UCE file (based on 7 different Panagrolaimus sp.)

####Tools used:
• 

**1. Data Curation and Quality**

1.1 Download and unzip the provided files
```
wget URL
tar -xvf filename.tar
```

1.2 Calculate and compare MD5sums. They should resemble the MD5sums that were provided by Novogene.
```
find . -type f -name '*.fq.gz' | while read -r file; do
echo "current file: $file"
md5sum "$file" >> md5sums.txt
done
```
1.3 Create MultiQC files to investigate the quality of the raw files.
```
multiqc *.fq.gz
```
1.4 Trim the sequences. In this case, 1 is always the forward, 2 the reverse file.
```
fastp -i ID_1.fq.gz -I ID_2.fq.gz -o ID_trimmed_1.fq.gz -O ID_trimmed_2.fq.gz
```
or for multiple files:
```
for dir in EPT_*; do
  (cd "$dir" && while read entry; do
    fastp -i $entry"_1.fq.gz" -I $entry"_2.fq.gz" -o ./$entry"_trimmed_1.fq.gz" -O ./$entry"_trimmed_2.fq.gz"; done < list_files);
done
```
1.5 Repeat step 1.3 (MultiQC) on the trimmed files.

**2. Reference-based Genome Alignments (cross-species mapping against of Acrobeles emmatus against Acrobeloides tricornis)**

2.1 Mapping with bwa-mem2

2.1.1 Create an index folder of the reference genome
```
bwa-mem2 index PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa
```
2.1.2 Map the sequences against the reference genome
```
bwa-mem2 mem -t 32 PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa ID_trimmed_1.fq.gz ID_trimmed_2.fq.gz -o ID_out.sam
```
2.1.3 Convert sam files to sorted bam files
```
samtools view -S -b ID_out.sam > ID_out.bam
samtools sort ID_out.bam > ID_sorted.bam
```
2.1.4 Check the mapping accuracy with samtools
```
samtools flagstat ID_sorted.bam
```
2.1.5 Check the mapping accuracy with Qualimap
```
qualimap bamqc -bam ID_sorted.bam -gff busco_nematodes_all.gff
````
2.2 Mapping with nextgenmap

2.2.1 Map the sequences against the reference genome
```
ngm -r PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa -1 ID_trimmed_1.fq.gz -2 ID_trimmed_2.fq.gz -o ID_out.sam -t 64
```
2.2.2 Check the mapping accuracy with samtools and Qualimap similar to steps 2.1.4 and 2.1.5

**3. De-novo Genome Assembly**

3.1 Concatenate the sequences of all 45 samples.
```
cat *_trimmed_1.fq.gz >> EPT_trimmed_1.fq.gz
cat *_trimmed_2.fq.gz >> EPT_trimmed_2.fq.gz
```
3.2 Downsample the sequences to a desired fraction. This is done to reduce the runtime of used tools.
```
seqtk EPT_trimmed_1.fq.gz 0.2 > sub20_1.fq.gz
seqtk EPT_trimmed_2.fq.gz 0.2 > sub20_2.fq.gz
```
3.3 SPAdes 

3.3.1 Assemble the concatenated sequences with SPAdes
```
spades.py -meta --only-assembler -t 100 -m 800 -1 sub20_1.fq.gz -2 sub20_2.fq.gz -o assembly_20_percent
```
3.3.2 Check the assembly accuracy with seqkit. The contigs.fa file can be found within the output folder of the assembly.
```
seqkit stats contigs.fa
```

3.3.3 Check BUSCO completeness. The lineage nematoda_odb10 is downloaded by the program automatically and does not need to be installed locally.
```
busco -i contigs.fa -o busco_EPT_odb10 -l nematoda_odb10 -m genome -c 32
```

3.4 Platanus_allee

3.4.1 Assemble the concatenated sequences with Platanus_allee
```
./platanus_allee assemble -o Platanus_all.out -f sub50_1.fq.gz sub50_2.fq.gz -m 800 -t 64
```
3.4.2 Check the assembly accurary with seqkit similar to step 3.3.2


**4. Blobtools**

4.1 Blast the contig file from the De-novo Genome Assembly against the nucleotide database
```
blastn -task megablast -query /home/jcaroval/02.TrimmedData/assembly_20_percent/contigs.fa -db /import/lragioneri/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -out /home/jcaroval/02.TrimmedData/blastn/megablast.out
```
4.2 Run minimap to receive a coverage file
```
minimap2 -ax sr -t 16 ./assembly_20_percent/contigs.fa sub20_1.fq.gz sub20_2.fq.gz | samtools sort -@16 -O BAM -o contigs.fasta.bam
```
4.3 Run Blobtools
```
blobtools create --fasta ./assembly_20_percent/contigs.fa --hits ./blastn/megablast.out --cov contigs.fasta.bam --taxid 55786 --taxdump ./taxdump blob_EPT
```
4.4 Visualize Blobplot
```
blobtools view --remote --view blob blob_EPT
```
