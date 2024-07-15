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
md5sum "$file" >> md5sums.txt;
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
qualimap bamqc \
-bam ID_sorted.bam \
-gff busco_nematodes_all.gff
````
2.2 Mapping with nextgenmap

2.2.1 Map the sequences against the reference genome
```
ngm \
-r PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa \
-1 ID_trimmed_1.fq.gz \
-2 ID_trimmed_2.fq.gz \
-o ID_out.sam \
-t 64
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
spades.py \
-meta \
--only-assembler \
-t 100 \
-m 800 \
-1 sub20_1.fq.gz \
-2 sub20_2.fq.gz \
-o assembly_20_percent
```
3.3.2 Check the assembly accuracy with seqkit. The contigs.fa file can be found within the output folder of the assembly.
```
seqkit stats contigs.fa
```

3.3.3 Check BUSCO completeness. The lineage nematoda_odb10 is downloaded by the program automatically and does not need to be installed locally.
```
busco \
-i contigs.fa \
-o busco_EPT_odb10 \
-l nematoda_odb10 \
-m genome \
-c 32
```

3.4 Platanus_allee

3.4.1 Assemble the concatenated sequences with Platanus_allee
```
./platanus_allee assemble \
-o Platanus_all.out \
-f sub50_1.fq.gz sub50_2.fq.gz \
-m 800 \
-t 64
```
3.4.2 Check the assembly accurary with seqkit similar to step 3.3.2


**4. Blobtools**

4.1 Blast the contig file from the De-novo Genome Assembly against the nucleotide database
```
blastn \
-task megablast \
-query /home/jcaroval/02.TrimmedData/assembly_20_percent/contigs.fa \
-db /import/lragioneri/nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-out /home/jcaroval/02.TrimmedData/blastn/megablast.out
```
4.2 Run minimap to receive a coverage file
```
minimap2 -ax sr -t 16 ./assembly_20_percent/contigs.fa sub20_1.fq.gz sub20_2.fq.gz | samtools sort -@16 -O BAM -o contigs.fasta.bam
```
4.3 Run Blobtools
```
blobtools create \
--fasta ./assembly_20_percent/contigs.fa \
--hits ./blastn/megablast.out \
--cov contigs.fasta.bam \
--taxid 55786 \
--taxdump ./taxdump blob_EPT
```
4.4 Visualize Blobplot
```
blobtools view \
--remote \
--view blob blob_EPT
```

**5. Mapping against UCEs**

This pipeline follows this protocol: https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#aligning-uce-loci

5.1 Assemble the data. Here, only a subset of 6 samples was assembled.
```
phyluce_assembly_assemblo_spades \
--conf assembly.conf \
--output spades-assemblies \
--core 64 \
--mem 800
```
5.2 Match contigs to probes
```
phyluce_assembly_match_contigs_to_probes \
--contigs spades-assemblies/contigs \
--probes Panagrolaimus1-v1-master-probe-list-DUPE-SCREENED.fasta \
--output UCEs --min-coverage 30
```
5.3 Get match counts
```
phyluce_assembly_get_match_counts \
--locus-db UCEs/probe.matches.sqlite \
--taxon-list-config taxon-set.conf \
--taxon-group 'all' \
--incomplete-matrix \
--output taxon-sets/all-taxa-incomplete.conf
```
5.4 Get fastas from match counts
```
phyluce_assembly_get_fastas_from_match_counts \
--contigs ../../spades-assemblies/contigs \
--locus-db ../../UCEs/probe.matches.sqlite \
--match-count-output all-taxa-incomplete.conf \
--output all-taxa-incomplete.fasta \
--incomplete-matrix all-taxa-incomplete.incomplete \
--log-path log
```
5.5 Get individual fasta files from all-taxa-incomplete.fasta
```
phyluce_assembly_explode_get_fastas_file \
--input all-taxa-incomplete.fasta \
--output exploded-fastas \
--by-taxon
```
5.6 Get summary stats on the individual fasta files
```
for i in exploded-fastas/*.fasta;
do
  phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```
5.7 Investigate the amount of contigs
```
seqkit stats all-taxa-incomplete.fasta
```

5.8 Remove duplicates from all-taxa-incomplete.fasta
```
seqkit sort –n all-taxa-incomplete-length.fasta > all-taxa-incomplete-sorted.fasta
cut -d ‘_’ -f 1 all-taxa-incomplete-sorted.fasta > all-taxa-incomplete-sorted-shortened.fasta
seqkit rmdup -n < all-taxa-incomplete-sortened-shortened.fasta > all-taxa-incomplete-no-dups.fasta
```
5.9 Check the amount of contigs again
```
seqkit stats all-taxa-incomplete-no-dups.fasta
```
5.10 Map individual sequences against the reference set with bwa-mem2

5.10.1 Create index for reference set
```
bwa-mem2 index all-taxa-incomplete-no-dups.fasta
```
5.10.2 Map against reference set
```
bwa-mem2 mem \
-t 32 \
./10.UCE_index/all-taxa-incomplete-no-dups.fasta \
./02.TrimmedData/ID/ID_trimmed_1.fq.gz \
./02.TrimmedData/ID/ID_trimmed_2.fq.gz \
-o 09.UCE/ID_UCEs_out.sam
```
or for multiple files:
```
for INDIVIDUAL_DIR in ./02.TrimmedData/EPT_*; do
  if [-d "$INDIVIDUAL_DIR" ]; then
      INDIVIDUAL=$(basename "$INDIVIDUAL_DIR")

      READ1="$INDIVIDUAL_DIR/*_1.fq.gz"
      READ2="$INDIVIDUAL_DIR/*_2.fq.gz"

      if ls $READ1 1> /dev/null 2>&1 && ls $READ2 1> /dev/null 2>&1; then
        OUTPUT_FILE="./09.UCE/bwa-mem2-UCE/${INDIVIDUAL}_UCEs_out.sam"

          bwa-mem2 mem -t 32 ./10.UCE_index/all-taxa-incomplete-no-dups.fasta $READ1 $READ2 -o $OUTPUT_FILE
          echo "Finished bwa-mem2 for $INDIVIDUAL."
      else
        echo "Skipped $INDIVIDUAL due to missing read files."
      fi
  fi
done
```

5.10.3 Set a flag-tag for each individual
```
SAM_DIR="/home/jcaroval/09.UCE/bwa-mem2-UCE"
OUTPUT_DIR="/home/jcaroval/09.UCE/bwa-mem2-UCE"
for input_sam in ${SAM_DIR}/*.sam; do
  filename=$(basename -- "$input_sam")
  sample=$(echo $filename | cut -d'_' -f1,2)
  output_sam="${OUTPUT_DIR}/${sample}_RG.sam"

  picard AddOrReplaceReadGroups \
       I=$input_sam \
       O=$output_sam \
       RGID=$sample \
       RGLB=lib${sample} \
       RGPL=illumina \
       RGPU=unit${sample} \
       RGSM=$sample \
       VALIDATION_STRINGENCY=LENIENT

  echo "Read Group for $sample added successfully."
```


5.10.4 Convert SAM to BAM
```
ls *_RG.sam | sed 's/_RG.sam//g' > list-XX 
while read f; do 
samtools view -b $f"_RG.sam" > $f"_RG.bam" ;
done < list-XX
```

5.10.5 Quality check of the mappings using samtools flagstat
```
for bam_file in *.bam
do
    sample_id=$(basename "$bam_file" .bam)
    echo "Sample ID: $sample_id" >> flagstat.out
    samtools flagstat "$bam_file" >> flagstat.out
    echo "" >> flagstat.out
done
```
or using samtools stats
```
for bam_file in *.bam
do
    sample_id=$(basename "$bam_file" .bam)
    echo "Sample ID: $sample_id" >> samtools.out
    samtools stats "$bam_file" | grep ^SN | cut -f 2- >> samtools.out
    echo "" >> samtools.out
done
```

5.11 Sort the BAM file
```
for file in *_RG.bam; do
    sample_id=$(basename "$file" .bam)
    sorted_file="${file}_sorted.bam"
    samtools sort "$file" -o "$sorted_file"
```

5.12 Extract only the mapped reads from the BAM files
```
for sorted_file in *_sorted.bam; do
    sample_id=$(basename "$sorted_file" .bam)
    output="${sorted_file}_mapped.bam"
    samtools view -q 30 -F 4 -b "$sorted_file" > "$output"
done
```

5.13 Repeat QC with flagstat to make sure that only mapped reads are included (100 % mapped reads)
```
for bam_file in *_RG_sorted_mapped.bam
do
    sample_id=$(basename "$bam_file" _RG_sorted_mapped.bam)
    echo "Sample ID: $sample_id" >> flagstat_mapped.out
    samtools flagstat "$bam_file" >> flagstat_mapped.out
    echo "" >> flagstat_mapped.out
done
```

5.14 Merge all BAM files into one
```
samtools merge -@ 64 -r merged.bam *_sorted_mapped.bam
```

5.15 Index the merged BAM file
```
samtools index merged.bam
```

5.14 Visualize the mapping in Tablet
```
./tablet /home/jcaroval/09.UCE/bwa-mem2-UCE/merged.bam /home/jcaroval/10.UCE_index/all-taxa-incomplete-no-dups.fasta
```
or on IVG.

5.15 Calculate the coverage
```
samtools depth -a merged.bam > coverage.txt
```

5.16 Split the coverage.txt file for each UCE and plot the coverage
```
coverage_file="coverage.txt"
output_dir="./coverage_no_slidingwindow"
awk '{ print > "'${output_dir}'/"$1"_coverage.txt" }' "${coverage_file}"
for uce_coverage_file in "${output_dir}"/*_coverage.txt; do
    uce_name=$(basename -- "${uce_coverage_file}" "_coverage.txt")
    gnuplot <<EOF
    set terminal png size 800,600
    set output '${output_dir}/${uce_name}_coverage_plot.png'
    set xlabel 'Position'
    set ylabel 'Coverage'
    plot '${uce_coverage_file}' using 2:3 with lines title '${uce_name}'
EOF
done
```
then, investigate the plots

5.17 (optional) Apply a sliding window approach to smoothen the data
```
coverage_file="coverage.txt"

output_dir="./output"
smoothed_output_dir="./smoothed_output"

awk '{ print > "'${output_dir}'/"$1"_coverage.txt" }' "${coverage_file}"

window_size=50
step_size=5

for uce_coverage_file in "${output_dir}"/*_coverage.txt; do
    uce_name=$(basename -- "${uce_coverage_file}" "_coverage.txt")
    
    python - <<EOF
import pandas as pd
import os

input_file = "${uce_coverage_file}"
output_file = "${smoothed_output_dir}/${uce_name}_smoothed_coverage.txt"

df = pd.read_csv(input_file, sep="\t", header=None, names=["UCE", "Position", "Coverage"])

window_size = ${window_size}
step_size = ${step_size}
smoothed_data = df['Coverage'].rolling(window=window_size, min_periods=1).mean().iloc[::step_size]
smoothed_positions = df['Position'].iloc[::step_size]

smoothed_df = pd.DataFrame({'Position': smoothed_positions, 'Smoothed_Coverage': smoothed_data})
smoothed_df.to_csv(output_file, sep="\t", index=False, header=["Position", "Smoothed_Coverage"])
EOF

    gnuplot <<EOF
    set terminal png size 800,600
    set output '${smoothed_output_dir}/${uce_name}_coverage_plot.png'
    set xlabel 'Position'
    set ylabel 'Coverage'
    plot '${uce_coverage_file}' using 2:3 with lines title 'Original Coverage', \
         '${smoothed_output_dir}/${uce_name}_smoothed_coverage.txt' using 1:2 with lines title 'Smoothed Coverage'
EOF
```

5.18 Retrieve coverage and abundance
```
ls *merged.bam | parallel -j 5 'samtools depth {} > {}.depth'
```

5.19 Sort the UCEs by their average coverage
```
input_file="merged.bam.depth"
output_file="merged_average_coverages.txt"
awk '{ sum[$1] += $3; count[$1]++ } END { for (uce in sum) print uce, sum[uce] / count[uce] }' "${input_file}" > "${output_file}"
```
and sort them by value
```
sort -n -k2 "merged_average_coverages.txt" > "merged_sorted_average_coverages.txt"
```
Here, a threshold was set to an average coverage of 360. This threshold was chosen because a minimum threshold of 8X per individual was desired. With 45 individuals, it can be assumed that UCEs with a coverage below 360 do not lead to sufficient results.

5.20 Remove the low-coverage UCEs from the files.

```
awk '$2 < 360 {print $1}' "merged_average_coverages.txt" > "uce_list_below_360.txt"
```
… for the merged file
```
samtools view -H merged.bam > header.sam
samtools view -@ 64 merged.bam | grep -v -F -f uce_list_below_360.txt | cat header.sam - | samtools view -b -@ 64 > merged_filtered.bam

rm header.sam
```
… and for all individual files
```
uce_list="uce_list_below_360.txt"
bam_dir="."
num_threads=64
for bam_file in "$bam_dir"/*_RG_sorted_mapped.bam; do
    sample_id=$(basename "$bam_file" .bam)
    output="${bam_dir}/${sample_id}_filtered.bam"

    samtools view -H "$bam_file" > header.sam
    samtools view -@ "$num_threads" "$bam_file" | grep -v -F -f "$uce_list" | cat header.sam - | samtools view -b -@ "$num_threads" > "$output"

    rm header.sam
done
```

5.21 Count the UCEs to verify that the correct number of UCEs had been removed

… for the merged file
```
samtools view -@ 64 "merged_filtered.bam" | awk '{print $3}' | sort | uniq | wc -l
```
… and for all individual files
```
uce_list="uce_list_below_360.txt"
bam_dir="."
num_threads=64
result_file="uce_counts.txt"

echo -e "Sample_ID\tBefore_Filtering\tAfter_Filtering" > "$result_file"

for bam_file in "$bam_dir"/*_RG_sorted_mapped.bam; do
    sample_id=$(basename "$bam_file" .bam)
    output="${bam_dir}/${sample_id}_filtered.bam"
    uce_count_before=$(samtools view -@ "$num_threads" "$bam_file" | awk '{print $3}' | sort | uniq | wc -l)
    uce_count_after=$(samtools view -@ "$num_threads" "$output" | awk '{print $3}' | sort | uniq | wc -l)
    echo -e "${sample_id}\t${uce_count_before}\t${uce_count_after}" >> "$result_file"
done
```

5.22 Extract flanking regions (defined by coverage 8X ≤ N ≤ 50X)
