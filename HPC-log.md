<h1> Aim </h1>  
<h1>Data</h1> 
<ul> • 45 Illumina short-read (PE150) samples (<i>Acrobeles emmatus</i>) </ul> 
<ul> • 1 Reference genome (<i>Acrobeloides tricornis</i>) </ul> 
<ul> • 1 UCE probe list (containing probes for 1606 different uces) </ul> 
Note: The probes were designed by Laura Villegas based on 7 different <i>Panagrolaimus</i> sp.

<h1>Tools</h1> 
<ul> • <tt>bcftools (V1.18)</tt> | Danecek et al., 2021</ul>
<ul> • <tt>BLAST+ (V2.14.1)</tt> | Camacho et al., 2009</ul>     
<ul> • <tt>blobtoolkit (V4.3.5)</tt> | Challis et al, 2019</ul>     
<ul> • <tt>BUSCO (V5.7.1)</tt> | Manni et al., 2021</ul> 
<ul> • <tt>bwa-mem2 (V2.2.1)</tt> | Vasimuddin et al., 2019</ul> 
<ul> • <tt>Conda (V23.3.1/V4.9.2)</tt></ul> 
<ul> • <tt>fastp (V0.23.2)</tt> | Chen et al., 2018</ul> 
<ul> • <tt>IGV (V2.16.2)</tt> | Robinson et al., 2011</ul> 
<ul> • <tt>minimap2 (V2.26-r1175)</tt> | Li, 2018</ul> 
<ul> • <tt>MultiQC (V1.0.dev0)</tt> | Ewels et al., 2016</ul> 
<ul> • <tt>nextgenmap (V0.5.5)</tt> | Sedlazeck et al., 2013</ul> 
<ul> • <tt>Phyluce (V1.7.3)</tt> | Faircloth, 2016</ul> 
<ul> • <tt>Picard (V2.27.5)</tt> | Broad Institute, 2019</ul> 
<ul> • <tt>Pixy (V1.2.10.beta2)</tt> | Korunes and Samuk, 2021</ul> 
<ul> • <tt>Platanus-allee (V2.2.2)</tt> | Kajitani et al., 2019</ul> 
<ul> • <tt>QualiMap (V2.2.2-dev)</tt> | García-Alcalde, 2012</ul> 
<ul> • <tt>samtools (V1.18)</tt> | Danecek et al., 2021</ul> 
<ul> • <tt>seqkit (V2.5.1)</tt> | Shen et al., 2016</ul> 
<ul> • <tt>seqtk (V1.2-r94)</tt> | Li, 2012</ul> 
<ul> • <tt>SPAdes (V3.13.0)</tt> | Prjibelski et al., 2020</ul> 
<ul> • <tt>Tablet (V1.21.02.08)</tt> | Milne et al., 2013</ul> 


<h1>1. Data Curation and Quality</h1> 

<h2>1.1 Download and unzip the provided files</h2> 

```
wget URL
tar -xvf filename.tar
```

<h2>1.2 Calculate and compare MD5sums.</h2> 

They MD5sums resemble those that were provided by Novogene.

```
find . -type f -name '*.fq.gz' | while read -r file; do
echo "current file: $file"
md5sum "$file" >> md5sums.txt;
done
```
<h2> 1.3 Create MultiQC files</h2> 

This is done to investigate the quality of the raw files before trimming.
```
multiqc *.fq.gz
```
<h2> 1.4 Trim the sequences. </h2> 

In this case, 1 is always the forward, 2 the reverse file.
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
<h2> 1.5 Repeat step 1.3 (MultiQC) on the trimmed files.</h2> 

```
multiqc *_trimmed_*.fq.gz
```

<h1> 2. Reference-based Genome Alignments </h1> 

Since there is no reference genome of <i>Acrobeles emmatus</i> yet, a cross-species mapping against of <i>Acrobeles emmatus</i> against <i>Acrobeloides tricornis</i> was proceeded.

<h2> 2.1 Mapping with bwa-mem2</h2> 

<h3> 2.1.1 Create an index folder of the reference genome</h3> 

```
bwa-mem2 index PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa
```
<h3>2.1.2 Map the sequences against the reference genome</h3>

```
bwa-mem2 mem -t 32 PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa ID_trimmed_1.fq.gz ID_trimmed_2.fq.gz -o ID_out.sam
```
<h3>2.1.3 Convert sam files to sorted bam files</h3>

```
samtools view -S -b ID_out.sam > ID_out.bam
samtools sort ID_out.bam > ID_sorted.bam
```
<h3>2.1.4 Check the mapping accuracy with samtools</h3>

```
samtools flagstat ID_sorted.bam
```
<h3>2.1.5 Check the mapping accuracy with Qualimap</h3>

```
qualimap bamqc \
-bam ID_sorted.bam \
-gff busco_nematodes_all.gff
````
<h2>2.2 Mapping with nextgenmap</h2>

<h3>2.2.1 Map the sequences against the reference genome</h3>

```
ngm \
-r PAP2217_hifi_ont_default.asm.bp.p_ctg.decont.fa \
-1 ID_trimmed_1.fq.gz \
-2 ID_trimmed_2.fq.gz \
-o ID_out.sam \
-t 64
```
<h3>2.2.2 Check the mapping accuracy with samtools and Qualimap similar to steps 2.1.4 and 2.1.5</h3>

<h1>3. De-novo Genome Assembly</h1>

<h2>3.1 Concatenate the sequences of all 45 samples.</h2>

```
cat *_trimmed_1.fq.gz >> EPT_trimmed_1.fq.gz
cat *_trimmed_2.fq.gz >> EPT_trimmed_2.fq.gz
```
<h2>3.2 Downsample the sequences to a desired fraction. This is done to reduce the runtime of used tools.</h2>

```
seqtk EPT_trimmed_1.fq.gz 0.2 > sub20_1.fq.gz
seqtk EPT_trimmed_2.fq.gz 0.2 > sub20_2.fq.gz
```
</h2>3.3 SPAdes</h2> 

<h3>3.3.1 Assemble the concatenated sequences with SPAdes</h3>

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
<h3>3.3.2 Check the assembly accuracy with seqkit. The contigs.fa file can be found within the output folder of the assembly.</h3>

```
seqkit stats contigs.fa
```

<h3>3.3.3 Check BUSCO completeness. The lineage nematoda_odb10 is downloaded by the program automatically and does not need to be installed locally.</h3>

```
busco \
-i contigs.fa \
-o busco_EPT_odb10 \
-l nematoda_odb10 \
-m genome \
-c 32
```

<h2>3.4 Platanus_allee</h2>

<h3>3.4.1 Assemble the concatenated sequences with Platanus_allee</h3>

```
./platanus_allee assemble \
-o Platanus_all.out \
-f sub50_1.fq.gz sub50_2.fq.gz \
-m 800 \
-t 64
```
<h3>3.4.2 Check the assembly accurary with seqkit similar to step 3.3.2</h3>


<h1>4. Blobtools</h1>

<h2>4.1 Blast the contig file from the De-novo Genome Assembly against the nucleotide database</h2>

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
<h2>4.2 Run minimap to receive a coverage file</h2>

```
minimap2 -ax sr -t 16 ./assembly_20_percent/contigs.fa sub20_1.fq.gz sub20_2.fq.gz | samtools sort -@16 -O BAM -o contigs.fasta.bam
```
<h2>4.3 Run Blobtools</h2>

```
blobtools create \
--fasta ./assembly_20_percent/contigs.fa \
--hits ./blastn/megablast.out \
--cov contigs.fasta.bam \
--taxid 55786 \
--taxdump ./taxdump blob_EPT
```
<h2>4.4 Visualize Blobplot</h2>

```
blobtools view \
--remote \
--view blob blob_EPT
```

<h1>5. Mapping against UCEs</h1>

This pipeline follows this protocol: https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#aligning-uce-loci

<h2>5.1 Assemble the data. Here, only a subset of 6 samples was assembled.</h2>

```
phyluce_assembly_assemblo_spades \
--conf assembly.conf \
--output spades-assemblies \
--core 64 \
--mem 800
```
<h2>5.2 Match contigs to probes</h2>

```
phyluce_assembly_match_contigs_to_probes \
--contigs spades-assemblies/contigs \
--probes Panagrolaimus1-v1-master-probe-list-DUPE-SCREENED.fasta \
--output UCEs --min-coverage 30
```
<h2>5.3 Get match counts</h2>

```
phyluce_assembly_get_match_counts \
--locus-db UCEs/probe.matches.sqlite \
--taxon-list-config taxon-set.conf \
--taxon-group 'all' \
--incomplete-matrix \
--output taxon-sets/all-taxa-incomplete.conf
```
<h2>5.4 Get fastas from match counts</h2>

```
phyluce_assembly_get_fastas_from_match_counts \
--contigs ../../spades-assemblies/contigs \
--locus-db ../../UCEs/probe.matches.sqlite \
--match-count-output all-taxa-incomplete.conf \
--output all-taxa-incomplete.fasta \
--incomplete-matrix all-taxa-incomplete.incomplete \
--log-path log
```
<h2>5.5 Get individual fasta files from all-taxa-incomplete.fasta</h2>

```
phyluce_assembly_explode_get_fastas_file \
--input all-taxa-incomplete.fasta \
--output exploded-fastas \
--by-taxon
```
<h2>5.6 Get summary stats on the individual fasta files</h2>

```
for i in exploded-fastas/*.fasta;
do
  phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```
<h2>5.7 Investigate the amount of contigs</h2>

```
seqkit stats all-taxa-incomplete.fasta
```

<h2>5.8 Remove duplicates from all-taxa-incomplete.fasta</h2>

```
seqkit sort –n all-taxa-incomplete-length.fasta > all-taxa-incomplete-sorted.fasta
cut -d ‘_’ -f 1 all-taxa-incomplete-sorted.fasta > all-taxa-incomplete-sorted-shortened.fasta
seqkit rmdup -n < all-taxa-incomplete-sortened-shortened.fasta > all-taxa-incomplete-no-dups.fasta
```
<h2>5.9 Check the amount of contigs again</h2>

```
seqkit stats all-taxa-incomplete-no-dups.fasta
```
<h2>5.10 Map individual sequences against the reference set with bwa-mem2</h2>

<h3>5.10.1 Create index for reference set</h3>

```
bwa-mem2 index all-taxa-incomplete-no-dups.fasta
```
<h3>5.10.2 Map against reference set</h3>

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

<h3>5.10.3 Set a flag-tag for each individual</h3>

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


<h3>5.10.4 Convert SAM to BAM</h3>

```
ls *_RG.sam | sed 's/_RG.sam//g' > list-XX 
while read f; do 
samtools view -b $f"_RG.sam" > $f"_RG.bam" ;
done < list-XX
```

<h3>5.10.5 Quality check of the mappings using samtools flagstat</h3>

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

<h2>5.11 Sort the BAM file</h2>

```
for file in *_RG.bam; do
    sample_id=$(basename "$file" .bam)
    sorted_file="${file}_sorted.bam"
    samtools sort "$file" -o "$sorted_file"
```

<h2>5.12 Extract only the mapped reads from the BAM files</h2>

```
for sorted_file in *_sorted.bam; do
    sample_id=$(basename "$sorted_file" .bam)
    output="${sorted_file}_mapped.bam"
    samtools view -q 30 -F 4 -b "$sorted_file" > "$output"
done
```

<h2>5.13 Repeat QC with flagstat to make sure that only mapped reads are included (100 % mapped reads)</h2>

```
for bam_file in *_RG_sorted_mapped.bam
do
    sample_id=$(basename "$bam_file" _RG_sorted_mapped.bam)
    echo "Sample ID: $sample_id" >> flagstat_mapped.out
    samtools flagstat "$bam_file" >> flagstat_mapped.out
    echo "" >> flagstat_mapped.out
done
```

<h2>5.14 Merge all BAM files into one</h2>

```
samtools merge -@ 64 -r merged.bam *_sorted_mapped.bam
```

<h2>5.15 Index the merged BAM file</h2>

```
samtools index merged.bam
```

<h2>5.14 Visualize the mapping in Tablet</h2>

```
./tablet /home/jcaroval/09.UCE/bwa-mem2-UCE/merged.bam /home/jcaroval/10.UCE_index/all-taxa-incomplete-no-dups.fasta
```
or on IVG.

<h2>5.15 Calculate the coverage</h2>

```
samtools depth -a merged.bam > coverage.txt
```

<h2>5.16 Split the coverage.txt file for each UCE and plot the coverage</h2>

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

<h2>5.17 (optional) Apply a sliding window approach to smoothen the data</h2>

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

<h2>5.18 Retrieve coverage and abundance</h2>

```
ls *merged.bam | parallel -j 5 'samtools depth {} > {}.depth'
```

<h2>5.19 Sort the UCEs by their average coverage</h2>

```
input_file="merged.bam.depth"
output_file="merged_average_coverages.txt"
awk '{ sum[$1] += $3; count[$1]++ } END { for (uce in sum) print uce, sum[uce] / count[uce] }' "${input_file}" > "${output_file}"
```
and sort them by value
```
sort -n -k2 "merged_average_coverages.txt" > "merged_sorted_average_coverages.txt"
```

<h2>5.20 Remove the low-coverage UCEs from the files.</h2>

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

<h2>5.21 Count the UCEs to verify that the correct number of UCEs had been removed</h2>

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

<h1>6. pixy</h1>

This pipeline continues from step 5.13 with the sorted and mapped BAM files.

<h2>6.1 Create bam-list for the EPT</h2>

```
ls *_RG_sorted_mapped.bam > EPT_all_bam_list.txt
```

<h2>6.2 Create pileup files per sampling spot</h2>

```
bcftools mpileup -f ../../../10.UCE_index/all-taxa-incomplete-no-dups.fasta -b EPT_all_bam_list.txt --threads 64 | bcftools call -m -Oz -f GQ -o EPT_all_mpileup_bcftools.vcf --threads 64
```

<h2>6.3 Filter coverage range</h2>

```
bcftools filter -i 'INFO/DP>=8 && INFO/DP<=50' -Oz -o EPT_all_filtered.vcf EPT_all_mpileup_bcftools.vcf --threads 64
```

<h2>6.4 Relocate the DP content</h2>

This had to be done as the default position of DP within the dataset could not be found by pixy.

```
input_vcf="EPT_all_filtered.vcf.gz"
output_vcf="EPT_all_filtered_modified.vcf.gz"
# Create a temporary file
temp_file=$(mktemp)
# Process the VCF file
while read -r line; do
    if [[ $line == \#\#* ]]; then
        # Meta-information lines
        echo "$line" >> "$temp_file"
    elif [[ $line == \#CHROM* ]]; then
        # Header line
        IFS=$'\t' read -r -a headers <<< "$line"
        echo "$line" >> "$temp_file"
    else
        # Data lines
        IFS=$'\t' read -r -a columns <<< "$line"
        info_field=${columns[7]}
        format_field=${columns[8]}
        sample_fields=("${columns[@]:9}")
        # Extract DP from INFO field
        dp_value="."
        new_info_field=""
        IFS=';' read -r -a info_items <<< "$info_field"
        for item in "${info_items[@]}"; do
            if [[ $item == DP=* ]]; then
                dp_value=${item#DP=}
            else
                if [ -n "$new_info_field" ]; then
                    new_info_field+=";"
                fi
                new_info_field+="$item"
            fi
        done
        # Modify FORMAT and sample fields
        if [[ $format_field != DP ]]; then
            format_field+=":DP"
        fi
        for i in "${!sample_fields[@]}"; do
            sample_fields[$i]+=":$dp_value"
        done
        # Write the modified line
        columns[7]=$new_info_field
        columns[8]=$format_field
        columns=("${columns[@]:0:9}" "${sample_fields[@]}")
        echo -e "${columns[*]}" | tr ' ' '\t' >> "$temp_file"
    fi
done < <(zcat "$input_vcf")
# Move the temporary file to the output file
bgzip -c "$temp_file" > "$output_vcf"
rm "$temp_file"
echo "Modified VCF file saved to $output_vcf"
```


<h2>6.5 Zip and index the vcf</h2>

```
bgzip EPT_all_filtered_modified.vcf #optional, as this was already done in the relocation-script 
```

… and run tabix
```
tabix EPT_all_filtered_modified.vcf.gz
```

<h2>6.6 Create populations file</h2>

This file needs to contain all sample IDs in column 1 and their regarding sampling spots in column 2

<h2>6.7 Calculate within population nucleotide diversity (pi).</h2>

The window_size is set to 10000 so that it always calculates for the whole UCE.

```
pixy --stats pi --vcf EPT_all_filtered_modified.vcf.gz --populations ../populations_file.txt --window_size 10000 --n_cores 64 --output_prefix windowsize10000
```

<h2>6.8 Calculate between population nucleotide divergence (dxy)</h2>

```
pixy --stats dxy --vcf EPT_all_filtered_modified.vcf.gz --populations ../populations_file.txt --window_size 10000 --n_cores 64 --output_prefix windowsize10000
```
