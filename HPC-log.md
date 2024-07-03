####Content:
• 45 Illumina Short-Read (PE150) Samples (Acrobeles emmatus)
• 1 Reference Genome (Acrobeloides tricornis)
• 1 UCE file (based on 7 different Panagrolaimus sp.)

####Tools used:
• 

####Data Curation and Quality
1. Calculate and compare MD5sums. They should resemble the MD5sums that were provided by Novogene.
```
find . -type f -name '*.fq.gz' | while read -r file; do
echo "current file: $file"
md5sum "$file" >> md5sums.txt
done
```
2. Create MultiQC files to investigate the quality of the raw files.
```
multiqc *.fq.gz
```
3. Trim the sequences. In this case, 1 is always the forward, 2 the reverse file.
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
4. Repeat step 2 (MultiQC) on the trimmed files.

####Reference-based Genome Alignments
1. 

