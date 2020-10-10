*******
# Probiotics-relieve-human-stress-and-anxiety-project
Probiotic consumption relieved human stress and anxiety symptoms via modulating the gut microbiota and neuroactive potential
******

## Bioinformatics software and key codes for metagenomic sequencing analysis 
### 1. ***For quality control and de_hosting: KneadData and Bowtie2***
```
klab_metaqc qc -s data.list -t human -j 10 -o qc_result_folder_rmhuman -f # klab_metaqc is a simple process combining the above two tools
```

### 2. ***For assembling: Megahit***
```
ls -d *| parallel -j 5 megahit -m 0.5 --min-contig-len 200 -t 10 --out-dir {}_Output_ass --out-prefix {} -1 {}/{}.rmhost.r1.fq.gz -2 {}/{}.rmhost.r2.fq.gz ::: * # parallel code
```

### 3. ***For genome Binning: MaxBin2 and MetaBAT2***
```
ls -d Sample* | parallel -j 3 metawrap binning -a Assembly_result/{}_Output_ass/{}_Output_contigs_min2000.fasta -o {}/{}.bins -t 16 --metabat2 --maxbin2 --concoct -l 2000 {}/{}.rmhost_1.fastq {}/{}.rmhost_2.fastq # parallel code
```

### 4. ***For bin_refinement: MetaWRAP***
```
ls -d *bin | parallel -j 5 metawrap bin_refinement -t 5 -m 200 -c 50 -x 10 -A {}/maxbin2_bins -B {}/metabat2_bins -o {}/metawrap # parallel code
````

### 5. ***For genome de_replicate: dRep***
```
cat bins_80_5_checkm.summary | awk 'BEGIN{OFS=",";}{print $1 ".fa", $2, $3}' > ./bins_80_5_checkm_info
sed 's/bin.fa,/genome,/g' bins_80_5_checkm_info > bins_80_5_checkm_info.csv
dRep dereplicate dereplicate_result -pa 0.95 -sa 0.95 -g ./bins_80_5/* -p 32 --genomeInfo bins_80_5_checkm_info.csv
```

### 6. ***For genome quality evaluate: checkM***
```
checkm lineage_wf -x fa -t 32 --pplacer_threads 32 bins_80_5 bins_80_5_checkm
summarize_checkm.py bins_80_5_checkm/storage/bin_stats_ext.tsv > bins_80_5_checkm.summary
```

### 7. ***For gene predicting and genome protein annotation: Prodigal and diamond(blastp)***
```
ls -d S* | parallel -j 4 bwa index -a bwtsw {}/{}_contigs_min2000.fasta
ls -d S* | parallel -j 4 bwa mem -t 10 {}/{}_Output_contigs_min2000.fasta ../../0.1_clean_data/all_sample/{}/{}.rmhost_1.fastq.gz ../../0.1_clean_data/all_sample/{}/{}.rmhost_2.fastq.gz -o {}/{}.sam 
for i in `ls -d S*`; do samtools view -bS -@ 12 ${i}/${i}.sam | samtools sort - -o ${i}/${i}.sort.bam; done 
diamond blastp --threads 32 --max-target-seqs 10 --db  /nvmessdnode3/opt/database/uniport/uniprot_trembl_sport.dmnd --query all_combine.faa --outfmt 6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out all_combine.dia
```



### 8. ***For map raw reads to the scaffolds and identify of neuroactive compounds: BBMap***
````
diamond blastp --threads 32 --max-target-seqs 10 --db  /nvmessdnode3/opt/database/uniport/uniprot_trembl_sport.dmnd --query all_combine.faa --outfmt 6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out all_combine.dia
####java -jar -Xmx500g /nvmessdnode3/opt/software/omixer/omixer-rpm-1.1/omixer-rpm-1.1.jar -d /nvmessdnode3/opt/software/omixer/omixer-rpm-1.1/new_pipeline_452/new_pathway.f  -i bin_kegg_matrix -t 20 -o omixer_out  -e 2 -c 0.66
```
