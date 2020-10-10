*******
# Probiotics-relieve-human-stress-and-anxiety-project
Probiotic consumption relieved human stress and anxiety symptoms via modulating the gut microbiota and neuroactive potential
******

## Bioinformatics software and key codes for metagenomic sequencing analysis 
### ***For quality control and de_hosting: KneadData and Bowtie2***
```
klab_metaqc qc -s data.list -t human -j 10 -o qc_result_folder_rmhuman -f # klab_metaqc is a simple process combining the above two tools
```

### ***For assembling: Megahit***
```
ls -d *| parallel -j 5 megahit -m 0.5 --min-contig-len 200 -t 10 --out-dir {}_Output_ass --out-prefix {} -1 {}/{}.rmhost.r1.fq.gz -2 {}/{}.rmhost.r2.fq.gz ::: * # parallel code
```

