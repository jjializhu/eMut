######
###   for example (BM89_M5)
###################################  Monopogen  ###########################################################################################
######     1) mutation detection
path="/cluster2/huanglab/jzhu/app/Monopogen" # where Monopogen is downloaded
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

python  ${path}/src/Monopogen.py  preProcess -b /cluster2/huanglab/jzhu/app/Monopogen/test/bam.lst \
 -o /cluster2/huanglab/jzhu/2-AML/Monopogen  -a ${path}/apps -t 8

python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 8   -r  /cluster2/huanglab/jzhu/app/Monopogen/test/region.lst \
    -p  /cluster2/huanglab/jzhu/app/Monopogen/example/ \
    -g  /cluster2/huanglab/jzhu/app/Monopogen/example/chr20_2Mb.hg38.fa   -m 3 -s all  -o /cluster2/huanglab/jzhu/2-AML/Monopogen

python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  /cluster2/huanglab/jzhu/app/Monopogen/test/region.lst  -t 50 \
    -i  /cluster2/huanglab/jzhu/2-AML/Monopogen  -l  /cluster2/huanglab/jzhu/app/Monopogen/example/CB_7K.maester_scRNA.csv   -s featureInfo     \
    -g /cluster2/huanglab/jzhu/app/Monopogen/example/GRCh38.chr20.fa

python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  /cluster2/huanglab/jzhu/app/Monopogen/test/region.lst  -t 22  -w 10MB \
    -i /cluster2/huanglab/jzhu/2-AML/Monopogen  -l /cluster2/huanglab/jzhu/app/Monopogen/example/CB_7K.maester_scRNA.csv -s cellScan     \
    -g /cluster2/huanglab/jzhu/app/Monopogen/example/GRCh38.chr20.fa


######      2) combine SNV for each sample
# SNV_extract.ipynb


######      2) mutations annotation
sort -k1,1 -k2n /cluster2/huanglab/jzhu/2-AML/Monopogen/summary/SNV.vcf > /cluster2/huanglab/jzhu/2-AML/Monopogen/summary/SNV.sorted.vcf
vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol \
    --numbers --domains --gene_phenotype --canonical --protein --biotype --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length \
    --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir /cluster2/huanglab/jzhu/ReferenceGenome/VEP_GRCh38/ \
    --fasta /cluster2/huanglab/jzhu/ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa  \
    --input_file /cluster2/huanglab/jzhu/2-AML/Monopogen/summary/SNV.sorted.vcf \
    --output_file /cluster2/huanglab/jzhu/2-AML/Monopogen/summary/SNV.vep.vcf \
    --polyphen b --af --af_1kg --af_esp --regulatory --fork 10 --force_overwrite --buffer_size  10000

perl /cluster2/huanglab/jzhu/app/vcf2maf-1.6.21/vcf2maf.pl --input-vcf /cluster2/huanglab/jzhu/2-AML/Monopogen/summary/SNV.vep.vcf \
                --output-maf /cluster2/huanglab/jzhu/2-AML/Monopogen/summary/SNV.vep.maf \
                --vep-path /cluster2/huanglab/jzhu/app/miniconda3/envs/VEP/bin/ \
                --vep-data /cluster2/huanglab/jzhu/ReferenceGenome/VEP_GRCh38/ \
                --ref-fasta /cluster2/huanglab/jzhu/ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
                --ncbi-build GRCh38 --vep-overwrite --inhibit-vep 

#####################################  GATK Mutect2    ##############################################################################################
#####      1) mutation detection
python /cluster2/huanglab/jzhu/2-AML/code1/GATK/1.run_GATK.py \
    -b /cluster/huanglab/mzhu/projects/aml/BM89_ATAC/BM89ATAC_dedup.bam \
    -o /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5 \
    -p /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/peaks.bed \
    -g /cluster/huanglab/mzhu/reference/UCSC/chrom_sizes/hg38.chrom.sizes_flt.bed

#####     2) mutation annotation(VEP)
cd /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5/vcf/single_filter_vcf/
for i in `ls *.filter.vcf`; do echo ${i%%.*} >> /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5/cells.txt; done
nohup python /cluster2/huanglab/jzhu/2-AML/code1/GATK/2.mutation_annotation.py \
    -input  /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5/cells.txt \
    -mafDir /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5/maf \
    -vepDir /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5/vep \
    -vcfDir /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/BM89_M5/vcf/single_filter_vcf > /dev/null &

#####     3)combine maf file of single-cell for each sample
sample="BM89_M5"
cd /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/$sample/maf
for i in `ls *.maf`; do sed -n '3,$p' $i >> /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/$sample/result.maf; done
cat /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/header.maf /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/$sample/result.maf > /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/$sample.maf
rm /cluster2/huanglab/jzhu/2-AML/2.GATK/all-scATAC/$sample/result.maf
