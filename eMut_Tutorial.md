
# Tutorial for eMut

## Introduction
The eMut, an integrated pipeline for detecting, imputing, and characterizing non-coding mutations in CREs with functional consequences at the single-cell level.

![image](https://github.com/jjializhu/eMut/blob/main/Figures/eMut_workflow.png)


## workflow
Briefly, eMut consists of two main modules: mutation detection and functional interpretation. <br />
Step1. Mutation detection: eMut detects mutations in each cell by implementing methods such as Monopogen or GATK by single-cell chromatin accessibility data. <br />
Step2. Mutation imputation (optional): Given the sparse of scATAC-seq data, we further imputed candidate mutated cells by network propagation using mutated cells (seed cells) in cell-cell similarity graph. <br />
Step3. Functional interpretation: <br />
1) recognize cell type-specific or lineage-specific mutations; 
2) identify hypermutated CREs with significant excess of mutations to characterize potentially important enhancers; 
3) predict the effects of mutations on transcription factor motifs (loss or gain);
4) compare target gene expression changes between mutated cells (or samples) and wild-type. 

### Step1. Mutation detection
Input: Bam-format files of scATAC-seq data<br>
As an example, the Monopogen implementation of somatic mutation prediction for a single sample (without matched normal samples), as well as the simultaneous detection of somatic and germline mutations based on GATK Mutect2, are shown here.
```
/eMut/1.run_GATK.py
/eMut/2.mutation_annotation.py
/eMut/Mutation_detection.sh
```

### Step2. Mutation imputation (optional)
Input: 
(1) mutation profile: mutation-by-cell matrix; <br>
(2) scATAC-seq data: peak-by-cell matrix or knnGraph; <br>

```r
library(graphics)
library(ggforce)
library(scales)
library(ggpubr)
library(ggplot2)
library(pbapply)
source("/eMut/R/function/functions.R")

load("./TileMatrix.Rdata")
load("./SNVMat.Rdata")

TRS.list<-SNVImputation(countMatrix=NULL,
                        knnGraph=mutualknn30,
                        SNVMatrix=mat,
                        mutations=row.names(mat),
                        numk=30,
                        queryCell_cutoff=5,
                        ncors=5)
```

### Step3. functional interpertaion
#### (1) cell type enrichment
Input: 
(1) mutation profile: mutation-by-cell matrix (raw or imputed); <br>
(2) scATAC-seq data: ArchR or signac object; <br>
Here is an example demonstration with an ArchR object .
```r
library(ArchR)
source("/eMut/R/function/functions.R")

proj<- loadArchRProject(path = "./ArchR", force = FALSE, showLogo = FALSE)
cells<-row.names(getCellColData(proj))
cellTypes<-proj$NamedClust
```
#####  cell type enrichment for mutation profile (raw)
```r
mut.type<-pblapply(row.names(mat),function(x){
    mutCells<-colnames(mat)[which(mat[x,]=="0/1" | mat[x,]=="1/1")]
    result<-cellTypeEnrich(cells,
                         cellTypes,
                         mutCells,
                         mutCells_cutoff=50,
                         cellType_cutoff=50)
    if(!is.null(result)){
      result$mut<-rep(x,nrow(result))
      result  
    }
}) %>% rbindlist()
mut.type$p<-p.adjust(mut.type$p)
```

##### cell type enrichment for mutation profile (imputed)
```r
mut.type<-pblapply(names(TRS.list),function(x){
    mutCells<-row.names(TRS.list[[x]])[TRS.list[[x]]$true_cell_top_idx==TRUE]
    result<-cellTypeEnrich(cells,
                         cellTypes,
                         mutCells,
                         mutCells_cutoff=50,
                         cellType_cutoff=50)
    if(!is.null(result)){
      result$mut<-rep(x,nrow(result))
      result  
    }
}) %>% rbindlist()
mut.type$p<-p.adjust(mut.type$p)
```

####  (2) hyperMutated CREs
Modified activeWGS XXXXXXXXXXXXXXX(method)
Input: 
(1) mutation profile: annotated mutation file(VCf/maf format); <br>
(2) peak file: The genomic location of chromatin accessible region; <br>

```r
library(ActiveDriverWGS)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
source("/eMut/R/function/ADWGS_test.r")
source("/eMut/R/function/activeDriverWGS.r")
source("/eMut/R/function/fix_all_results.r")

##  load  mutation profile
mut.df<-data.table::fread(
    file = "./AML10.maf",
    sep = "\t",stringsAsFactors = FALSE,verbose = FALSE,data.table = TRUE,
    showProgress = TRUE,header = TRUE,fill = TRUE,
    skip =1,quote = "")
mut.df<-mut.df[mut.df$FILTER=="PASS",]
mut.df<-mut.df[,c("Chromosome","Start_Position","End_Position",
                             "Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
colnames(mut.df)<-c("chr","pos1","pos2","ref","alt","patient")
mut.df<-mut.[df$chr!="chrM",]

##  load peak
peaks.df<-read.table("./peaks.bed",header=F,quote="",sep="\t")
colnames(peaks.df)<-c("chr","start","end","id")
openRegions<-GenomicRanges::GRanges(seqnames=peaks.df$chr, 
                                    IRanges::IRanges(peaks.df$start,peaks.df$end))
##  identify hypermutated CREs
hyperMut<-ActiveDriverWGS(mutations = mut.df,
                         elements = peaks.df,
                         ref_genome = "hg38",
                         mc.cores=4,
                         window_size=1000000,  ### window
                         detect_depleted_mutations=FALSE,
                         openRegions = openRegions,
                         recovery.dir=paste0("./tmp/",x))

```
####  (3) prediction of TF binding motif change
Input: 
(1) mutation profile: mutation file(VCf/maf format); <br>

```r
library(motifbreakR)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

##   change mutation format
mut<-data.frame(Chromosome=SNV.df$CHROM,Start_Position=SNV.df$POS-1,End_Position=SNV.df$POS,
                names=gsub(";",":",SNV.df$ID),score=rep(0,nrow(SNV.df)),strand=rep("+",nrow(SNV.df)))
mut<-mut %>% distinct(names,.keep_all = TRUE)
write.table(mut,file="./Monopogen/summary/SNVsForMotifBreakR.bed",sep="\t",
            col.names = F,row.names = F,quote=F)

##   motif change prediction
data(motifbreakR_motif)
ENCODE<-subset (motifbreakR_motif, dataSource=="ENCODE-motif" & organism=="Hsapiens")
snps.mb.frombed <- snps.from.file(file = "./Monopogen/summary/SNVsForMotifBreakR.bed",
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")

SNV.motifs <- motifbreakR(snpList = snps.mb.frombed, filterp = TRUE,
                       pwmList = ENCODE,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam(5))
```

####  (4) Comparsion of target gene expression between mutated samples(cells) and wild-type
Input: 
(1) mutation profile: Combined mutations for all samples (VCf/maf format, ); <br>
(2) 

```r


```




