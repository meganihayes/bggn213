Class14
================
Megan Hayes
February 27, 2019

``` r
source("http://bioconductor.org/biocLite.R")
```

    ## Bioconductor version 3.4 (BiocInstaller 1.24.0), ?biocLite for help

    ## A new version of Bioconductor is available after installing the most
    ##   recent version of R; see http://bioconductor.org/install

``` r
biocLite()
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.3 (2017-03-06).

    ## Old packages: 'ade4', 'backports', 'boot', 'checkmate', 'class',
    ##   'cluster', 'codetools', 'colorspace', 'data.table', 'evaluate',
    ##   'foreign', 'Hmisc', 'lattice', 'MASS', 'Matrix', 'mgcv',
    ##   'RcppArmadillo', 'RCurl', 'rpart', 'scales', 'stringi', 'stringr',
    ##   'survival', 'tibble', 'utf8', 'xfun', 'XML'

``` r
biocLite("DESeq2")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.3 (2017-03-06).

    ## Installing package(s) 'DESeq2'

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/ww/dtnfg13n74n9q1bnfgwys6540000gn/T//Rtmpp1WNZE/downloaded_packages

    ## Old packages: 'ade4', 'backports', 'boot', 'checkmate', 'class',
    ##   'cluster', 'codetools', 'colorspace', 'data.table', 'evaluate',
    ##   'foreign', 'Hmisc', 'lattice', 'MASS', 'Matrix', 'mgcv',
    ##   'RcppArmadillo', 'RCurl', 'rpart', 'scales', 'stringi', 'stringr',
    ##   'survival', 'tibble', 'utf8', 'xfun', 'XML'

``` r
biocLite("DESeq2")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.3 (2017-03-06).

    ## Installing package(s) 'DESeq2'

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/ww/dtnfg13n74n9q1bnfgwys6540000gn/T//Rtmpp1WNZE/downloaded_packages

    ## Old packages: 'ade4', 'backports', 'boot', 'checkmate', 'class',
    ##   'cluster', 'codetools', 'colorspace', 'data.table', 'evaluate',
    ##   'foreign', 'Hmisc', 'lattice', 'MASS', 'Matrix', 'mgcv',
    ##   'RcppArmadillo', 'RCurl', 'rpart', 'scales', 'stringi', 'stringr',
    ##   'survival', 'tibble', 'utf8', 'xfun', 'XML'

``` r
counts <- read.csv("data/airway_scaledcounts.csv", stringsAsFactors = FALSE)
metadata <-  read.csv("data/airway_metadata.csv", stringsAsFactors = FALSE)
```

``` r
nrow(counts)
```

    ## [1] 38694

``` r
all(metadata$id == colnames(counts)[-1])
```

    ## [1] TRUE

``` r
# find the column names (ids) of the 'control' exps
control.ids <- metadata$id[metadata[,"dex"] =="control"]
head(counts[,control.ids])
```

    ##   SRR1039508 SRR1039512 SRR1039516 SRR1039520
    ## 1        723        904       1170        806
    ## 2          0          0          0          0
    ## 3        467        616        582        417
    ## 4        347        364        318        330
    ## 5         96         73        118        102
    ## 6          0          1          2          0

``` r
control <- counts[,control.ids]
rownames(control) <- counts$ensgene
head(control)
```

    ##                 SRR1039508 SRR1039512 SRR1039516 SRR1039520
    ## ENSG00000000003        723        904       1170        806
    ## ENSG00000000005          0          0          0          0
    ## ENSG00000000419        467        616        582        417
    ## ENSG00000000457        347        364        318        330
    ## ENSG00000000460         96         73        118        102
    ## ENSG00000000938          0          1          2          0

``` r
# find the column names (ids) of the 'treated' exps
treated.ids <- metadata$id[metadata[,"dex"] == "treated"]
treated <- counts[,treated.ids]
rownames(treated) <- counts$ensgene
head(treated)
```

    ##                 SRR1039509 SRR1039513 SRR1039517 SRR1039521
    ## ENSG00000000003        486        445       1097        604
    ## ENSG00000000005          0          0          0          0
    ## ENSG00000000419        523        371        781        509
    ## ENSG00000000457        258        237        447        324
    ## ENSG00000000460         81         66         94         74
    ## ENSG00000000938          0          0          0          0

``` r
## control means
control.mean <- rowSums(control)/ncol(control)
head(control.mean)
```

    ## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 
    ##          900.75            0.00          520.50          339.75 
    ## ENSG00000000460 ENSG00000000938 
    ##           97.25            0.75

``` r
## treated means
treated.mean <- rowMeans(treated)
head(treated.mean)
```

    ## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 
    ##          658.00            0.00          546.00          316.50 
    ## ENSG00000000460 ENSG00000000938 
    ##           78.75            0.00

``` r
## combine treated and control means
meancounts <- data.frame(control.mean, treated.mean)
colSums(meancounts)
```

    ## control.mean treated.mean 
    ##     23005324     22196524

``` r
## plot meancounts
plot(meancounts$control.mean, meancounts$treated.mean, log = "xy")
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0
    ## omitted from logarithmic plot

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0
    ## omitted from logarithmic plot

![](Class14_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
##log-log plot (adjust axis)
meancounts$log2fc <- log2(meancounts[,"treated.mean"]/meancounts[,"control.mean"])
head(meancounts)
```

    ##                 control.mean treated.mean      log2fc
    ## ENSG00000000003       900.75       658.00 -0.45303916
    ## ENSG00000000005         0.00         0.00         NaN
    ## ENSG00000000419       520.50       546.00  0.06900279
    ## ENSG00000000457       339.75       316.50 -0.10226805
    ## ENSG00000000460        97.25        78.75 -0.30441833
    ## ENSG00000000938         0.75         0.00        -Inf

``` r
#we need to filter ou the zero value entry genes
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)
head(zero.vals)
```

    ##                 row col
    ## ENSG00000000005   2   1
    ## ENSG00000004848  65   1
    ## ENSG00000004948  70   1
    ## ENSG00000005001  73   1
    ## ENSG00000006059 121   1
    ## ENSG00000006071 123   1

``` r
to.rm <- unique(zero.vals[,1])
#unique -> zeros don't get counted twice
mycounts <- meancounts[-to.rm,]
##-to.rm -> everything but meancount[-to.rm]
head(mycounts)
```

    ##                 control.mean treated.mean      log2fc
    ## ENSG00000000003       900.75       658.00 -0.45303916
    ## ENSG00000000419       520.50       546.00  0.06900279
    ## ENSG00000000457       339.75       316.50 -0.10226805
    ## ENSG00000000460        97.25        78.75 -0.30441833
    ## ENSG00000000971      5219.00      6687.50  0.35769358
    ## ENSG00000001036      2327.00      1785.75 -0.38194109

``` r
#how many genes above a designated threshold (in this case, log2fc +-2)
up.ind <- mycounts$log2fc > 2
down.ind <- mycounts$log2fc < (-2)
sum(up.ind)
```

    ## [1] 250

``` r
sum(down.ind)
```

    ## [1] 367

``` r
anno <- read.csv("annotables_grch38.csv")
```

``` r
head(anno)
```

    ##           ensgene entrez   symbol chr     start       end strand
    ## 1 ENSG00000000003   7105   TSPAN6   X 100627109 100639991     -1
    ## 2 ENSG00000000005  64102     TNMD   X 100584802 100599885      1
    ## 3 ENSG00000000419   8813     DPM1  20  50934867  50958555     -1
    ## 4 ENSG00000000457  57147    SCYL3   1 169849631 169894267     -1
    ## 5 ENSG00000000460  55732 C1orf112   1 169662007 169854080      1
    ## 6 ENSG00000000938   2268      FGR   1  27612064  27635277     -1
    ##          biotype
    ## 1 protein_coding
    ## 2 protein_coding
    ## 3 protein_coding
    ## 4 protein_coding
    ## 5 protein_coding
    ## 6 protein_coding
    ##                                                                                                  description
    ## 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
    ## 2                                                            tenomodulin [Source:HGNC Symbol;Acc:HGNC:17757]
    ## 3 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
    ## 4                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
    ## 5                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
    ## 6                          FGR proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:3697]

``` r
dim(anno)
```

    ## [1] 66531     9

``` r
head(mycounts)
```

    ##                 control.mean treated.mean      log2fc
    ## ENSG00000000003       900.75       658.00 -0.45303916
    ## ENSG00000000419       520.50       546.00  0.06900279
    ## ENSG00000000457       339.75       316.50 -0.10226805
    ## ENSG00000000460        97.25        78.75 -0.30441833
    ## ENSG00000000971      5219.00      6687.50  0.35769358
    ## ENSG00000001036      2327.00      1785.75 -0.38194109

``` r
##using the merge() function to add the annotation data to our mycounts results

mycounts.anno <- merge(x=mycounts, y=anno, by.x="row.names", by.y = "ensgene")
head(mycounts.anno)
```

    ##         Row.names control.mean treated.mean      log2fc entrez   symbol
    ## 1 ENSG00000000003       900.75       658.00 -0.45303916   7105   TSPAN6
    ## 2 ENSG00000000419       520.50       546.00  0.06900279   8813     DPM1
    ## 3 ENSG00000000457       339.75       316.50 -0.10226805  57147    SCYL3
    ## 4 ENSG00000000460        97.25        78.75 -0.30441833  55732 C1orf112
    ## 5 ENSG00000000971      5219.00      6687.50  0.35769358   3075      CFH
    ## 6 ENSG00000001036      2327.00      1785.75 -0.38194109   2519    FUCA2
    ##   chr     start       end strand        biotype
    ## 1   X 100627109 100639991     -1 protein_coding
    ## 2  20  50934867  50958555     -1 protein_coding
    ## 3   1 169849631 169894267     -1 protein_coding
    ## 4   1 169662007 169854080      1 protein_coding
    ## 5   1 196651878 196747504      1 protein_coding
    ## 6   6 143494811 143511690     -1 protein_coding
    ##                                                                                                  description
    ## 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
    ## 2 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
    ## 3                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
    ## 4                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
    ## 5                                                     complement factor H [Source:HGNC Symbol;Acc:HGNC:4883]
    ## 6                                          fucosidase, alpha-L- 2, plasma [Source:HGNC Symbol;Acc:HGNC:4008]

``` r
##nsource("http://bioconductor.org/biocLite.R")
biocLite("AnnotationDbi")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.3 (2017-03-06).

    ## Installing package(s) 'AnnotationDbi'

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/ww/dtnfg13n74n9q1bnfgwys6540000gn/T//Rtmpp1WNZE/downloaded_packages

    ## Old packages: 'ade4', 'backports', 'boot', 'checkmate', 'class',
    ##   'cluster', 'codetools', 'colorspace', 'data.table', 'evaluate',
    ##   'foreign', 'Hmisc', 'lattice', 'MASS', 'Matrix', 'mgcv',
    ##   'RcppArmadillo', 'RCurl', 'rpart', 'scales', 'stringi', 'stringr',
    ##   'survival', 'tibble', 'utf8', 'xfun', 'XML'

``` r
biocLite("org.Hs.eg.db")
```

    ## BioC_mirror: https://bioconductor.org

    ## Using Bioconductor 3.4 (BiocInstaller 1.24.0), R 3.3.3 (2017-03-06).

    ## Installing package(s) 'org.Hs.eg.db'

    ## installing the source package 'org.Hs.eg.db'

    ## Old packages: 'ade4', 'backports', 'boot', 'checkmate', 'class',
    ##   'cluster', 'codetools', 'colorspace', 'data.table', 'evaluate',
    ##   'foreign', 'Hmisc', 'lattice', 'MASS', 'Matrix', 'mgcv',
    ##   'RcppArmadillo', 'RCurl', 'rpart', 'scales', 'stringi', 'stringr',
    ##   'survival', 'tibble', 'utf8', 'xfun', 'XML'

``` r
library("AnnotationDbi")
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colnames,
    ##     do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff,
    ##     sort, table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     colMeans, colSums, expand.grid, rowMeans, rowSums

``` r
library("org.Hs.eg.db")
```

    ## 

``` r
columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
    ##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
    ## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
    ## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [25] "UNIGENE"      "UNIPROT"

``` r
mycounts$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(mycounts),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
head(mycounts)
```

    ##                 control.mean treated.mean      log2fc   symbol
    ## ENSG00000000003       900.75       658.00 -0.45303916   TSPAN6
    ## ENSG00000000419       520.50       546.00  0.06900279     DPM1
    ## ENSG00000000457       339.75       316.50 -0.10226805    SCYL3
    ## ENSG00000000460        97.25        78.75 -0.30441833 C1orf112
    ## ENSG00000000971      5219.00      6687.50  0.35769358      CFH
    ## ENSG00000001036      2327.00      1785.75 -0.38194109    FUCA2

``` r
library(DESeq2)
```

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

``` r
citation("DESeq2")
```

    ## 
    ##   Michael I Love, Wolfgang Huber and Simon Anders (2014):
    ##   Moderated estimation of fold change and dispersion for RNA-Seq
    ##   data with DESeq2. Genome Biology
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
    ##     author = {Michael I Love and Wolfgang Huber and Simon Anders},
    ##     year = {2014},
    ##     journal = {Genome Biology},
    ##     doi = {10.1186/s13059-014-0550-8},
    ##     volume = {15},
    ##     issue = {12},
    ##     pages = {550},
    ##   }

``` r
## build or DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
dds
```

    ## class: DESeqDataSet 
    ## dim: 38694 8 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(38694): ENSG00000000003 ENSG00000000005 ...
    ##   ENSG00000283120 ENSG00000283123
    ## rowData names(0):
    ## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    ## colData names(4): id dex celltype geo_id

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
## get our results
res <- results(dds)
summary(res, alpha=0.01)
```

    ## 
    ## out of 25258 with nonzero total read count
    ## adjusted p-value < 0.01
    ## LFC > 0 (up)     : 849, 3.4% 
    ## LFC < 0 (down)   : 584, 2.3% 
    ## outliers [1]     : 212, 0.84% 
    ## low counts [2]   : 9918, 39% 
    ## (mean count < 10)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
plot(res$log2FoldChange, -log(res$padj))
```

![](Class14_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
mycols <- rep("gray", nrow(res))
mycols[abs(res$log2FoldChange) > 2] = "red"
subset.inds <- res$padj < 0.01 & abs(res$log2FoldChange) > 2
mycols[subset.inds] <- "blue"

plot(res$log2FoldChange, -log(res$padj), col=mycols)
```

![](Class14_files/figure-markdown_github/unnamed-chunk-28-1.png)
