## Introduction

In this section, identifier (IDs) mapping is performed on an example
transcriptomics data set, which was original annotated using HGNC
symbols. The dataset has been preprocessed already, for details see step
1 and 2 of the multi-omics workflow at:
<https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/transcriptomics_analysis>
. We map the HGNC symbols to Entrez Gene and Ensembl IDs, since tools
downstream of this step require different input formats for the IDs.

We use two tools for this mapping; first we use org.Hs.eg.db
\[<doi:10.18129/B9.bioc.org.Hs.eg.db>\] and AnnotationDbi; second we use
BridgeDb \[<doi:10.18129/B9.bioc.BridgeDbR>\].

## R environment setup

``` r
# empty the R environment
rm (list = ls())
# check if libraries are already installed, otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"org.Hs.eg.db" %in% installed.packages()) BiocManager::install("org.Hs.eg.db")  
if(!"AnnotationDbi" %in% installed.packages()) BiocManager::install("AnnotationDbi")
if(!"BridgeDbR" %in% installed.packages()) BiocManager::install("BridgeDbR")
if(!"dplyr" %in% installed.packages()) install.packages("dplyr")
if(!"rmarkdown" %in% installed.packages())install.packages("rmarkdown") 
if(!"data.table" %in% installed.packages())install.packages("data.table")
if(!"knitr" %in% installed.packages())install.packages("knitr")

 

#load installed libraries
suppressPackageStartupMessages({
  library(rstudioapi) # interface for interacting with RStudio IDE with R code.
  library(org.Hs.eg.db) #This is the organism annotation package ("org") for Homo sapiens ("Hs"), organized as an AnnotationDbi   package ("db"), using Entrez Gene IDs ("eg") as primary key.
  library(AnnotationDbi) # for connecting and querying annotation databases
  library(BridgeDbR)
  library(dplyr)
  library(rmarkdown)
  library(data.table)
  library(knitr)
})

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset

The data will be read for the disease on two biopsy locations

``` r
#we have four datasets, two different disorders and two different biopsy locations
#Reading the files for CD disorder
(filenames <- list.files("data", pattern = "CD", full.names = TRUE)) #Printing the file names
```

    ## [1] "data/table_CD_Ileum_vs_nonIBD_Ileum.tab"  
    ## [2] "data/table_CD_Rectum_vs_nonIBD_Rectum.tab"

``` r
dataset_CD <- lapply(filenames, read.delim) #Reading the files

#filter out  unused columns, we select geneSymbol, log2FC and pvalue and
#merge two dataset of two locations into one datafile
dataset_CD <- merge(dataset_CD[[1]] %>% select("X", "FoldChange", "padj"), #Merging and sub-setting the columns required
                    dataset_CD[[2]] %>% select("X", "FoldChange", "padj"), 
                    by = "X", all = TRUE)

#Reading the files for UC disorder
(filenames <- list.files("data", pattern = "UC", full.names = TRUE)) #Printing the file names
```

    ## [1] "data/table_UC_Ileum_vs_nonIBD_Ileum.tab"  
    ## [2] "data/table_UC_Rectum_vs_nonIBD_Rectum.tab"

``` r
dataset_UC <- lapply(filenames, read.delim) #Reading the files

#filter out  unused columns, we select geneSymbol, log2FC and pvalue and
#merge two dataset of two locations into one datafile per disorder
dataset_UC <- merge(dataset_UC[[1]] %>% select("X", "FoldChange", "padj"), #Merging and sub-setting the columns required
                    dataset_UC[[2]] %>% select("X", "FoldChange", "padj"), 
                    by = "X", all = TRUE)

#change column names
colnames(dataset_CD) <- colnames(dataset_UC) <- c("GeneSymbol", "log2FC_ileum", "pvalue_ileum", "log2FC_rectum", "pvalue_rectum")

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD"))) # removing variables that are not required
```

## Converting hgnc gene symbols to the corresponding Entrez (NCBI) gene IDs (org.Hs.eg.db)

``` r
#converting gene symbols to entrez ID since these are required for the enrichR function
hs <- org.Hs.eg.db #This object is a simple mapping of Entrez Gene identifier

# check if all the gene symbols are the same in both datasets
all(dataset_CD$GeneSymbol == dataset_UC$GeneSymbol) 
```

    ## [1] TRUE

``` r
##same gene symbols, so mapping based on one of the files 
entrezID <- AnnotationDbi::select(hs, keys = dataset_CD$GeneSymbol, 
            columns = c("ENTREZID", "SYMBOL"), 
            keytype = "SYMBOL")

# checking the one-to-multiple mappings
all(table(entrezID$SYMBOL) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
entrezID_doubles_Hs <- length(table(entrezID$SYMBOL) [table(entrezID$SYMBOL) > 1])
# entrezID_doubles_Hs <- names(table(entrezID$SYMBOL)[table(entrezID$SYMBOL) > 1])
# entrezID [entrezID$SYMBOL %in% entrezID_doubles_Hs, ] # If you want to see which genes have multiple Entrez Gene identifier

#filter out double gene symbols
entrezID <- entrezID %>% distinct(entrezID$SYMBOL, .keep_all = TRUE)

# add entrezIDs for each gene symbol in the dataset
dataset_CD$ENTREZ.ID_org.Hs <- entrezID$ENTREZID [match(dataset_CD$GeneSymbol, entrezID$SYMBOL)] 
dataset_UC$ENTREZ.ID_org.Hs <- entrezID$ENTREZID [match(dataset_UC$GeneSymbol, entrezID$SYMBOL)] 

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "hs"))) # removing variables that are not required
```

## Converting hgnc gene symbols to the corresponding Ensembl IDs (org.Hs.eg.db)

``` r
#converting gene symbols to Ensembl ID since these are required for the Cytoscape multiomics visualization
ensemblID <- AnnotationDbi::select(hs, keys = dataset_CD$GeneSymbol, 
            columns = c("ENSEMBL", "SYMBOL"), 
            keytype = "SYMBOL")

# checking the one-to-multiple mappings
all(table(ensemblID$SYMBOL) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
ensemblID_doubles_Hs <- length(table(ensemblID$SYMBOL)[table(ensemblID$SYMBOL) > 1])
# length(na.omit(ensemblID$SYMBOL)) - length(unique(na.omit(ensemblID$SYMBOL)))
# ensemblID_doubles_Hs <- names(table(ensemblID$SYMBOL)[table(ensemblID$SYMBOL) > 1])
# ensemblID %>% filter(SYMBOL %in% ensemblID_doubles_Hs) %>% arrange(SYMBOL) %>% head() # If you want to see which genes have multiple Ensembl ID

#filter out double gene symbols
ensemblID <- ensemblID %>% distinct(ensemblID$SYMBOL, .keep_all = TRUE)
# add entrezIDs for each gene symbol in the dataset
dataset_CD$Ensembl.ID_org.Hs <- ensemblID$ENSEMBL [match(dataset_CD$GeneSymbol, ensemblID$SYMBOL)] 
dataset_UC$Ensembl.ID_org.Hs <- ensemblID$ENSEMBL [match(dataset_UC$GeneSymbol, ensemblID$SYMBOL)] 

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs"))) # removing variables that are not required
```

## Converting hgnc gene symbols to the corresponding Entrez (NCBI) gene IDs (BridgeDb)

please look up
[here](https://bridgedb.github.io/pages/system-codes.html) for the
correct system code for each datasource

``` r
#loading human derby database and converting gene symbols to entrez ID since these are required for the enrichR function
dbLocation = "data/Hs_Derby_Ensembl_105.bridge"
mapper <- loadDatabase(dbLocation)
input <- data.frame(source = rep("H", length(dataset_CD$GeneSymbol)),
                     identifier = dataset_CD$GeneSymbol)
entrezID <- maps(mapper = mapper, input, target = "L") 

# checking the one-to-multiple mappings
all(table(entrezID$identifier) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
entrezID_doubles_BridgeDb <- length(table(entrezID$identifier) [table(entrezID$identifier) > 1])
# entrezID_doubles_BridgeDb <- names(table(entrezID$identifier)[table(entrezID$identifier) > 1])
# entrezID [entrezID$identifier %in% entrezID_doubles_BridgeDb, ] # If you want to see which genes have multiple Entrez Gene identifier
#filter out double gene symbols
entrezID <- entrezID %>% distinct(entrezID$identifier, .keep_all = TRUE)

# add entrezIDs for each gene symbol in the dataset
dataset_CD$ENTREZ.ID_BridgeDb <- entrezID$mapping [match(dataset_CD$GeneSymbol, entrezID$identifier)] 
dataset_UC$ENTREZ.ID_BridgeDb <- entrezID$mapping [match(dataset_UC$GeneSymbol, entrezID$identifier)] 

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs", "mapper", "input", "entrezID_doubles_BridgeDb"))) # removing variables that are not required
```

## Converting hgnc gene symbols to the corresponding Ensembl IDs (BridgeDb)

``` r
#converting gene symbols to Ensembl ID since these are required for the Cytoscape multiomics visualization
ensemblID <- maps(mapper = mapper, input, target = "En") 

# checking the one-to-multiple mappings
all(table(ensemblID$identifier) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
ensemblID_doubles_BridgeDb <- length(table(ensemblID$identifier) [table(ensemblID$identifier) > 1])
# ensemblID_doubles_BridgeDb <- names(table(ensemblID$identifier)[table(ensemblID$identifier) > 1])
# ensemblID [ensemblID$identifier %in% ensemblID_doubles_BridgeDb, ] # If you want to see which genes have multiple Ensembl ID
#filter out double gene symbols
ensemblID <- ensemblID %>% distinct(ensemblID$identifier, .keep_all = TRUE)

# add ensemblIDs for each gene symbol in the dataset
dataset_CD$Ensembl.ID_BridgeDb <- ensemblID$mapping [match(dataset_CD$GeneSymbol, ensemblID$identifier)] 
dataset_UC$Ensembl.ID_BridgeDb <- ensemblID$mapping [match(dataset_UC$GeneSymbol, ensemblID$identifier)] 

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs", "mapper", "input", "entrezID_doubles_BridgeDb", "ensemblID_doubles_BridgeDb"))) # removing variables that are not required
```

## Using BridgeDb for secondary to primary mapping of hgnc gene symbols

### a. mapping the hgnc symbols to hgnc IDs

HGNC recycles the hgnc gene symbols which means there are identifiers
(symbols) that are secondary for a certain entity (e.g the primary hgnc
symbol is FYN and the secondary hgnc symbol is SLK), while the same
identifiers might be used as a primary symbol for another entity (SLK is
a primary symbol for another gene). This makes the mapping based on hgnc
gene symbols challenging.
`Here We assumed when a hgnc gene symbol has a hgnc ID, it is primary and we only map those without a hgnc ID.`

``` r
#converting gene symbols to hgnc ID 
hgncID <- maps(mapper = mapper, input, target = "Hac") 

# checking the one-to-multiple mappings
all(table(hgncID$identifier) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] TRUE

``` r
# add HGNC id for each gene symbol in the dataset
dataset_CD$HGNC.ID_BridgeDb <- hgncID$mapping[match(dataset_CD$GeneSymbol, hgncID$identifier)]
dataset_UC$HGNC.ID_BridgeDb <- hgncID$mapping[match(dataset_UC$GeneSymbol, hgncID$identifier)]

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs", "entrezID_doubles_BridgeDb", "ensemblID_doubles_BridgeDb"))) # removing variables that are not required
```

### b. mapping secondary hgnc symbols to primary hgnc symbols

Mapping hgnc symbols without a HGNC id mapped

``` r
# loading the HGNC secondary id derby database
dbLocation = "derbyDB/hgnc_all_secIds.bridge"
mapper <- loadDatabase(dbLocation)
#Subset hgnc gene symbols with no hgnc ID
input <- dataset_CD$GeneSymbol [is.na(dataset_CD$HGNC.ID_BridgeDb)] 
input <- data.frame(source = rep("H", length(input)),
                    identifier = input)
#converting secondary gene symbols to primary gene symbols 
hgnc <- maps(mapper, input, "H") %>% 
    filter(isPrimary == "T") # Keeping only rows where the mapping is annotated as primary id (defined in BridgeDb java  library)

# checking the one-to-multiple mappings
all(table(hgnc$identifier) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
hgnc_doubles_PriID_BridgeDb <- length(table(hgnc$identifier) [table(hgnc$identifier) > 1])
# hgnc_doubles_PriID_BridgeDb <- names(table(hgnc$identifier)[table(hgnc$identifier) > 1])
# hgnc [hgnc$identifier %in% hgnc_doubles_PriID_BridgeDb, ] # If you want to see which genes have multiple primary id
```

Checking the reason that the secondary symbol mapped to multiple primary
symbols

``` r
hgnc [hgnc$identifier == "AGPAT9", ]
```

    ##            source identifier target mapping isPrimary
    ## H:GPAT3:T       H     AGPAT9      H   GPAT3         T
    ## H:LPCAT1:T      H     AGPAT9      H  LPCAT1         T

``` r
#            source identifier target mapping isPrimary
# H:GPAT3:T       H     AGPAT9      H   GPAT3         T
# H:LPCAT1:T      H     AGPAT9      H  LPCAT1         T
```

In this example where we manually checked: GPAT3 is the Previous symbol
(this field displays any symbols that were previously HGNC-approved
nomenclature.) LPCAT1 is the Alias symbol (Alternative symbols that have
been used to refer to the gene. Aliases may be from literature, from
other databases or may be added to represent membership of a gene
group.)

``` r
hgnc [hgnc$identifier == "U3", ]
```

    ##              source identifier target  mapping isPrimary
    ## H:SNORD3P3:T      H         U3      H SNORD3P3         T
    ## H:SNORD3A:T       H         U3      H  SNORD3A         T
    ## H:SNORD3P4:T      H         U3      H SNORD3P4         T
    ## H:SNORD3P1:T      H         U3      H SNORD3P1         T
    ## H:SNORD3F:T       H         U3      H  SNORD3F         T

This is an example of using the same alias for 5 different hgnc genes

The issue of one-to-multiple mappings could not be easily fixed (needs
manual check) and emphasizes on the importance of using unique,
persistent, resolvable identifiers.

We then checked whether the primary symbols mapped are already present
in the gene expression data

``` r
all(!hgnc$mapping %in% dataset_CD$GeneSymbol)#True if the mapped id doesn't annotate another gene and is not present in the dataset, otherwise FALSE
```

    ## [1] FALSE

Some of the primary symbols are already presented in the dataset (used
as primary symbol for another gene). This shows another limitation of
using gene symbols as identifiers, and to avoid the duplicated gene
symbols we only mapped those symbols that are not presented in the
dataset.

``` r
hgnc <- hgnc [!hgnc$mapping %in% dataset_CD$GeneSymbol,] 
#filter out double gene symbols
hgnc <- hgnc %>% distinct(hgnc$identifier, .keep_all = TRUE)

# add primary (current) hgnc gene symbol each gene symbol in the dataset
dataset_CD$Current_GeneSymbol <- hgnc$mapping[match(dataset_CD$GeneSymbol, hgnc$identifier)]
dataset_CD$Current_GeneSymbol [is.na(dataset_CD$Current_GeneSymbol)] = dataset_CD$GeneSymbol [is.na(dataset_CD$Current_GeneSymbol)]

dataset_UC$Current_GeneSymbol <- hgnc$mapping[match(dataset_UC$GeneSymbol, hgnc$identifier)]
dataset_UC$Current_GeneSymbol [is.na(dataset_CD$Current_GeneSymbol)] = dataset_UC$GeneSymbol [is.na(dataset_CD$Current_GeneSymbol)]

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs", "entrezID_doubles_BridgeDb", "ensemblID_doubles_BridgeDb", "hgnc_doubles_PriID_BridgeDb"))) # removing variables that are not required
```

## Converting `primary` hgnc gene symbols to the corresponding Entrez (NCBI) gene IDs (BridgeDb)

please look up
[here](https://bridgedb.github.io/pages/system-codes.html) for the
correct system code for each datasource

``` r
#loading human derby database and converting primary gene symbols to entrez ID since these are required for the enrichR function
dbLocation = "data/Hs_Derby_Ensembl_105.bridge"
mapper <- loadDatabase(dbLocation)
input <- data.frame(source = rep("H", length(dataset_CD$Current_GeneSymbol)),
                     identifier = dataset_CD$Current_GeneSymbol)
entrezID <- maps(mapper = mapper, input, target = "L") 

# checking the one-to-multiple mappings
all(table(entrezID$identifier) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
entrezID_doubles_PriID_BridgeDb <- length(table(entrezID$identifier) [table(entrezID$identifier) > 1])
# entrezID_doubles_PriID_BridgeDb <- names(table(entrezID$identifier)[table(entrezID$identifier) > 1])
# entrezID [entrezID$identifier %in% entrezID_doubles_PriID_BridgeDb, ] # If you want to see which genes have multiple Entrez Gene identifier
#filter out double gene symbols
entrezID <- entrezID %>% distinct(entrezID$identifier, .keep_all = TRUE)

# add entrezIDs for each gene symbol in the dataset
dataset_CD$ENTREZ.ID_BridgeDb_PriID <- entrezID$mapping [match(dataset_CD$GeneSymbol, entrezID$identifier)] 
dataset_UC$ENTREZ.ID_BridgeDb_PriID <- entrezID$mapping [match(dataset_UC$GeneSymbol, entrezID$identifier)] 

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs", "entrezID_doubles_BridgeDb", "ensemblID_doubles_BridgeDb", "hgnc_doubles_PriID_BridgeDb", "entrezID_doubles_PriID_BridgeDb", "input", "mapper"))) # removing variables that are not required
```

## Converting `primary` hgnc gene symbols to the corresponding Ensembl IDs (BridgeDb)

``` r
#converting gene symbols to Ensembl ID since these are required for the Cytoscape multiomics visualization
ensemblID <- maps(mapper = mapper, input, target = "En") 

# checking the one-to-multiple mappings
all(table(ensemblID$identifier) == 1) #True if there is no one-to-multiple mapping, otherwise FALSE
```

    ## [1] FALSE

``` r
#store one-to-multiple mapping info
ensemblID_doubles_PriID_BridgeDb <- length(table(ensemblID$identifier) [table(ensemblID$identifier) > 1])
# ensemblID_doubles_PriID_BridgeDb <- names(table(ensemblID$identifier)[table(ensemblID$identifier) > 1])
# ensemblID [ensemblID$identifier %in% ensemblID_doubles_PriID_BridgeDb, ] # If you want to see which genes have multiple Ensembl ID
#filter out double gene symbols
ensemblID <- ensemblID %>% distinct(ensemblID$identifier, .keep_all = TRUE)

# add ensemblIDs for each gene symbol in the dataset
dataset_CD$Ensembl.ID_BridgeDb_PriID <- ensemblID$mapping [match(dataset_CD$GeneSymbol, ensemblID$identifier)] 
dataset_UC$Ensembl.ID_BridgeDb_PriID <- ensemblID$mapping [match(dataset_UC$GeneSymbol, ensemblID$identifier)] 

rm(list = setdiff(ls(), c("dataset_UC", "dataset_CD", "entrezID_doubles_Hs", "ensemblID_doubles_Hs", "entrezID_doubles_BridgeDb", "ensemblID_doubles_BridgeDb", "hgnc_doubles_PriID_BridgeDb", "entrezID_doubles_PriID_BridgeDb", "ensemblID_doubles_PriID_BridgeDb"))) # removing variables that are not required
```

##Mapping stats:

    ## [1] "The total number of HGNC Symbol in the transcriptomics dataset is: 17670"

|                                                                     | org.Hs | BridgeDb | PrimaryID_BridgeDb |
|:---------------------------------------------|-----:|------:|-------------:|
| The total number of unique Entrez IDs                               |  15022 |    14133 |                  0 |
| The total number of missing mappings for HGNC Symbol to Entrez IDs  |   2648 |     3537 |                  0 |
| The total number of one-to-many mappings for Entrez IDs             |      3 |       76 |                 79 |
| The total number of unique Ensembl IDs                              |  14596 |    15009 |              15009 |
| The total number of missing mappings for HGNC Symbol to Ensembl IDs |   3074 |     2661 |               2661 |
| The total number of one-to-many mappings for Ensembl IDs            |    835 |        6 |                  7 |

    ## 
    ## To cite package 'org.Hs.eg.db' in publications use:
    ## 
    ##   Marc Carlson (2021). org.Hs.eg.db: Genome wide annotation for Human.
    ##   R package version 3.14.0.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {org.Hs.eg.db: Genome wide annotation for Human},
    ##     author = {Marc Carlson},
    ##     year = {2021},
    ##     note = {R package version 3.14.0},
    ##   }
    ## 
    ## ATTENTION: This citation information has been auto-generated from the
    ## package DESCRIPTION file and may need manual editing, see
    ## 'help("citation")'.

    ## 
    ## To cite package 'BridgeDbR' in publications use:
    ## 
    ##   Chris Leemans, Egon Willighagen, Denise Slenter, Anwesha Bohler, Lars
    ##   Eijssen and Tooba Abbassi-Daloii (2022). BridgeDbR: Code for using
    ##   BridgeDb identifier mapping framework from within R. R package
    ##   version 2.7.2. https://github.com/bridgedb/BridgeDbR
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {BridgeDbR: Code for using BridgeDb identifier mapping framework from within
    ## R},
    ##     author = {Chris Leemans and Egon Willighagen and Denise Slenter and Anwesha Bohler and Lars Eijssen and Tooba Abbassi-Daloii},
    ##     year = {2022},
    ##     note = {R package version 2.7.2},
    ##     url = {https://github.com/bridgedb/BridgeDbR},
    ##   }

##Save data, print session info and remove large datasets:

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22000)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.39           data.table_1.14.2    rmarkdown_2.14      
    ##  [4] dplyr_1.0.9          BridgeDbR_2.7.2      rJava_1.0-6         
    ##  [7] org.Hs.eg.db_3.14.0  AnnotationDbi_1.56.2 IRanges_2.28.0      
    ## [10] S4Vectors_0.32.4     Biobase_2.54.0       BiocGenerics_0.40.0 
    ## [13] rstudioapi_0.13     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] KEGGREST_1.34.0        tidyselect_1.1.2       xfun_0.31             
    ##  [4] purrr_0.3.4            vctrs_0.4.1            generics_0.1.2        
    ##  [7] htmltools_0.5.2        yaml_2.3.5             utf8_1.2.2            
    ## [10] blob_1.2.3             rlang_1.0.2            pillar_1.7.0          
    ## [13] glue_1.6.2             DBI_1.1.2              bit64_4.0.5           
    ## [16] GenomeInfoDbData_1.2.7 lifecycle_1.0.1        stringr_1.4.0         
    ## [19] zlibbioc_1.40.0        Biostrings_2.62.0      memoise_2.0.1         
    ## [22] evaluate_0.15          fastmap_1.1.0          GenomeInfoDb_1.30.1   
    ## [25] curl_4.3.2             fansi_1.0.3            highr_0.9             
    ## [28] Rcpp_1.0.8.3           BiocManager_1.30.18    cachem_1.0.6          
    ## [31] XVector_0.34.0         bit_4.0.4              png_0.1-7             
    ## [34] digest_0.6.29          stringi_1.7.6          cli_3.2.0             
    ## [37] tools_4.1.2            bitops_1.0-7           magrittr_2.0.3        
    ## [40] RCurl_1.98-1.6         RSQLite_2.2.14         tibble_3.1.7          
    ## [43] crayon_1.5.1           pkgconfig_2.0.3        ellipsis_0.3.2        
    ## [46] assertthat_0.2.1       httr_1.4.3             R6_2.5.1              
    ## [49] compiler_4.1.2
