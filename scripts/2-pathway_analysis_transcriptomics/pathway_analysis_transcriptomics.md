## Introduction

In this section, we will map the significantly changed transcripts to
pathways from the WikiPathways and Reactome collection, in order to
investigate the potential difference in statistical calculations
depending on the input identifiers (IDs). We will work with three IDs:
HGNC symbols, Entrez (NCBI) gene IDs and Ensembl. The first was present
in the original dataset, while the second and third are generated with
different mapping tools (org.Hs and BridgeDbR). We will use the SPARQL
endpoint interface to compare our dataset to the content of both pathway
databases. Other approaches exist (e.g. gProfiler), however not all of
these approaches allow performing their analysis with different
identifiers.

## R environment setup

``` r
# empty the R environment
rm (list = ls())

# check if libraries are already installed > otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"dplyr" %in% installed.packages()) install.packages("dplyr")
if(!"data.table" %in% installed.packages()) install.packages("data.table")
if(!"knitr" %in% installed.packages()) install.packages("knitr")

#loading installed libraries
suppressPackageStartupMessages({
  library(rstudioapi) # interface for interacting with RStudio IDE with R code.
  library(dplyr)
  library(data.table)
  library(knitr)
  })

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset and subsetting genes with significant changes

The data from step 1 will be read in and further processed. For a fairer
comparison, we report how many rows per ID type contain a significantly
changed gene.

``` r
#we have four datasets, two different disorders and two different biopsy locations
dataset_CD <- read.delim("../1-identifier_mapping_transcriptomics/results/IDMapping_CD")
dataset_UC <- read.delim("../1-identifier_mapping_transcriptomics/results/IDMapping_UC")

#The following selection criteria will be used to determine the significantly changed transcripts Fold change = 1.5, log2FC = 0.58 and p.value < 0.05.
#CD dataset
## Significant rows
sig.rows.ileum_CD <- which((dataset_CD$log2FC_ileum >= 0.58 | dataset_CD$log2FC_ileum <= -0.58) & dataset_CD$pvalue_ileum < 0.05) #ileum location
sig.rows.rectum_CD <- which((dataset_CD$log2FC_rectum >= 0.58 | dataset_CD$log2FC_rectum <= -0.58) & dataset_CD$pvalue_rectum < 0.05) #rectum location
##HGNC_Symbols:
sig.genes.ileum_HGNC_CD <- dataset_CD$GeneSymbol[sig.rows.ileum_CD] #ileum location
sig.genes.rectum_HGNC_CD <- dataset_CD$GeneSymbol[sig.rows.rectum_CD] #rectum location
## org.Hs
### Entrez ID:
sig.genes.ileum_ENTREZ_Hs_CD <- na.omit(dataset_CD$ENTREZ.ID_org.Hs[sig.rows.ileum_CD]) #ileum location
sig.genes.rectum_ENTREZ_Hs_CD <- na.omit(dataset_CD$ENTREZ.ID_org.Hs[sig.rows.rectum_CD]) #rectum location
### Ensembl IDs:
sig.genes.ileum_Ensembl_Hs_CD <- na.omit(dataset_CD$Ensembl.ID_org.Hs[sig.rows.ileum_CD]) #ileum location
sig.genes.rectum_Ensembl_Hs_CD <- na.omit(dataset_CD$Ensembl.ID_org.Hs[sig.rows.rectum_CD]) #rectum location
## BridgeDb
### Entrez ID:
sig.genes.ileum_ENTREZ_BridgeDb_CD <- na.omit(dataset_CD$ENTREZ.ID_BridgeDb[sig.rows.ileum_CD]) #ileum location
sig.genes.rectum_ENTREZ_BridgeDb_CD <- na.omit(dataset_CD$ENTREZ.ID_BridgeDb[sig.rows.rectum_CD]) #rectum location
### Ensembl IDs:
sig.genes.ileum_Ensembl_BridgeDb_CD <- na.omit(dataset_CD$Ensembl.ID_BridgeDb[sig.rows.ileum_CD]) #ileum location
sig.genes.rectum_Ensembl_BridgeDb_CD <- na.omit(dataset_CD$Ensembl.ID_BridgeDb[sig.rows.rectum_CD]) #rectum location
## Primary id mapping,BridgeDb,
### primary HGNC Symbols:
sig.genes.ileum_HGNC_PriID_BridgeDb_CD <- dataset_CD$Current_GeneSymbol[sig.rows.ileum_CD] #ileum location
sig.genes.rectum_HGNC_PriID_BridgeDb_CD <- dataset_CD$Current_GeneSymbol[sig.rows.rectum_CD] #rectum location
### Entrez ID:
sig.genes.ileum_ENTREZ_PriID_BridgeDb_CD <- na.omit(dataset_CD$ENTREZ.ID_BridgeDb_PriID[sig.rows.ileum_CD]) #ileum location
sig.genes.rectum_ENTREZ_PriID_BridgeDb_CD <- na.omit(dataset_CD$ENTREZ.ID_BridgeDb_PriID[sig.rows.rectum_CD]) #rectum location
### Ensembl IDs:
sig.genes.ileum_Ensembl_PriID_BridgeDb_CD <- na.omit(dataset_CD$Ensembl.ID_BridgeDb_PriID[sig.rows.ileum_CD]) #ileum location
sig.genes.rectum_Ensembl_PriID_BridgeDb_CD <- na.omit(dataset_CD$Ensembl.ID_BridgeDb_PriID[sig.rows.rectum_CD]) #rectum location

#UC dataset
## Significant rows
sig.rows.ileum_UC <- which((dataset_UC$log2FC_ileum >= 0.58 | dataset_UC$log2FC_ileum <= -0.58) & dataset_UC$pvalue_ileum < 0.05) #ileum location
sig.rows.rectum_UC <- which((dataset_UC$log2FC_rectum >= 0.58 | dataset_UC$log2FC_rectum <= -0.58) & dataset_UC$pvalue_rectum < 0.05) #rectum location
##HGNC_Symbols:
sig.genes.ileum_HGNC_UC <- dataset_UC$GeneSymbol[sig.rows.ileum_UC] #ileum location
sig.genes.rectum_HGNC_UC <- dataset_UC$GeneSymbol[sig.rows.rectum_UC] #rectum location
## org.Hs
### Entrez ID:
sig.genes.ileum_ENTREZ_Hs_UC <- na.omit(dataset_UC$ENTREZ.ID_org.Hs[sig.rows.ileum_UC]) #ileum location
sig.genes.rectum_ENTREZ_Hs_UC <- na.omit(dataset_UC$ENTREZ.ID_org.Hs[sig.rows.rectum_UC]) #rectum location
### Ensembl IDs:
sig.genes.ileum_Ensembl_Hs_UC <- na.omit(dataset_UC$Ensembl.ID_org.Hs[sig.rows.ileum_UC]) #ileum location
sig.genes.rectum_Ensembl_Hs_UC <- na.omit(dataset_UC$Ensembl.ID_org.Hs[sig.rows.rectum_UC]) #rectum location
## BridgeDb
### Entrez ID:
sig.genes.ileum_ENTREZ_BridgeDb_UC <- na.omit(dataset_UC$ENTREZ.ID_BridgeDb[sig.rows.ileum_UC]) #ileum location
sig.genes.rectum_ENTREZ_BridgeDb_UC <- na.omit(dataset_UC$ENTREZ.ID_BridgeDb[sig.rows.rectum_UC]) #rectum location
### Ensembl IDs:
sig.genes.ileum_Ensembl_BridgeDb_UC <- na.omit(dataset_UC$Ensembl.ID_BridgeDb[sig.rows.ileum_UC]) #ileum location
sig.genes.rectum_Ensembl_BridgeDb_UC <- na.omit(dataset_UC$Ensembl.ID_BridgeDb[sig.rows.rectum_UC]) #rectum location
## Primary id mapping,BridgeDb,
### primary HGNC Symbols:
sig.genes.ileum_HGNC_PriID_BridgeDb_UC <- dataset_UC$Current_GeneSymbol[sig.rows.ileum_UC] #ileum location
sig.genes.rectum_HGNC_PriID_BridgeDb_UC <- dataset_UC$Current_GeneSymbol[sig.rows.rectum_UC] #rectum location
### Entrez ID:
sig.genes.ileum_ENTREZ_PriID_BridgeDb_UC <- na.omit(dataset_UC$ENTREZ.ID_BridgeDb_PriID[sig.rows.ileum_UC]) #ileum location
sig.genes.rectum_ENTREZ_PriID_BridgeDb_UC <- na.omit(dataset_UC$ENTREZ.ID_BridgeDb_PriID[sig.rows.rectum_UC]) #rectum location
### Ensembl IDs:
sig.genes.ileum_Ensembl_PriID_BridgeDb_UC <- na.omit(dataset_UC$Ensembl.ID_BridgeDb_PriID[sig.rows.ileum_UC]) #ileum location
sig.genes.rectum_Ensembl_PriID_BridgeDb_UC <- na.omit(dataset_UC$Ensembl.ID_BridgeDb_PriID[sig.rows.rectum_UC]) #rectum location
```

##Sign. Mapping stats:

``` r
##HGNC_originalData
MappingStats <- data.table(`  ` =  c("The total number of genes in the transcriptomics dataset" ,
                                     "The total number of significant genes with HGNC Symbol IDs for CD (ileum location)",
                                     "The total number of significant genes with HGNC Symbol IDs for CD (rectum location)",
                                     "The total number of significant genes with HGNC Symbol IDs for UC (ileum location)",
                                     "The total number of significant genes with HGNC Symbol IDs for UC (rectum location)",
                                     "The total number of unique Entrez IDs in the transcriptomics dataset" ,
                                     "The total number of significant genes with Entrez IDs for CD (ileum location)",
                                     "The total number of significant genes with Entrez IDs for CD (rectum location)",
                                     "The total number of significant genes with Entrez IDs for UC (ileum location)",
                                     "The total number of significant genes with Entrez IDs for UC (rectum location)",
                                     "The total number of unique Ensembl IDs in the transcriptomics dataset" ,
                                     "The total number of significant genes with Ensembl IDs for CD (ileum location)",
                                     "The total number of significant genes with Ensembl IDs for CD (rectum location)",
                                     "The total number of significant genes with Ensembl IDs for UC (ileum location)",
                                     "The total number of significant genes with Ensembl IDs for UC (rectum location)"),
                           org.Hs = c(nrow (dataset_CD), 
                                      length(sig.genes.ileum_HGNC_CD), length(sig.genes.rectum_HGNC_CD), 
                                      length(sig.genes.ileum_HGNC_UC), length(sig.genes.rectum_HGNC_UC),
                                      nrow (dataset_CD[na.omit(dataset_CD$ENTREZ.ID_org.Hs),]), 
                                      length(sig.genes.ileum_ENTREZ_Hs_CD), length(sig.genes.rectum_ENTREZ_Hs_CD), 
                                      length(sig.genes.ileum_ENTREZ_Hs_UC), length(sig.genes.rectum_ENTREZ_Hs_UC),
                                      nrow (dataset_CD[na.omit(dataset_CD$Ensembl.ID_org.Hs),]), 
                                      length(sig.genes.ileum_Ensembl_Hs_CD), length(sig.genes.rectum_Ensembl_Hs_CD), 
                                      length(sig.genes.ileum_Ensembl_Hs_UC), length(sig.genes.rectum_Ensembl_Hs_UC)),
                           BridgeDb = c(nrow (dataset_CD), 
                                      length(sig.genes.ileum_HGNC_CD), length(sig.genes.rectum_HGNC_CD), 
                                      length(sig.genes.ileum_HGNC_UC),  length(sig.genes.rectum_HGNC_UC),
                                      nrow (dataset_CD[na.omit(dataset_CD$ENTREZ.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_ENTREZ_BridgeDb_CD), length(sig.genes.rectum_ENTREZ_BridgeDb_CD), 
                                      length(sig.genes.ileum_ENTREZ_BridgeDb_UC), length(sig.genes.rectum_ENTREZ_BridgeDb_UC),
                                      nrow (dataset_CD[na.omit(dataset_CD$Ensembl.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_Ensembl_BridgeDb_CD), length(sig.genes.rectum_Ensembl_BridgeDb_CD), 
                                      length(sig.genes.ileum_Ensembl_BridgeDb_UC), length(sig.genes.rectum_Ensembl_BridgeDb_UC)),
                           PrimaryID_BridgeDb = c(c(nrow (dataset_CD), 
                                      length(sig.genes.ileum_HGNC_PriID_BridgeDb_CD), length(sig.genes.rectum_HGNC_PriID_BridgeDb_CD), 
                                      length(sig.genes.ileum_HGNC_PriID_BridgeDb_UC), length(sig.genes.rectum_HGNC_PriID_BridgeDb_UC),
                                      nrow (dataset_CD[na.omit(dataset_CD$ENTREZ.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_ENTREZ_PriID_BridgeDb_CD), length(sig.genes.rectum_ENTREZ_PriID_BridgeDb_CD), 
                                      length(sig.genes.ileum_ENTREZ_PriID_BridgeDb_UC), length(sig.genes.rectum_ENTREZ_PriID_BridgeDb_UC),
                                      nrow (dataset_CD[na.omit(dataset_CD$Ensembl.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_Ensembl_PriID_BridgeDb_CD), length(sig.genes.rectum_Ensembl_PriID_BridgeDb_CD), 
                                      length(sig.genes.ileum_Ensembl_PriID_BridgeDb_UC),  length(sig.genes.rectum_Ensembl_PriID_BridgeDb_UC))))
kable(MappingStats)
```

|                                                                                     | org.Hs | BridgeDb | PrimaryID_BridgeDb |
|:------------------------------------------------|----:|------:|-----------:|
| The total number of genes in the transcriptomics dataset                            |  17670 |    17670 |              17670 |
| The total number of significant genes with HGNC Symbol IDs for CD (ileum location)  |   1896 |     1896 |               1896 |
| The total number of significant genes with HGNC Symbol IDs for CD (rectum location) |   1339 |     1339 |               1339 |
| The total number of significant genes with HGNC Symbol IDs for UC (ileum location)  |    121 |      121 |                121 |
| The total number of significant genes with HGNC Symbol IDs for UC (rectum location) |   6809 |     6809 |               6809 |
| The total number of unique Entrez IDs in the transcriptomics dataset                |  15022 |    14133 |              14133 |
| The total number of significant genes with Entrez IDs for CD (ileum location)       |   1697 |     1558 |               1558 |
| The total number of significant genes with Entrez IDs for CD (rectum location)      |   1195 |     1114 |               1114 |
| The total number of significant genes with Entrez IDs for UC (ileum location)       |     99 |       97 |                 97 |
| The total number of significant genes with Entrez IDs for UC (rectum location)      |   5870 |     5593 |               5593 |
| The total number of unique Ensembl IDs in the transcriptomics dataset               |  14596 |    15009 |              15009 |
| The total number of significant genes with Ensembl IDs for CD (ileum location)      |   1665 |     1695 |               1695 |
| The total number of significant genes with Ensembl IDs for CD (rectum location)     |   1171 |     1195 |               1195 |
| The total number of significant genes with Ensembl IDs for UC (ileum location)      |     97 |       99 |                 99 |
| The total number of significant genes with Ensembl IDs for UC (rectum location)     |   5771 |     5864 |               5864 |

``` r
rm(list = setdiff(ls(), grep ("dataset_UC|dataset_CD|sig.genes", ls(), value = TRUE))) # removing variables that are not required
```

Find pathways for each dataset, based on signifcantly changed genes and
different IDs.

``` r
if(!"SPARQL" %in% installed.packages()) install.packages("SPARQL")
# alternative way to intall this package if you get this warning "package ‘SPARQL’ is not available for this version of R"
if(!"SPARQL" %in% installed.packages()) install.packages("https://cran.r-project.org/src/contrib/Archive/SPARQL/SPARQL_1.16.tar.gz", repos = NULL, type="source")
library(SPARQL)

##Connect to Endpoint WikiPathways
endpointwp <- "https://sparql.wikipathways.org/sparql"
## 1. Query metadata:
queryMetadata <-
"SELECT DISTINCT ?dataset (str(?titleLit) as ?title) ?date ?license 
WHERE {
   ?dataset a void:Dataset ;
   dcterms:title ?titleLit ;
   dcterms:license ?license ;
   pav:createdOn ?date .
 }"
#below code should be performed first to handle the ssl certificate error
options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
resultsMetadata <- SPARQL(endpointwp, queryMetadata, curl_args = list(useragent = R.version.string))
showresultsMetadata <- resultsMetadata$results
remove(queryMetadata, resultsMetadata)

#For now, filter out Reactome PWs due to visualization issues in Cytoscape.
item1 = "PREFIX hgnc: <https://identifiers.org/hgnc.symbol/>
PREFIX entrez: <https://identifiers.org/ncbigene/>
PREFIX ensembl: <https://identifiers.org/ensembl/>
PREFIX cur: <http://vocabularies.wikipathways.org/wp#Curation:>
select distinct ?pathwayRes (str(?wpid) as ?pathway) (str(?title) as ?pathwayTitle) (count(distinct ?geneId) AS ?GenesInPWs) 
where {
VALUES ?geneId {"
item2 = "}
  ?gene dcterms:isPartOf ?pathwayRes ;"
item3_HGNC= "
    wp:bdbHgncSymbol ?geneId ."
item3_ENTREZ= "
     wp:bdbEntrezGene ?geneId ."
item3_Ensembl= "
    wp:bdbEnsembl ?geneId ."
item4= "
  ?pathwayRes a wp:Pathway ;
              wp:organismName 'Homo sapiens' ; 
    dcterms:identifier ?wpid ;
    dc:title ?title .
    
  #?pathwayRes wp:ontologyTag cur:Reactome_Approved . 
  ?pathwayRes wp:ontologyTag cur:AnalysisCollection .   
}

ORDER BY DESC(?GenesInPWs)"

##Split significant genes into list of max. 300 (or x sections if that's easier), otherwise SPARQL endpoint trows a 414 error --> Ensembl has longest ID, 220 entries seems to be the max.
##Merge the content of the split content back together for the output of the PW Analysis.

ls(pattern = "sig.genes")
```

    ##  [1] "sig.genes.ileum_Ensembl_BridgeDb_CD"       
    ##  [2] "sig.genes.ileum_Ensembl_BridgeDb_UC"       
    ##  [3] "sig.genes.ileum_Ensembl_Hs_CD"             
    ##  [4] "sig.genes.ileum_Ensembl_Hs_UC"             
    ##  [5] "sig.genes.ileum_Ensembl_PriID_BridgeDb_CD" 
    ##  [6] "sig.genes.ileum_Ensembl_PriID_BridgeDb_UC" 
    ##  [7] "sig.genes.ileum_ENTREZ_BridgeDb_CD"        
    ##  [8] "sig.genes.ileum_ENTREZ_BridgeDb_UC"        
    ##  [9] "sig.genes.ileum_ENTREZ_Hs_CD"              
    ## [10] "sig.genes.ileum_ENTREZ_Hs_UC"              
    ## [11] "sig.genes.ileum_ENTREZ_PriID_BridgeDb_CD"  
    ## [12] "sig.genes.ileum_ENTREZ_PriID_BridgeDb_UC"  
    ## [13] "sig.genes.ileum_HGNC_CD"                   
    ## [14] "sig.genes.ileum_HGNC_PriID_BridgeDb_CD"    
    ## [15] "sig.genes.ileum_HGNC_PriID_BridgeDb_UC"    
    ## [16] "sig.genes.ileum_HGNC_UC"                   
    ## [17] "sig.genes.rectum_Ensembl_BridgeDb_CD"      
    ## [18] "sig.genes.rectum_Ensembl_BridgeDb_UC"      
    ## [19] "sig.genes.rectum_Ensembl_Hs_CD"            
    ## [20] "sig.genes.rectum_Ensembl_Hs_UC"            
    ## [21] "sig.genes.rectum_Ensembl_PriID_BridgeDb_CD"
    ## [22] "sig.genes.rectum_Ensembl_PriID_BridgeDb_UC"
    ## [23] "sig.genes.rectum_ENTREZ_BridgeDb_CD"       
    ## [24] "sig.genes.rectum_ENTREZ_BridgeDb_UC"       
    ## [25] "sig.genes.rectum_ENTREZ_Hs_CD"             
    ## [26] "sig.genes.rectum_ENTREZ_Hs_UC"             
    ## [27] "sig.genes.rectum_ENTREZ_PriID_BridgeDb_CD" 
    ## [28] "sig.genes.rectum_ENTREZ_PriID_BridgeDb_UC" 
    ## [29] "sig.genes.rectum_HGNC_CD"                  
    ## [30] "sig.genes.rectum_HGNC_PriID_BridgeDb_CD"   
    ## [31] "sig.genes.rectum_HGNC_PriID_BridgeDb_UC"   
    ## [32] "sig.genes.rectum_HGNC_UC"

``` r
for (gene_list in ls(pattern = "sig.genes")){
  sig.genes = get (gene_list)
  IDsource = gsub ("_.*", "", gsub ("^[^_]+_", "", gene_list))
  item3 <- get (ls (pattern = paste0("item3_", IDsource)))
  if (IDsource == "HGNC") query <- paste0("hgnc:", sig.genes)
  if (IDsource == "Ensembl") query <- paste0("ensembl:", sig.genes)
  if (IDsource == "ENTREZ") query <- paste0("entrez:", sig.genes)

  split_query <- split(query, ceiling(seq_along(query) / 220))
  
  showresults_CombinePWs <- c()
  for (i in 1:length (split_query)) {
    string <- paste(split_query[[i]], collapse=' ')
    query_CombinePWs <- paste(item1, string, item2, item3, item4)
    results_CombinePWs <- SPARQL(endpointwp, query_CombinePWs, curl_args = list(useragent = R.version.string))
    showresults_CombinePWs <- rbind (showresults_CombinePWs, results_CombinePWs$results)
  }
  outputFile = paste0 ("results/CombinePWs", gsub ("sig.genes", "", gene_list))
  showresults_CombinePWs %>% 
    group_by(pathwayRes, pathway, pathwayTitle) %>% 
    summarise(GenesInPWs = sum(GenesInPWs)) %>%
    write.table(outputFile, 
            sep = "\t" , quote = FALSE, row.names = FALSE)
  assign(paste0 ("CombinePWs", gsub ("sig.genes", "", gene_list)), showresults_CombinePWs)
  rm (sig.genes, IDsource, item3, query, split_query, showresults_CombinePWs, i, string, query_CombinePWs, results_CombinePWs, outputFile)
}
```

##Pathway Mapping stats:

``` r
##TODO: the MappingStats table should be checked.

##HGNC_originalData
MappingStats <- data.table(`  ` =  c("The total number of genes in the transcriptomics dataset" ,
                                     "The total number of significant genes with HGNC Symbol IDs for CD (ileum location)",
                                     "The total number of pathways with HGNC Symbol IDs for CD (ileum location)",
                                     "The total number of significant genes with HGNC Symbol IDs for CD (rectum location)",
                                     "The total number of pathways with HGNC Symbol IDs for CD (rectum location)",
                                     "The total number of significant genes with HGNC Symbol IDs for UC (ileum location)",
                                     "The total number of pathways with HGNC Symbol IDs for UC (ileum location)",
                                     "The total number of significant genes with HGNC Symbol IDs for UC (rectum location)",
                                     "The total number of pathways with HGNC Symbol IDs for UC (rectum location)",
                                     
                                     "The total number of significant genes with primary HGNC Symbol IDs for CD (ileum location)",
                                     "The total number of pathways with primary HGNC Symbol IDs for CD (ileum location)",
                                     "The total number of significant genes with primary HGNC Symbol IDs for CD (rectum location)",
                                     "The total number of pathways with primary HGNC Symbol IDs for CD (rectum location)",
                                     "The total number of significant genes with primary HGNC Symbol IDs for UC (ileum location)",
                                     "The total number of pathways with primary HGNC Symbol IDs for UC (ileum location)",
                                     "The total number of significant genes with primary HGNC Symbol IDs for UC (rectum location)",
                                     "The total number of pathways with primary HGNC Symbol IDs for UC (rectum location)",
                                     
                                     "The total number of unique Entrez IDs in the transcriptomics dataset" ,
                                     "The total number of significant genes with Entrez IDs for CD (ileum location)",
                                     "The total number of pathways with Entrez IDs for CD (ileum location)",
                                     "The total number of significant genes with Entrez IDs for CD (rectum location)",
                                     "The total number of pathways with Entrez IDs for CD (rectum location)",
                                     "The total number of significant genes with Entrez IDs for UC (ileum location)",
                                     "The total number of pathways with Entrez IDs for UC (ileum location)",
                                     "The total number of significant genes with Entrez IDs for UC (rectum location)",
                                     "The total number of pathways with Entrez IDs for UC (rectum location)",
                                     "The total number of unique Ensembl IDs in the transcriptomics dataset" ,
                                     "The total number of significant genes with Ensembl IDs for CD (ileum location)",
                                     "The total number of pathways with Ensembl IDs for CD (ileum location)",
                                     "The total number of significant genes with Ensembl IDs for CD (rectum location)",
                                     "The total number of pathways with Ensembl IDs for CD (rectum location)",
                                     "The total number of significant genes with Ensembl IDs for UC (ileum location)",
                                     "The total number of pathways with Ensembl IDs for UC (ileum location)",
                                     "The total number of significant genes with Ensembl IDs for UC (rectum location)",
                                     "The total number of pathways with Ensembl IDs for UC (rectum location)"),
                           org.Hs = c(nrow (dataset_CD), 
                                      length(sig.genes.ileum_HGNC_CD), length(unique(CombinePWs.ileum_HGNC_CD$pathway)),
                                      length(sig.genes.rectum_HGNC_CD), length(unique(CombinePWs.rectum_HGNC_CD$pathway)),
                                      length(sig.genes.ileum_HGNC_UC), length(unique(CombinePWs.ileum_HGNC_UC$pathway)), 
                                      length(sig.genes.rectum_HGNC_UC), length(unique(CombinePWs.rectum_HGNC_UC$pathway)),
                                      "-", "-", "-", "-", "-", "-", "-", "-",
                                      nrow (dataset_CD[na.omit(dataset_CD$ENTREZ.ID_org.Hs),]), 
                                      length(sig.genes.ileum_ENTREZ_Hs_CD), length(unique(CombinePWs.ileum_ENTREZ_Hs_CD$pathway)),
                                      length(sig.genes.rectum_ENTREZ_Hs_CD), length(unique(CombinePWs.rectum_ENTREZ_Hs_CD$pathway)),
                                      length(sig.genes.ileum_ENTREZ_Hs_UC), length(unique(CombinePWs.ileum_ENTREZ_Hs_UC$pathway)),
                                      length(sig.genes.rectum_ENTREZ_Hs_UC), length(unique(CombinePWs.rectum_ENTREZ_Hs_UC$pathway)),
                                      nrow (dataset_CD[na.omit(dataset_CD$Ensembl.ID_org.Hs),]), 
                                      length(sig.genes.ileum_Ensembl_Hs_CD), length(unique(CombinePWs.ileum_Ensembl_Hs_CD$pathway)),
                                      length(sig.genes.rectum_Ensembl_Hs_CD), length(unique(CombinePWs.rectum_Ensembl_Hs_CD$pathway)),
                                      length(sig.genes.ileum_Ensembl_Hs_UC), length(unique(CombinePWs.ileum_Ensembl_Hs_UC$pathway)),
                                      length(sig.genes.rectum_Ensembl_Hs_UC), length(unique(CombinePWs.rectum_Ensembl_Hs_UC$pathway))),
                           BridgeDb = c(nrow (dataset_CD), 
                                      length(sig.genes.ileum_HGNC_CD), length(unique(CombinePWs.ileum_HGNC_CD$pathway)),
                                      length(sig.genes.rectum_HGNC_CD), length(unique(CombinePWs.rectum_HGNC_CD$pathway)),
                                      length(sig.genes.ileum_HGNC_UC), length(unique(CombinePWs.ileum_HGNC_UC$pathway)),
                                      length(sig.genes.rectum_HGNC_UC), length(unique(CombinePWs.rectum_HGNC_UC$pathway)),
                                      "-", "-", "-", "-", "-", "-", "-", "-",
                                      nrow (dataset_CD[na.omit(dataset_CD$ENTREZ.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_ENTREZ_BridgeDb_CD), length(unique(CombinePWs.ileum_ENTREZ_BridgeDb_CD$pathway)),
                                      length(sig.genes.rectum_ENTREZ_BridgeDb_CD), length(unique(CombinePWs.rectum_ENTREZ_BridgeDb_CD$pathway)),
                                      length(sig.genes.ileum_ENTREZ_BridgeDb_UC), length(unique(CombinePWs.ileum_ENTREZ_BridgeDb_UC$pathway)),
                                      length(sig.genes.rectum_ENTREZ_BridgeDb_UC), length(unique(CombinePWs.rectum_ENTREZ_BridgeDb_UC$pathway)),
                                      nrow (dataset_CD[na.omit(dataset_CD$Ensembl.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_Ensembl_BridgeDb_CD), length(unique(CombinePWs.ileum_Ensembl_BridgeDb_CD$pathway)),
                                      length(sig.genes.rectum_Ensembl_BridgeDb_CD), length(unique(CombinePWs.rectum_Ensembl_BridgeDb_CD$pathway)),
                                      length(sig.genes.ileum_Ensembl_BridgeDb_UC), length(unique(CombinePWs.ileum_Ensembl_BridgeDb_UC$pathway)),
                                      length(sig.genes.rectum_Ensembl_BridgeDb_UC), length(unique(CombinePWs.rectum_Ensembl_BridgeDb_UC$pathway))),
                           PrimaryID_BridgeDb = c(nrow (dataset_CD),
                                      length(sig.genes.ileum_HGNC_CD), length(unique(CombinePWs.ileum_HGNC_CD$pathway)),
                                      length(sig.genes.rectum_HGNC_CD), length(unique(CombinePWs.rectum_HGNC_CD$pathway)),
                                      length(sig.genes.ileum_HGNC_UC), length(unique(CombinePWs.ileum_HGNC_UC$pathway)),
                                      length(sig.genes.rectum_HGNC_UC), length(unique(CombinePWs.rectum_HGNC_UC$pathway)),
                                      nrow (dataset_CD),
                                      length(sig.genes.ileum_HGNC_PriID_BridgeDb_CD), length(unique(CombinePWs.ileum_HGNC_PriID_BridgeDb_CD$pathway)),
                                      length(sig.genes.rectum_HGNC_PriID_BridgeDb_CD), length(unique(CombinePWs.rectum_HGNC_PriID_BridgeDb_CD$pathway)),
                                      length(sig.genes.ileum_HGNC_PriID_BridgeDb_UC), length(unique(CombinePWs.ileum_HGNC_PriID_BridgeDb_UC$pathway)),
                                      length(sig.genes.rectum_HGNC_PriID_BridgeDb_UC), length(unique(CombinePWs.rectum_HGNC_PriID_BridgeDb_UC$pathway)),
                                      nrow (dataset_CD[na.omit(dataset_CD$ENTREZ.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_ENTREZ_PriID_BridgeDb_CD), length(unique(CombinePWs.ileum_ENTREZ_PriID_BridgeDb_CD$pathway)),
                                      length(sig.genes.rectum_ENTREZ_PriID_BridgeDb_CD), length(unique(CombinePWs.rectum_ENTREZ_PriID_BridgeDb_CD$pathway)),
                                      length(sig.genes.ileum_ENTREZ_PriID_BridgeDb_UC), length(unique(CombinePWs.ileum_ENTREZ_PriID_BridgeDb_UC$pathway)),
                                      length(sig.genes.rectum_ENTREZ_PriID_BridgeDb_UC), length(unique(CombinePWs.rectum_ENTREZ_PriID_BridgeDb_UC$pathway)),
                                      nrow (dataset_CD[na.omit(dataset_CD$Ensembl.ID_BridgeDb),]), 
                                      length(sig.genes.ileum_Ensembl_PriID_BridgeDb_CD), length(unique(CombinePWs.ileum_Ensembl_PriID_BridgeDb_CD$pathway)),
                                      length(sig.genes.rectum_Ensembl_PriID_BridgeDb_CD), length(unique(CombinePWs.rectum_Ensembl_PriID_BridgeDb_CD$pathway)),
                                      length(sig.genes.ileum_Ensembl_PriID_BridgeDb_UC), length(unique(CombinePWs.ileum_Ensembl_PriID_BridgeDb_UC$pathway)),
                                      length(sig.genes.rectum_Ensembl_PriID_BridgeDb_UC), length(unique(CombinePWs.rectum_Ensembl_PriID_BridgeDb_UC$pathway))))
kable(MappingStats)
```

|                                                                                             | org.Hs | BridgeDb | PrimaryID_BridgeDb |
|:-------------------------------------------------|:----|:-----|-----------:|
| The total number of genes in the transcriptomics dataset                                    | 17670  | 17670    |              17670 |
| The total number of significant genes with HGNC Symbol IDs for CD (ileum location)          | 1896   | 1896     |               1896 |
| The total number of pathways with HGNC Symbol IDs for CD (ileum location)                   | 611    | 611      |                611 |
| The total number of significant genes with HGNC Symbol IDs for CD (rectum location)         | 1339   | 1339     |               1339 |
| The total number of pathways with HGNC Symbol IDs for CD (rectum location)                  | 593    | 593      |                593 |
| The total number of significant genes with HGNC Symbol IDs for UC (ileum location)          | 121    | 121      |                121 |
| The total number of pathways with HGNC Symbol IDs for UC (ileum location)                   | 59     | 59       |                 59 |
| The total number of significant genes with HGNC Symbol IDs for UC (rectum location)         | 6809   | 6809     |               6809 |
| The total number of pathways with HGNC Symbol IDs for UC (rectum location)                  | 728    | 728      |                728 |
| The total number of significant genes with primary HGNC Symbol IDs for CD (ileum location)  | \-     | \-       |              17670 |
| The total number of pathways with primary HGNC Symbol IDs for CD (ileum location)           | \-     | \-       |               1896 |
| The total number of significant genes with primary HGNC Symbol IDs for CD (rectum location) | \-     | \-       |                611 |
| The total number of pathways with primary HGNC Symbol IDs for CD (rectum location)          | \-     | \-       |               1339 |
| The total number of significant genes with primary HGNC Symbol IDs for UC (ileum location)  | \-     | \-       |                594 |
| The total number of pathways with primary HGNC Symbol IDs for UC (ileum location)           | \-     | \-       |                121 |
| The total number of significant genes with primary HGNC Symbol IDs for UC (rectum location) | \-     | \-       |                  3 |
| The total number of pathways with primary HGNC Symbol IDs for UC (rectum location)          | \-     | \-       |               6809 |
| The total number of unique Entrez IDs in the transcriptomics dataset                        | 15022  | 14133    |                145 |
| The total number of significant genes with Entrez IDs for CD (ileum location)               | 1697   | 1558     |              14133 |
| The total number of pathways with Entrez IDs for CD (ileum location)                        | 611    | 611      |               1558 |
| The total number of significant genes with Entrez IDs for CD (rectum location)              | 1195   | 1114     |                611 |
| The total number of pathways with Entrez IDs for CD (rectum location)                       | 593    | 593      |               1114 |
| The total number of significant genes with Entrez IDs for UC (ileum location)               | 99     | 97       |                593 |
| The total number of pathways with Entrez IDs for UC (ileum location)                        | 59     | 59       |                 97 |
| The total number of significant genes with Entrez IDs for UC (rectum location)              | 5870   | 5593     |                 59 |
| The total number of pathways with Entrez IDs for UC (rectum location)                       | 728    | 728      |               5593 |
| The total number of unique Ensembl IDs in the transcriptomics dataset                       | 14596  | 15009    |                728 |
| The total number of significant genes with Ensembl IDs for CD (ileum location)              | 1665   | 1695     |              15009 |
| The total number of pathways with Ensembl IDs for CD (ileum location)                       | 611    | 611      |               1695 |
| The total number of significant genes with Ensembl IDs for CD (rectum location)             | 1171   | 1195     |                611 |
| The total number of pathways with Ensembl IDs for CD (rectum location)                      | 593    | 593      |               1195 |
| The total number of significant genes with Ensembl IDs for UC (ileum location)              | 97     | 99       |                593 |
| The total number of pathways with Ensembl IDs for UC (ileum location)                       | 59     | 59       |                 99 |
| The total number of significant genes with Ensembl IDs for UC (rectum location)             | 5771   | 5864     |                 59 |
| The total number of pathways with Ensembl IDs for UC (rectum location)                      | 728    | 728      |               5864 |
| The total number of genes in the transcriptomics dataset                                    | 17670  | 17670    |                728 |