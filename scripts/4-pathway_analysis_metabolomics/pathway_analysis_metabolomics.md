## Introduction

In this workflow, we link the metabolites of interest to pathway data
from WikiPathways, based on their HMDB and ChEBI identifiers.

## R environment setup

``` r
# empty the R environment
rm (list = ls())

# check if libraries are already installed, otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"SPARQL" %in% installed.packages()) install.packages("SPARQL")
if(!"dplyr" %in% installed.packages()) install.packages("dplyr")
if(!"data.table" %in% installed.packages()) install.packages("data.table")
if(!"knitr" %in% installed.packages()) install.packages("knitr")
if(!"reshape2" %in% installed.packages()) install.packages("reshape2")
if(!"ggplot2" %in% installed.packages()) install.packages("ggplot2")

#loading installed libraries
suppressPackageStartupMessages({
  library(rstudioapi)
  library(SPARQL)
  library(dplyr)
  library(data.table)
  library(knitr)
  library(reshape2)
  library(ggplot2)
})

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset and creating the identifier lists for the pathway analysis

``` r
#We have two datasets (CD and UC disorders)
mbx_dataset_CD <- read.delim("../3-identifier_mapping_metabolomics/results/mbx_IDMapping_CD")
mbx_dataset_UC <- read.delim("../3-identifier_mapping_metabolomics/results/mbx_IDMapping_UC")

#make list of metabolites for the pathway analysis
### HMDB IDs:
sig.metabolites.HMDB_CD <- na.omit(unique(mbx_dataset_CD$HMDBID)) #CD
sig.metabolites.HMDB_UC <- na.omit(unique(mbx_dataset_UC$HMDBID)) #UC
## BridgeDb
### ChEBI IDs:
sig.metabolites.ChEBI_BridgeDb_CD <- na.omit(unique(mbx_dataset_CD$ChEBI_BridgeDb)) #CD
sig.metabolites.ChEBI_BridgeDb_UC <- na.omit(unique(mbx_dataset_UC$ChEBI_BridgeDb)) #UC
## Primary id mapping, BridgeDb,
### primary HMDB IDs:
sig.metabolites.HMDB_PriID_BridgeDb_CD <- na.omit(unique(mbx_dataset_CD$Current_HMDBID)) #CD
sig.metabolites.HMDB_PriID_BridgeDb_UC <- na.omit(unique(mbx_dataset_UC$Current_HMDBID)) #UC
### ChEBI IDs:
sig.metabolites.ChEBI_PriID_BridgeDb_CD <- na.omit(unique(mbx_dataset_CD$ChEBI_PriID_BridgeDb)) #CD
sig.metabolites.ChEBI_PriID_BridgeDb_UC <- na.omit(unique(mbx_dataset_UC$ChEBI_PriID_BridgeDb)) #UC
```

## Find pathways for each dataset, based on different IDs.

``` r
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
options(RCurlOptions = list(cainfo = paste0(tempdir(), "/cacert.pem" ), ssl.verifypeer = FALSE))
resultsMetadata <- SPARQL(endpointwp, queryMetadata, curl_args = list(useragent = R.version.string))
showresultsMetadata <- resultsMetadata$results
remove(queryMetadata, resultsMetadata)

#For now, filter out Reactome PWs due to visualization issues in Cytoscape.
item1 = "PREFIX ch: <https://identifiers.org/hmdb/>
PREFIX ce: <https://identifiers.org/chebi/CHEBI:>
PREFIX cur: <http://vocabularies.wikipathways.org/wp#Curation:>
select distinct ?pathwayRes (str(?wpid) as ?pathway) (str(?title) as ?pathwayTitle) (count(distinct ?metaboliteID) AS ?BiomarkersInPWs) (count(distinct ?metaboliteDatanode) AS ?TotalMetabolitesinPW) where {
VALUES ?metaboliteID {"
item2 = "}
 
 ?metaboliteDatanode    a wp:Metabolite ;
                        dcterms:isPartOf ?pathwayRes .
 
 ?datanode  dcterms:isPartOf ?pathwayRes ;   "
item3_HMDB= "
    wp:bdbHmdb  ?metaboliteID ."
item3_ChEBI= "
     wp:bdbChEBI ?metaboliteID ." 

item4=
 " ?pathwayRes a wp:Pathway ;
             wp:organismName 'Homo sapiens' ; 
            dcterms:identifier ?wpid ;
            dc:title ?title .

  #?pathwayRes wp:ontologyTag cur:Reactome_Approved . 
  ?pathwayRes wp:ontologyTag cur:AnalysisCollection .           
}
ORDER BY DESC(?BiomarkersInPWs)"


##Split significant metabolites into list of max. 220 entries, to avoid SPARQL endpoint trowing a 414 error. 
##Merge the content of the split content back together for the output of the PW Analysis.

for (metabolite_list in ls(pattern = "sig.metabolites")){
  sig.metabolites = get (metabolite_list)
  IDsource = gsub ("_.*", "", gsub (".*\\.", "", metabolite_list))
  item3 <- get (ls (pattern = paste0("item3_", IDsource)))
  if (IDsource == "HMDB") query <- paste0("ch:", sig.metabolites)
  # if (IDsource == "ChEBI") query <- paste0("CHEBI:", sig.metabolites)
  if (IDsource == "ChEBI") query <- gsub ("CHEBI", "ce", sig.metabolites)
  # if (IDsource == "ChEBI") query <- sig.metabolites

  split_query <- split(query, ceiling(seq_along(query) / 220))
  
  showresults_CombinePWs <- c()
  for (i in 1:length (split_query)) {
    string <- paste(split_query[[i]], collapse=' ')
    query_CombinePWs <- paste(item1, string, item2, item3, item4)
    results_CombinePWs <- SPARQL(endpointwp, query_CombinePWs, curl_args = list(useragent = R.version.string))
    showresults_CombinePWs <- rbind (showresults_CombinePWs, results_CombinePWs$results)
  }
  outputFile = paste0 ("results/CombinePWs", gsub ("sig.metabolites", "", metabolite_list), ".txt")
  (showresults_CombinePWs <- showresults_CombinePWs %>% 
      group_by(pathwayRes, pathway, pathwayTitle, TotalMetabolitesinPW) %>% 
      summarise(BiomarkersInPWs = sum(BiomarkersInPWs)) %>%
      mutate(probabilities = dhyper(BiomarkersInPWs, TotalMetabolitesinPW, (length(query) - BiomarkersInPWs), length(query), log = FALSE)) %>% # Calculate hypergeometric density p-value for all pathways.
      arrange(desc(BiomarkersInPWs), probabilities)) %>%
    write.table(outputFile, sep = "\t" , quote = FALSE, row.names = FALSE)
  
  colnames(showresults_CombinePWs)[colnames(showresults_CombinePWs) == "BiomarkersInPWs"] <- paste0 ("BiomarkersInPWs", gsub ("sig.metabolites", "", metabolite_list))
  colnames(showresults_CombinePWs)[colnames(showresults_CombinePWs) == "probabilities"] <- paste0 ("probabilities", gsub ("sig.metabolites", "", metabolite_list))
  
  assign(paste0 ("CombinePWs", gsub ("sig.metabolites", "", metabolite_list)), showresults_CombinePWs)
  rm (sig.metabolites, IDsource, item3, query, split_query, showresults_CombinePWs, i, string, query_CombinePWs, results_CombinePWs, outputFile)
}
```

##Pathway Mapping stats:

``` r
MappingStats <- data.table(`  ` =  c("#significant metabolites with HMDB IDs for CD (primary in PrimaryID_BridgeDb)" ,
                                     "#pathways with HMDB IDs for CD (primary in PrimaryID_BridgeDb)",
                                     "#significant metabolites with HMDB IDs for UC (primary in PrimaryID_BridgeDb)",
                                     "#pathways with HMDB IDs for UC (primary in PrimaryID_BridgeDb)",

                                     "#significant metabolites with ChEBI IDs for CD",
                                     "#pathways with ChEBI IDs for CD",
                                     "#significant metabolites with ChEBI IDs for UC",
                                     "#pathways with ChEBI IDs for UC"),
                           BridgeDb = c(length(sig.metabolites.HMDB_CD), 
                                        length(unique(CombinePWs.HMDB_CD$pathway)),
                                        length(sig.metabolites.HMDB_UC), 
                                        length(unique(CombinePWs.HMDB_UC$pathway)),
                                        
                                        length(sig.metabolites.ChEBI_BridgeDb_CD),
                                        length(unique(CombinePWs.ChEBI_BridgeDb_CD$pathway)),
                                        length(sig.metabolites.ChEBI_BridgeDb_UC),
                                        length(unique(CombinePWs.ChEBI_BridgeDb_UC$pathway))),
                           PrimaryID_BridgeDb = c(length(sig.metabolites.HMDB_PriID_BridgeDb_CD),
                                                  length(unique(CombinePWs.HMDB_PriID_BridgeDb_CD$pathway)),
                                                  length(sig.metabolites.HMDB_PriID_BridgeDb_UC),
                                                  length(unique(CombinePWs.HMDB_PriID_BridgeDb_UC$pathway)),
                                                  
                                                  length(sig.metabolites.ChEBI_PriID_BridgeDb_CD),
                                                  length(unique(CombinePWs.ChEBI_PriID_BridgeDb_CD$pathway)),
                                                  length(sig.metabolites.ChEBI_PriID_BridgeDb_UC),
                                                  length(unique(CombinePWs.ChEBI_PriID_BridgeDb_UC$pathway))))
kable(MappingStats)
```

|                                                                               | BridgeDb | PrimaryID_BridgeDb |
|:---------------------------------------------------|------:|-------------:|
| #significant metabolites with HMDB IDs for CD (primary in PrimaryID_BridgeDb) |      438 |                436 |
| #pathways with HMDB IDs for CD (primary in PrimaryID_BridgeDb)                |      229 |                229 |
| #significant metabolites with HMDB IDs for UC (primary in PrimaryID_BridgeDb) |      437 |                435 |
| #pathways with HMDB IDs for UC (primary in PrimaryID_BridgeDb)                |      229 |                229 |
| #significant metabolites with ChEBI IDs for CD                                |      363 |                366 |
| #pathways with ChEBI IDs for CD                                               |      228 |                228 |
| #significant metabolites with ChEBI IDs for UC                                |      362 |                365 |
| #pathways with ChEBI IDs for UC                                               |      228 |                228 |

## Combining pathways for the vidualization

``` r
CombinePWs_CD <- Reduce(function(x, y) merge(x, y, all = T), lapply(ls(pattern = "Combin.*CD"), get), accumulate = F)
CombinePWs_UC <- Reduce(function(x, y) merge(x, y, all = T), lapply(ls(pattern = "Combin.*UC"), get), accumulate = F)

#CD
CombinePWs_CD_toPlot <- CombinePWs_CD [, !grepl("probabilities|pathwayRes|pathwayTitle|TotalMetabolitesinPW", colnames (CombinePWs_CD))] %>% select (pathway, BiomarkersInPWs.HMDB_CD, BiomarkersInPWs.ChEBI_BridgeDb_CD, BiomarkersInPWs.HMDB_PriID_BridgeDb_CD, BiomarkersInPWs.ChEBI_PriID_BridgeDb_CD)
colnames(CombinePWs_CD_toPlot) <- gsub ("BiomarkersInPWs.|_CD", "", colnames(CombinePWs_CD_toPlot))
all_the_same <- apply(CombinePWs_CD_toPlot[-c(1, 2)], 1, function(x) all(x == x[1]))
#Keep the pathways with different number of BiomarkersInPWs across different mappings
CombinePWs_CD_toPlot <- CombinePWs_CD_toPlot %>%
  filter(!all_the_same) %>%
  reshape2::melt() %>%
  mutate(pathway = paste0(pathway, ":", CombinePWs_CD$pathwayTitle[match(pathway, CombinePWs_CD$pathway)]),
         value = ifelse (is.na(value), 0, value))
ggplot(CombinePWs_CD_toPlot, aes(x = variable, y = pathway, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(y = NULL, x = NULL) +
  geom_text(aes(label = value), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none") +
  ggtitle("CD")
```

![](pathway_analysis_metabolomics_files/figure-markdown_github/Comparison-1.png)

``` r
#UC
CombinePWs_UC_toPlot <- CombinePWs_UC [, !grepl("probabilities|pathwayRes|pathwayTitle|TotalMetabolitesinPW", colnames (CombinePWs_UC))] %>% select (pathway, BiomarkersInPWs.HMDB_UC, BiomarkersInPWs.ChEBI_BridgeDb_UC, BiomarkersInPWs.HMDB_PriID_BridgeDb_UC, BiomarkersInPWs.ChEBI_PriID_BridgeDb_UC)
colnames(CombinePWs_UC_toPlot) <- gsub ("BiomarkersInPWs.|_UC", "", colnames(CombinePWs_UC_toPlot))
#Keep the pathways with different number of BiomarkersInPWs across different mappings
all_the_same <- apply(CombinePWs_UC_toPlot[-c(1, 2)], 1, function(x) all(x == x[1]))
CombinePWs_UC_toPlot <- CombinePWs_UC_toPlot %>%
  filter(!all_the_same) %>%
  reshape2::melt() %>%
  mutate(pathway = paste0(pathway, ":", CombinePWs_UC$pathwayTitle[match(pathway, CombinePWs_UC$pathway)]),
         value = ifelse (is.na(value), 0, value))
ggplot(CombinePWs_UC_toPlot, aes(x = variable, y = pathway, fill = value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(y = NULL, x = NULL) +
  geom_text(aes(label = value), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none") +
  ggtitle("UC")
```

![](pathway_analysis_metabolomics_files/figure-markdown_github/Comparison-2.png)