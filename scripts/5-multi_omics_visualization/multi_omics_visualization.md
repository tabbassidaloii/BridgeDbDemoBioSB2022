## Introduction

In this script, we visualize multi-omics data in Cytoscape through the
WikiPathways App for Cytoscape. Through this app, pathways from
WikiPathways, Reactome and LIPID MAPS are available. Identifiers from
these three pathway model databases are unified and harmonized to
Ensembl and ChEBI IDs, to allow for multi-omics data visualization. The
transcriptomics data and metabolites data are combined to automatically
visualize their respective log2FC and p-values.

## R environment setup

``` r
#Empty the R environment
rm (list = ls())
#Check if libraries are already installed, otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"RCy3" %in% installed.packages()) BiocManager::install("RCy3")
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"RColorBrewer" %in% installed.packages()) install.packages("RColorBrewer")
if(!"dplyr" %in% installed.packages()) install.packages("dplyr")

#Load installed libraries
suppressPackageStartupMessages({
  library(rstudioapi)
  library(RCy3) #Connect cytoscape via R
  library(rWikiPathways) #Get pathways from WikiPathways
  library(RColorBrewer)
  library(dplyr)
})

#Set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing data from step 1 (Transcriptomics identifier mapping) and 3 (Metabolomics identifier mapping)

``` r
##TODO: start workflow on data in file from this folder (in case step 1 and 3 are not run yet?)

#Transcriptomics and take a subset of data, which is relevant for the multiomics-visualization:
tSet_ileum_CD <- read.delim("../1-identifier_mapping_transcriptomics/results/IDMapping_CD") %>%
  select(Ensembl.ID_PriID_BridgeDb, log2FC_ileum, pvalue_ileum) %>%
  rename(Identifier = Ensembl.ID_PriID_BridgeDb, log2FC = log2FC_ileum, pvalue = pvalue_ileum) %>%
  filter(!is.na(Identifier))
tSet_rectum_CD <- read.delim("../1-identifier_mapping_transcriptomics/results/IDMapping_CD") %>%
  select(Ensembl.ID_PriID_BridgeDb, log2FC_rectum, pvalue_rectum) %>%
  rename(Identifier = Ensembl.ID_PriID_BridgeDb, log2FC = log2FC_rectum, pvalue = pvalue_rectum) %>%
  filter(!is.na(Identifier))
tSet_ileum_UC <- read.delim("../1-identifier_mapping_transcriptomics/results/IDMapping_UC") %>%
  select(Ensembl.ID_PriID_BridgeDb, log2FC_ileum, pvalue_ileum) %>%
  rename(Identifier = Ensembl.ID_PriID_BridgeDb, log2FC = log2FC_ileum, pvalue = pvalue_ileum) %>%
  filter(!is.na(Identifier))
tSet_rectum_UC <- read.delim("../1-identifier_mapping_transcriptomics/results/IDMapping_UC") %>%
  select(Ensembl.ID_PriID_BridgeDb, log2FC_rectum, pvalue_rectum) %>%
  rename(Identifier = Ensembl.ID_PriID_BridgeDb, log2FC = log2FC_rectum, pvalue = pvalue_rectum) %>%
  filter(!is.na(Identifier))

#Metabolomics and take a subset of data, which is relevant for the multiomics-visualization:
mSet_CD <- read.delim("../3-identifier_mapping_metabolomics/results/mbx_IDMapping_CD") %>%
  select(ChEBI_PriID_BridgeDb, log2FC, pvalue) %>%
  rename(Identifier = ChEBI_PriID_BridgeDb) %>%
  filter(!is.na(Identifier))
mSet_UC <- read.delim("../3-identifier_mapping_metabolomics/results/mbx_IDMapping_UC") %>%
  select(ChEBI_PriID_BridgeDb, log2FC, pvalue) %>%
  rename(Identifier = ChEBI_PriID_BridgeDb) %>%
  filter(!is.na(Identifier))

#Merge two OMICS types together
Combine_ileum_CD_multi <- rbind(tSet_ileum_CD, mSet_CD)
Combine_rectum_CD_multi <- rbind(tSet_rectum_CD, mSet_CD)
Combine_ileum_UC_multi <- rbind(tSet_ileum_UC, mSet_UC)
Combine_rectum_UC_multi <- rbind(tSet_rectum_UC, mSet_UC)
```

## Selecting relevant pathway based on step 2 (Transcriptomics pathway analysis) and 4 (Metabolomics pathway analysis)

Please run steps 2 and 4 first, in order to retrieve updated data for
the following step.

``` r
#Significant pathways (changes due to disease)
#Transcriptomics
CombinePWs_ileum_CD_transcriptomics <- read.delim("../2-pathway_analysis_transcriptomics/results/CombinePWs.ileum_CD.txt") %>%
  select(pathwayRes, pathway, pathwayTitle, TotalGenesinPW, GenesInPWs.ileum_Ensembl_PriID_BridgeDb_CD, probabilities.ileum_Ensembl_PriID_BridgeDb_CD)
colnames(CombinePWs_ileum_CD_transcriptomics) <- gsub (".ileum_Ensembl_PriID_BridgeDb_CD", "", colnames(CombinePWs_ileum_CD_transcriptomics))
colnames(CombinePWs_ileum_CD_transcriptomics)[colnames(CombinePWs_ileum_CD_transcriptomics) == "probabilities"] = "probabilities_transcriptome"

CombinePWs_rectum_CD_transcriptomics <- read.delim("../2-pathway_analysis_transcriptomics/results/CombinePWs.rectum_CD.txt") %>%
  select(pathwayRes, pathway, pathwayTitle, TotalGenesinPW, GenesInPWs.rectum_Ensembl_PriID_BridgeDb_CD, probabilities.rectum_Ensembl_PriID_BridgeDb_CD)
colnames(CombinePWs_rectum_CD_transcriptomics) <- gsub (".rectum_Ensembl_PriID_BridgeDb_CD", "", colnames(CombinePWs_rectum_CD_transcriptomics))
colnames(CombinePWs_rectum_CD_transcriptomics)[colnames(CombinePWs_rectum_CD_transcriptomics) == "probabilities"] = "probabilities_transcriptome"

CombinePWs_ileum_UC_transcriptomics <- read.delim("../2-pathway_analysis_transcriptomics/results/CombinePWs.ileum_UC.txt") %>%
  select(pathwayRes, pathway, pathwayTitle, TotalGenesinPW, GenesInPWs.ileum_Ensembl_PriID_BridgeDb_UC, probabilities.ileum_Ensembl_PriID_BridgeDb_UC)
colnames(CombinePWs_ileum_UC_transcriptomics) <- gsub (".ileum_Ensembl_PriID_BridgeDb_UC", "", colnames(CombinePWs_ileum_UC_transcriptomics))
colnames(CombinePWs_ileum_UC_transcriptomics)[colnames(CombinePWs_ileum_UC_transcriptomics) == "probabilities"] = "probabilities_transcriptome"

CombinePWs_rectum_UC_transcriptomics <- read.delim("../2-pathway_analysis_transcriptomics/results/CombinePWs.rectum_UC.txt") %>%
  select(pathwayRes, pathway, pathwayTitle, TotalGenesinPW, GenesInPWs.rectum_Ensembl_PriID_BridgeDb_UC, probabilities.rectum_Ensembl_PriID_BridgeDb_UC)
colnames(CombinePWs_rectum_UC_transcriptomics) <- gsub (".rectum_Ensembl_PriID_BridgeDb_UC", "", colnames(CombinePWs_rectum_UC_transcriptomics))
colnames(CombinePWs_rectum_UC_transcriptomics)[colnames(CombinePWs_rectum_UC_transcriptomics) == "probabilities"] = "probabilities_transcriptome"

#Metabolomics
CombinePWs_CD_metabolomics <- read.delim("../4-pathway_analysis_metabolomics/results/CombinePWs_CD_mbx.txt") %>% 
  select(pathwayRes, pathway, pathwayTitle, TotalMetabolitesinPW, BiomarkersInPWs.ChEBI_PriID_BridgeDb_CD, probabilities.ChEBI_PriID_BridgeDb_CD)
colnames(CombinePWs_CD_metabolomics) <- gsub (".ChEBI_PriID_BridgeDb_CD", "", colnames(CombinePWs_CD_metabolomics))
colnames(CombinePWs_CD_metabolomics)[colnames(CombinePWs_CD_metabolomics) == "probabilities"] = "probabilities_metabolomics"

CombinePWs_UC_metabolomics <- read.delim("../4-pathway_analysis_metabolomics/results/CombinePWs_UC_mbx.txt") %>% 
  select(pathwayRes, pathway, pathwayTitle, TotalMetabolitesinPW, BiomarkersInPWs.ChEBI_PriID_BridgeDb_UC, probabilities.ChEBI_PriID_BridgeDb_UC)
colnames(CombinePWs_UC_metabolomics) <- gsub (".ChEBI_PriID_BridgeDb_UC", "", colnames(CombinePWs_UC_metabolomics))
colnames(CombinePWs_UC_metabolomics)[colnames(CombinePWs_UC_metabolomics) == "probabilities"] = "probabilities_metabolomics"

#Interesting pathways 
#Transcriptomics
Int_PWs_ileum_CD_transcriptomics <- read.csv("../2-pathway_analysis_transcriptomics/results/Int_PWs_ileum_CD.txt")$x
Int_PWs_rectum_CD_transcriptomics <- read.csv("../2-pathway_analysis_transcriptomics/results/Int_PWs_rectum_CD.txt")$x
Int_PWs_ileum_UC_transcriptomics <- read.csv("../2-pathway_analysis_transcriptomics/results/Int_PWs_ileum_UC.txt")$x
Int_PWs_rectum_UC_transcriptomics <- read.csv("../2-pathway_analysis_transcriptomics/results/Int_PWs_rectum_UC.txt")$x

#Metabolomics
Int_PWs_CD_metabolomics <- read.csv("../4-pathway_analysis_metabolomics/results/Int_PWs_CD_mbx.txt")$x
Int_PWs_UC_metabolomics <- read.csv("../4-pathway_analysis_metabolomics/results/Int_PWs_UC_mbx.txt")$x

#Merge interesting pathways for both data types
Int_PWs_ileum_CD <- unique(c(Int_PWs_ileum_CD_transcriptomics, Int_PWs_CD_metabolomics))
Int_PWs_rectum_CD <- unique(c(Int_PWs_rectum_CD_transcriptomics, Int_PWs_CD_metabolomics))
Int_PWs_ileum_UC <- unique(c(Int_PWs_ileum_UC_transcriptomics, Int_PWs_UC_metabolomics))
Int_PWs_rectum_UC <- unique(c(Int_PWs_rectum_UC_transcriptomics, Int_PWs_UC_metabolomics))

##Merge two OMICS types together:
(CombinePWs_ileum_CD_multi <- merge (CombinePWs_CD_metabolomics, CombinePWs_ileum_CD_transcriptomics) %>% 
  arrange (desc(GenesInPWs), desc(BiomarkersInPWs))) %>% 
  write.table("results/CombinePWs_ileum_CD_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)
(Int_CombinePWs_ileum_CD_multi <- CombinePWs_ileum_CD_multi %>%
    filter(pathway %in% Int_PWs_ileum_CD)) %>% 
  write.table("results/Int_CombinePWs_ileum_CD_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)

(CombinePWs_rectum_CD_multi <- merge (CombinePWs_CD_metabolomics, CombinePWs_rectum_CD_transcriptomics) %>% 
  arrange (desc(GenesInPWs), desc(BiomarkersInPWs))) %>% 
  write.table("results/CombinePWs_rectum_CD_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)
(Int_CombinePWs_rectum_CD_multi <- CombinePWs_rectum_CD_multi %>%
    filter(pathway %in% Int_PWs_rectum_CD)) %>% 
  write.table("results/Int_CombinePWs_rectum_CD_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)

(CombinePWs_ileum_UC_multi <- merge (CombinePWs_UC_metabolomics, CombinePWs_ileum_UC_transcriptomics) %>% 
  arrange (desc(GenesInPWs), desc(BiomarkersInPWs))) %>% 
  write.table("results/CombinePWs_ileum_UC_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)
(Int_CombinePWs_ileum_UC_multi <- CombinePWs_ileum_UC_multi %>%
    filter(pathway %in% Int_PWs_ileum_UC)) %>% 
  write.table("results/Int_CombinePWs_ileum_UC_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)

(CombinePWs_rectum_UC_multi <- merge (CombinePWs_UC_metabolomics, CombinePWs_rectum_UC_transcriptomics) %>% 
  arrange (desc(GenesInPWs), desc(BiomarkersInPWs))) %>% 
  write.table("results/CombinePWs_rectum_UC_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)
(Int_CombinePWs_rectum_UC_multi <- CombinePWs_rectum_UC_multi %>%
    filter(pathway %in% Int_PWs_rectum_UC)) %>% 
  write.table("results/Int_CombinePWs_rectum_UC_multi.txt", sep = "\t" , quote = FALSE, row.names = FALSE)
```

## Importing pathway in Cytoscape

``` r
#Make sure to first launch Cytoscape outside of Rstudio, v.3.9.1 and the CyREST is installed
cytoscapePing()
#Check metadata of tool
cytoscapeVersionInfo()
```

    ##       apiVersion cytoscapeVersion 
    ##             "v1"          "3.9.1"

``` r
#Close all opened session before starting
closeSession(FALSE)
#Set up WikiPathways app in Cytoscape, v.3.3.10
if("WikiPathways" %in% commandsHelp("")) print("Success: the WikiPathways app is installed") else print("Warning: WikiPathways app is not installed. Please install the WikiPathways app before proceeding.")
```

    ## [1] "Available namespaces:"
    ## [1] "Warning: WikiPathways app is not installed. Please install the WikiPathways app before proceeding."

``` r
if(!"WikiPathways" %in% commandsHelp("")) installApp("WikiPathways")
```

    ## [1] "Available namespaces:"

``` r
#Pathway IDs to be visualized
pathway.id <- "WP4726"# Sphingolipid metabolism: integrated pathway --> Selected based on overlap between important PWs for both transcriptomics and metabolomic 
#Select a dataset to visualize on the pathway:
names_of_dataframes <- "Combine_rectum_UC_multi"
#Import pathways as pathway in cytoscape
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',pathway.id)) 
```

## Data upload with regular xRef ID mapping

``` r
#Select dataset to visualize
dataset <- get (names_of_dataframes)
#Remove duplicate rows
dataset <- dataset %>% distinct(Identifier, .keep_all = TRUE)
#Load data to the imported pathway in cytoscape by key column as omics.ID
loadTableData(table = "node", data = dataset, data.key.column = "Identifier", table.key.column = "XrefId")
```

    ## [1] "Success: Data loaded in defaultnode table"

## Visualization options

``` r
#New visual style is created
RCy3::copyVisualStyle("default", "pathwayStyle")
#Set new style as the current style
RCy3::setVisualStyle("pathwayStyle")
```

    ##                 message 
    ## "Visual Style applied."

``` r
#Set node dimensions as fixed sizes
RCy3::lockNodeDimensions(TRUE, style.name="pathwayStyle")

#Node shape mapping
RCy3::setNodeShapeMapping('Type', c('GeneProduct','Protein', 'Metabolite'), c('ELLIPSE','ELLIPSE','RECTANGLE'), style.name = "pathwayStyle")
```

    ## NULL

``` r
#Change node height
RCy3::setNodeHeightMapping('Type', c('GeneProduct','Protein', 'Metabolite'), c(23, 23, 25), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Change node width
RCy3::setNodeWidthMapping('Type', c('GeneProduct','Protein', 'Metabolite'), c(60, 60, 100), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Set node color based on log2FC for both genes and metabolites
node.colors <- c(rev(brewer.pal(3, "RdBu")))
setNodeColorMapping("log2FC", c(-1, 0, 1), node.colors, default.color = "#D3D3D3", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Set node border width and color based on p-value
#First we need to get all p-values from node table
pvalues <- getTableColumns(table = 'node', columns = 'pvalue')
pvalues <- na.omit(pvalues)
#Create a range for all sign. p-values, and one for all not significant.
significant_pvalues <- pvalues[(pvalues < 0.05)]
not.significant_pvalues <- pvalues[(pvalues >= 0.05)]
significant_pvalues.colors <- rep("#2e9d1d", length(significant_pvalues))
not.significant_pvalues.colors <- rep("#FFFFFF", length(not.significant_pvalues))

setNodeBorderWidthMapping('pvalue', table.column.values = NULL , c(6,6) , mapping.type = "c", style.name = "pathwayStyle")
```

    ## NULL

``` r
setNodeBorderColorMapping('pvalue', c(significant_pvalues,not.significant_pvalues), c(significant_pvalues.colors, not.significant_pvalues.colors), default.color = "#AAAAAA", mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Save output 
filename_multiomics <- paste0("results/", pathway.id, "_", names_of_dataframes, "_regular_visualization.png")
png.file <- file.path(getwd(), filename_multiomics)
exportImage(png.file, 'PNG', zoom = 500)
```

    ##                                                                                                                                           file 
    ## "D:\\bridgeDb\\BridgeDbDemoBioSB2022\\scripts\\5-multi_omics_visualization\\results\\WP4726_Combine_rectum_UC_multi_regular_visualization.png"

## Data upload with multiomics unified and harmonized mappings:

``` r
#Open a new pathway in Cytoscape for enhanced visualization.
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=', pathway.id)) 

#Get node table from imported pathway in cytoscape
ID.cols <- getTableColumns(table ="node", columns = c("XrefId", "Ensembl", "ChEBI"))
#Filter out rows which contain NA value for columns Ensembl and ChEBI
ID.cols <- ID.cols[!with(ID.cols, is.na(Ensembl) & is.na(ChEBI)),]
#If a row value in the Ensembl column is NA then replace it with ChEBI  
ID.cols$Ensembl <- ifelse(is.na(ID.cols$Ensembl), ID.cols$ChEBI, ID.cols$Ensembl)
#Use the only one column contains both Ensembl and ChEBI identifiers
ID.cols <- data.frame(ID.cols[, c(1, 2)])
#Change column name
colnames(ID.cols)[2] <- "omics.ID"
#Load all the multi-omics IDs in cytoscape by key column as XrefId
loadTableData(table = "node", data = ID.cols, data.key.column = "XrefId", table.key.column = "XrefId")
```

    ## [1] "Success: Data loaded in defaultnode table"

``` r
#Select dataset to visualize
dataset <- get (names_of_dataframes)
#Remove duplicate rows
dataset <- dataset %>% distinct(Identifier, .keep_all = TRUE)
#Load data to the imported pathway in cytoscape by key column as omics.ID
loadTableData(table = "node", data = dataset, data.key.column = "Identifier", table.key.column = "omics.ID")
```

    ## [1] "Success: Data loaded in defaultnode table"

## Visualization options

``` r
#New visual style is created
RCy3::copyVisualStyle("default", "pathwayStyle")
#Set new style as the current style
RCy3::setVisualStyle("pathwayStyle")
```

    ##                 message 
    ## "Visual Style applied."

``` r
#Set node dimensions as fixed sizes
RCy3::lockNodeDimensions(TRUE, style.name="pathwayStyle")

#Node shape mapping
RCy3::setNodeShapeMapping('Type', c('GeneProduct', 'Protein', 'Metabolite'), c('ELLIPSE', 'ELLIPSE', 'RECTANGLE'), style.name = "pathwayStyle")
```

    ## NULL

``` r
#Change node height
RCy3::setNodeHeightMapping('Type', c('GeneProduct', 'Protein', 'Metabolite'), c(23, 23, 25), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Change node width
RCy3::setNodeWidthMapping('Type', c('GeneProduct', 'Protein', 'Metabolite'), c(60, 60, 100), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Set node color based on log2FC for both genes and metabolites
node.colors <- c(rev(brewer.pal(3, "RdBu")))
setNodeColorMapping("log2FC", c(-1, 0, 1), node.colors, default.color = "#D3D3D3", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Set node border width and color based on p-value
#First we need to get all p-values from node table
pvalues <- getTableColumns(table = 'node', columns = 'pvalue')
pvalues <- na.omit(pvalues)
#Create a range for all sign. p-values, and one for all not significant.
significant_pvalues <- pvalues[(pvalues < 0.05)]
not.significant_pvalues <- pvalues[(pvalues >= 0.05)]
significant_pvalues.colors <- rep("#2e9d1d", length(significant_pvalues))
not.significant_pvalues.colors <- rep("#FFFFFF", length(not.significant_pvalues))

setNodeBorderWidthMapping('pvalue', table.column.values = NULL , c(6, 6) , mapping.type = "c", style.name = "pathwayStyle")
```

    ## NULL

``` r
setNodeBorderColorMapping('pvalue', c(significant_pvalues, not.significant_pvalues), 
                          c(significant_pvalues.colors, not.significant_pvalues.colors), default.color = "#AAAAAA", mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Save output 
filename_multiomics <- paste0("results/", pathway.id, "_", names_of_dataframes, "_omics_visualization.png")
png.file <- file.path(getwd(), filename_multiomics)
exportImage(png.file, 'PNG', zoom = 500)
```

    ##                                                                                                                                         file 
    ## "D:\\bridgeDb\\BridgeDbDemoBioSB2022\\scripts\\5-multi_omics_visualization\\results\\WP4726_Combine_rectum_UC_multi_omics_visualization.png"
