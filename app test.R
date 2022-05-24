library(QFeatures)
# library(conflicted)
library(msdata)
library(tidyverse)
library(limma)
library(msqrob2)
library(plotly)
library(gridExtra)
library(shinyWidgets)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(shiny)
library(shinydashboard)
library(shinyBS)
library(ggplot2)
library(MSnbase)
library(DT)
library(fdrtool)
library(MsCoreUtils)
library(MASS)
library(BiocParallel)
library(missForest)
library(clusterProfiler)
library(DEP)
library(grid)
library(ReactomePA)
library(clusterProfiler.dplyr)
library(tidyr)
library(GO.db)
library(reactome.db)
library(PFAM.db)
library(ComplexHeatmap)
library(reactlog)
# library(esquisse)
source("msqrob function.R")

#shiny.maxRequestSize=30*1024^2 

## ui ##-----------
ui <- dashboardPage(
  
  dashboardHeader(title = "Msqrob"),#** msqrob ui ##------------
  dashboardSidebar(	width = 335,
                    sidebarMenu(
                      convertMenuItem(
                        menuItem("Msqrob", tabName = "msqrob", icon = icon("th"),
                                 menuItem("Files", selected = TRUE,
                                          fluidRow(width = 9,
                                                   h5("Choose the file"),
                                                   fileInput("file",
                                                             "Your file",
                                                             accept=c('text/csv',
                                                                      'text/comma-separated-values,text/plain',
                                                                      '.csv')),
                                                   h5("Choose the type of file"),
                                                   radioButtons("file_type","peptides or protein groups",
                                                                choices=c("peptides","protein groups"),
                                                                selected = "peptides"))
                                 ),	
                                 menuItem("Columns",
                                          selectizeInput("intensity_colmun",
                                                         "Type of Intensity",
                                                         choices=c("LFQ","Intensity"),
                                                         selected = "LFQ"),
                                          uiOutput("ID_colmun"),
                                          uiOutput("Genename_colmun"),
                                          uiOutput("Control"),
                                          downloadButton("download","download report")
                                 ),
                                 menuItem("Correction options",
                                          radioButtons("imputation_type",
                                                       "Type of Imputation",
                                                       choices = c("RF", "MinProb","MinDet","bpca","min","MLE","zero","knn","QRILC",
                                                                   "No imputation"),
                                                       selected = "No imputation"),
                                          br(),br(),
                                          selectizeInput("filter_threshold",
                                                         "Filter Condition",
                                                         choices=c(0,1),selected=0),
                                          br(),
                                          radioButtons("fdr_type",
                                                       "Type of FDR correction",
                                                       choices = c("BH", "fdrtool"),
                                                       selected = "BH"),
                                          shinyBS::bsTooltip("fdr_correction", "Choose the method of pvalue adjustment", "right", options = list(container = "body"))
                                 ),
                                 actionButton("analyze", "Analyze"),
                                 shinyBS::bsTooltip("analyze", "Click on it to analyze your data", "right", options = list(container = "body")),
                                 tags$hr(),
                                 tags$style(type="text/css", "#downloadData {background-color:white;color: black;font-family: Source Sans Pro}"),
                                 uiOutput("downloadTable"),
                                 uiOutput("downloadButton")
                        ),"data_tab"),
                      convertMenuItem(#** Genelist tool menuItem ----
                                      menuItem("Genelist tool options", selected = FALSE, tabName = "genelist_tab", icon = icon("th")), "genelist_tab"
                      ),
                      ## annotation menuitem ##---------
                      convertMenuItem(
                        menuItem("Annotation options", selected = TRUE, tabName = "Annotation_tab", icon = icon("th"),
                                 #menuItem("Import from", selected = TRUE, 
                                 checkboxInput("import_from_genelist_tool_for_annotation", "Import from gene list tool", value = FALSE),
                                 #),
                                 #uiOutput("import_for_annotation"),
                                 #uiOutput("import_contrast_for_annotation"),
                                 # fluidRow(column(width = 5, checkboxInput("From_LFQ", "DEP-LFQ", value = FALSE)),
                                 #          column(width = 6, checkboxInput("From_RNAseq", "DEG-RNAseq", value = FALSE))),
                                 #uiOutput("text_input_for_annotation"),
                                 #uiOutput("organism_for_annotation"),
                                 uiOutput("genelist_tool_for_annotation"),
                                 uiOutput("text_input_for_annotation"),
                                 selectizeInput("organism_for_annotation",
                                                "Select organism",
                                                choices=c("human","mouse"),
                                                selected="human"),
                                 actionButton("analyze_for_annotation", "Analyze"),
                                 tags$hr(),
                                 tags$style(type="text/css", "#downloadannotation {background-color:white;color: black;font-family: Source Sans Pro}"),
                                 uiOutput("downloadannotation_for_output"),
                                 #shinyBS::bsTooltip("import_from_for_annotation", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                 #shinyBS::bsTooltip("import_for_annotation", "Choose genes that you want to do the gene annotation. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
                                 #shinyBS::bsTooltip("import_contrast_for_annotation", "The contrast that you want to do the gene annotation. Note that: [Any significant] represents the genes that are significant in at least one contrast", "right", options = list(container = "body")),          
                                 shinyBS::bsTooltip("text_input_for_annotation", "Paste gene list here", "right", options = list(container = "body")),
                                 shinyBS::bsTooltip("organism_for_annotation", "Select the organism", "right", options = list(container = "body")),
                                 shinyBS::bsTooltip("analyze_for_annotation", "Click on it to analyze", "right", options = list(container = "body")),
                                 shinyBS::bsTooltip("downloadannotation_for_output", "Click on it to download the annotation table", "right", options = list(container = "body"))
                        ), "Annotation_tab" 
                      ),
                      #** ORA menuItem ----
                      menuItem("ORA options", selected = TRUE, icon = icon("th"),
                               convertMenuItem(# GO analysis ----
                                               menuItem("GO options", selected = TRUE, tabName = "GO_tab", icon = icon("angle-double-right"),#bookmark bars book-open
                                                        # menuItem("Import from", selected = TRUE, 
                                                        checkboxInput("import_from_genelist_tool_for_go", "Import from gene list tool", value = FALSE),
                                                        # ),
                                                        # uiOutput("import_for_go"),
                                                        uiOutput("genelist_tool_for_go"),
                                                        uiOutput("text_input_for_go"),
                                                        # uiOutput("organism_for_go"),
                                                        # uiOutput("df_with_lg2fc"),
                                                        # textAreaInput(inputId = "text_input_for_go", label = "Please paste your gene list", placeholder = "TP53\nSOX2", rows = 8, width = "100%"),
                                                        selectizeInput("organism_for_go",
                                                                       "Select organism",
                                                                       choices=c("Human","Mouse"),
                                                                       selected="Human"),
                                                        checkboxInput("df_with_lg2fc", "If with log2 fold change", value = FALSE),
                                                        radioButtons("go_color",
                                                                     "colorBy",
                                                                     c("pvalue", "p.adjust"),
                                                                     selected = "p.adjust"),
                                                        actionButton("analyze_for_go", "Analyze"),
                                                        tags$hr(),
                                                        tags$style(type="text/css", "#downloadgo {background-color:white;color: black;font-family: Source Sans Pro}"),
                                                        uiOutput("downloadTable_go"),
                                                        uiOutput("downloadButton_go"),
                                                        # downloadButton("downloadgo", "Save table")
                                                        shinyBS::bsTooltip("import_from_for_go", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("import_for_go", "Choose genes that you want to do the GO analysis. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("import_contrast_for_go", "The contrast that you want to do the GO analysis. Note that: [Any significant] represents the genes that are significant in at least one contrast, and when you select [Any significant], you should not choose [If with log2 fold change]", "right", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("text_input_for_go", "Paste your gene list here", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("organism_for_go", "Select the organism", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("df_with_lg2fc", "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("analyze_for_go", "Click on it to analyze", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadTable_go", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadButton_go", "Click on it to download the result table", "right", options = list(container = "body"))
                                               ), "GO_tab"           
                               ),
                               convertMenuItem(#** KEGG analysis ----
                                               menuItem("KEGG options", selected = TRUE, tabName = "KEGG_tab", icon = icon("angle-double-right"),
                                                        checkboxInput("import_from_genelist_tool_for_kegg", "Import from gene list tool", value = FALSE),
                                                        uiOutput("genelist_tool_for_kegg"),
                                                        uiOutput("text_input_for_kegg"),
                                                        selectizeInput("organism_for_kegg",
                                                                       "Select organism",
                                                                       choices=c("Human","Mouse"),
                                                                       selected="Human"),
                                                        checkboxInput("df_with_lg2fc_for_kegg", "If with log2 fold change", value = FALSE),
                                                        radioButtons("kegg_color",
                                                                     "colorBy",
                                                                     c("pvalue", "p.adjust"),
                                                                     selected = "p.adjust"),
                                                        # menuItem("Import from", selected = TRUE, 
                                                        #          checkboxInput("import_from_for_kegg", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                                                        # ),
                                                        # uiOutput("import_for_kegg"),
                                                        # uiOutput("import_contrast_for_kegg"),
                                                        # uiOutput("text_input_for_kegg"),
                                                        # uiOutput("organism_for_kegg"),
                                                        # uiOutput("df_with_lg2fc_for_kegg"),
                                                        actionButton("analyze_for_kegg", "Analyze"),
                                                        tags$hr(),
                                                        tags$style(type="text/css", "#downloadkegg {background-color:white;color: black;font-family: Source Sans Pro}"),
                                                        uiOutput("downloadTable_kegg"),
                                                        uiOutput("downloadButton_kegg"),
                                                        # downloadButton("downloadkegg", "Save table")
                                                        shinyBS::bsTooltip("import_from_for_kegg", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("import_for_kegg", "Choose genes that you want to do the KEGG analysis. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("import_contrast_for_kegg", "The contrast that you want to do the KEGG analysis. Note that: [Any significant] represents the genes that are significant in at least one contrast, and when you select [Any significant], you should not choose [If with log2 fold change]", "right", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("text_input_for_kegg", "Paste your gene list here", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("organism_for_kegg", "Select the organism: hsa represents human, mmu represents mouse", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("df_with_lg2fc_for_kegg", "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("analyze_for_kegg", "Click on it to analyze", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadTable_kegg", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadButton_kegg", "Click on it to download the result table", "right", options = list(container = "body"))
                                               ), "KEGG_tab"      
                               ),
                               convertMenuItem(#** Reactome analysis ----
                                               menuItem("Reactome options", selected = TRUE, tabName = "Reactome_tab", icon = icon("angle-double-right"),
                                                        checkboxInput(inputId = "import_from_genelist_tool_for_reactome", "Import from gene list tool", value = FALSE),
                                                        uiOutput("genelist_tool_for_reactome"),
                                                        uiOutput("text_input_for_reactome"),
                                                        selectizeInput("organism_for_reactome",
                                                                       "Select organism",
                                                                       choices=c("Human","Mouse"),
                                                                       selected="Human"),
                                                        checkboxInput("df_with_lg2fc_for_reactome", "If with log2 fold change", value = FALSE),
                                                        radioButtons("reactome_color",
                                                                     "colorBy",
                                                                     c("pvalue", "p.adjust"),
                                                                     selected = "p.adjust"),
                                                        # menuItem("Import from", selected = TRUE, 
                                                        #          checkboxInput("import_from_for_reactome", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                                                        # ),
                                                        # uiOutput("import_for_reactome"),
                                                        # uiOutput("import_contrast_for_reactome"),
                                                        # uiOutput("text_input_for_reactome"),
                                                        # uiOutput("organism_for_reactome"),
                                                        # uiOutput("df_with_lg2fc_for_reactome"),
                                                        actionButton("analyze_for_reactome", "Analyze"),
                                                        tags$hr(),
                                                        tags$style(type="text/css", "#downloadreactome {background-color:white;color: black;font-family: Source Sans Pro}"),
                                                        uiOutput("downloadTable_reactome"),
                                                        uiOutput("downloadButton_reactome"),
                                                        # downloadButton("downloadreactome", "Save table")
                                                        shinyBS::bsTooltip("import_from_for_reactome", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("import_for_reactome", "Choose genes that you want to do the reactome analysis. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("import_contrast_for_reactome", "The contrast that you want to do the reactome analysis. Note that: [Any significant] represents the genes that are significant in at least one contrast, and when you select [Any significant], you should not choose [If with log2 fold change]", "right", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("text_input_for_reactome", "Paste your gene list here", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("organism_for_reactome", "Select the organism", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("df_with_lg2fc_for_reactome", "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("analyze_for_reactome", "Click on it to analyze", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadTable_reactome", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadButton_reactome", "Click on it to download the result table", "right", options = list(container = "body"))
                                               ), "Reactome_tab"            
                               )
                      ),
                      menuItem("GSEA options", selected = TRUE, icon = icon("th"),
                               convertMenuItem(#** gseGO analysis----
                                               menuItem("gseGO options", selected = TRUE, tabName = "gseGO_tab", icon = icon("angle-double-right"),
                                                        checkboxInput("import_from_genelist_tool_for_gsego", "DEP-LFQ or DEG-RNAseq", value = FALSE),
                                                        uiOutput("genelist_tool_for_gsego"),
                                                        uiOutput("text_input_for_gsego"),
                                                        selectizeInput("organism_for_gsego",
                                                                       "Select organism",
                                                                       choices=c("Human","Mouse"),
                                                                       selected="Human"),
                                                        radioButtons("gsego_color",
                                                                     "colorBy",
                                                                     c("pvalue", "p.adjust"),
                                                                     selected = "p.adjust"),
                                                        # menuItem("Import from", selected = TRUE, 
                                                        #          checkboxInput("import_from_for_gsego", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                                                        # ),
                                                        # uiOutput("import_for_gsego"),
                                                        # uiOutput("import_contrast_for_gsego"),
                                                        
                                                        # uiOutput("organism_for_gsego"),
                                                        actionButton("analyze_for_gsego", "Analyze"),
                                                        tags$hr(),
                                                        tags$style(type="text/css", "#downloadgsego {background-color:white;color: black;font-family: Source Sans Pro}"),
                                                        uiOutput("downloadTable_gsego"),
                                                        uiOutput("downloadButton_gsego"),
                                                        # downloadButton("downloadgsego", "Save table")
                                                        shinyBS::bsTooltip("import_from_for_gsego", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("import_for_gsego", "Choose genes that you want to do the gseGO analysis. DEP-LFQ: import from DEP-LFQ panel; DEG-RNAseq: import from DEG-RNAseq panel", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("import_contrast_for_gsego", "The contrast that you want to do the gseGO analysis", "right", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("text_input_for_gsego", "Paste your gene list with log2 fold change here", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("organism_for_gsego", "Select the organism", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("analyze_for_gsego", "Click on it to analyze", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadTable_gsego", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadButton_gsego", "Click on it to download the result table", "right", options = list(container = "body"))
                                               ), "gseGO_tab"            
                               ),
                               convertMenuItem(#** gseKEGG analysis----
                                               menuItem("gseKEGG options", selected = TRUE, tabName = "gseKEGG_tab", icon = icon("angle-double-right"),
                                                        checkboxInput("import_from_genelist_tool_for_gsekegg", "DEP-LFQ or DEG-RNAseq", value = FALSE),
                                                        uiOutput("genelist_tool_for_gsekegg"),
                                                        uiOutput("text_input_for_gsekegg"),
                                                        selectizeInput("organism_for_gsekegg",
                                                                       "Select organism",
                                                                       choices=c("Human","Mouse"),
                                                                       selected="Human"),
                                                        radioButtons("gsekegg_color",
                                                                     "colorBy",
                                                                     c("pvalue", "p.adjust"),
                                                                     selected = "p.adjust"),
                                                        # menuItem("Import from", selected = TRUE, 
                                                        #          checkboxInput("import_from_for_gsekegg", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                                                        # ),
                                                        # uiOutput("import_for_gsekegg"),
                                                        # uiOutput("import_contrast_for_gsekegg"),
                                                        # uiOutput("text_input_for_gsekegg"),
                                                        # uiOutput("organism_for_gsekegg"),
                                                        actionButton("analyze_for_gsekegg", "Analyze"),
                                                        tags$hr(),
                                                        tags$style(type="text/css", "#downloadgsekegg {background-color:white;color: black;font-family: Source Sans Pro}"),
                                                        uiOutput("downloadTable_gsekegg"),
                                                        uiOutput("downloadButton_gsekegg"),
                                                        # downloadButton("downloadgsekegg", "Save table")
                                                        shinyBS::bsTooltip("import_from_for_gsekegg", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("import_for_gsekegg", "Choose genes that you want to do the gseKEGG analysis. DEP-LFQ: import from DEP-LFQ panel; DEG-RNAseq: import from DEG-RNAseq panel", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("import_contrast_for_gsekegg", "The contrast that you want to do the gseKEGG analysis", "right", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("text_input_for_gsekegg", "Paste your gene list with log2 fold change here", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("organism_for_gsekegg", "Select the organism: hsa represents human, mmu represents mouse", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("analyze_for_gsekegg", "Click on it to analyze", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadTable_gsekegg", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadButton_gsekegg", "Click on it to download the result table", "right", options = list(container = "body"))
                                               ), "gseKEGG_tab"
                               ),
                               convertMenuItem(#** gseReactome analysis----
                                               menuItem("gseReactome options", selected = TRUE, tabName = "gseReactome_tab", icon = icon("angle-double-right"),
                                                        checkboxInput("import_from_genelist_tool_for_gsereactome", "DEP-LFQ or DEG-RNAseq", value = FALSE),
                                                        uiOutput("genelist_tool_for_gsereactome"),
                                                        uiOutput("text_input_for_gsereactome"),
                                                        selectizeInput("organism_for_gsereactome",
                                                                       "Select organism",
                                                                       choices=c("Human","Mouse"),
                                                                       selected="Human"),
                                                        radioButtons("gsereactome_color",
                                                                     "colorBy",
                                                                     c("pvalue", "p.adjust"),
                                                                     selected = "p.adjust"),
                                                        # menuItem("Import from", selected = TRUE, 
                                                        #          checkboxInput("import_from_for_gsereactome", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                                                        # ),
                                                        # uiOutput("import_for_gsereactome"),
                                                        # uiOutput("import_contrast_for_gsereactome"),
                                                        # uiOutput("text_input_for_gsereactome"),
                                                        # uiOutput("organism_for_gsereactome"),
                                                        actionButton("analyze_for_gsereactome", "Analyze"),
                                                        tags$hr(),
                                                        tags$style(type="text/css", "#downloadgsereactome {background-color:white;color: black;font-family: Source Sans Pro}"),
                                                        uiOutput("downloadTable_gsereactome"),
                                                        uiOutput("downloadButton_gsereactome"),
                                                        # downloadButton("downloadgsereactome", "Save table")
                                                        shinyBS::bsTooltip("import_from_for_gsereactome", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("import_for_gsereactome", "Choose genes that you want to do the gseReactome analysis. DEP-LFQ: import from DEP-LFQ panel; DEG-RNAseq: import from DEG-RNAseq panel", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("import_contrast_for_gsereactome", "The contrast that you want to do the gseReactome analysis", "right", options = list(container = "body")),          
                                                        shinyBS::bsTooltip("text_input_for_gsereactome", "Paste your gene list with log2 fold change here", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("organism_for_gsereactome", "Select the organism", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("analyze_for_gsereactome", "Click on it to analyze", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadTable_gsereactome", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                                                        shinyBS::bsTooltip("downloadButton_gsereactome", "Click on it to download the result table", "right", options = list(container = "body"))
                                               ), "gseReactome_tab"           
                               )
                      )
                    )),	
  ## dashboardBody ##---------------
  dashboardBody(	
    tabItems( ## tabItem msqrob ##--------
      tabItem(tabName = "msqrob",
              fluidRow(
                box(numericInput("padj","adj. P value",
                                 min = 0.0001, max = 0.1, value = 0.05),width = 2),
                box(numericInput("logFC","Log2 fold change",
                                 min = 0, max = 10, value = 1),width = 2),
                infoBoxOutput("significantBox")
              ),
              fluidRow(
                column(width = 7,
                       box(title = "Top Table",
                           # selectizeInput("table_sig","Show Significant protein",choices=c("TRUE", "FALSE"),width = 12),
                           box(uiOutput("select"), width = 6),
                           # box(uiOutput("exclude"), width = 6),
                           DT::dataTableOutput("table"), width = 12)
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Heatmap",
                                       fluidRow(width = 6,
                                                box(
                                                  selectizeInput("colorbar","colorbar",
                                                                 choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                                                 selected = c("RdBu"), multiple = FALSE), width = 6),
                                                box(
                                                  selectizeInput("heatmaptype","heatmaptype",
                                                                   choices=c("contrast", "centered"),selected=c("centered")
                                                                   ),width = 6),
                                                box(
                                                  selectizeInput("heatmapcluster","heatmapcluster",
                                                                   choices=c("euclidean", "maximum", "manhattan", "canberra",
                                                                             "binary", "minkowski", "pearson", "spearman", 
                                                                             "kendall", "gower"),selected=c("euclidean")
                                                                   ),width = 6),
                                                box(
                                                  numericInput("col_limit","Col limit",
                                                               min = 0, max = 10, value = 4),width = 6)
                                       ),
                                       fluidRow(width = 6,
                                                box(numericInput("Heatmap_Width",
                                                             "width",
                                                             min = 1, max = 30, value = 6)),
                                                box(numericInput("Heatmap_Height",
                                                             "height",
                                                             min = 1, max = 30, value = 6))),
                                       fluidRow(
                                         box(plotOutput("heatmap"),inline=T,width = 12),
                                         downloadButton('downloadHeatmap', 'Save heatmap')
                                       )
                                       ),
                              tabPanel(title = "Volcano plot",
                                       fluidRow(
                                         box(uiOutput("volcano_cond"), width = 9),
                                         box(numericInput("fontsize",
                                                          "Font size",
                                                          min = 0, max = 8, value = 4),
                                             
                                             width = 3)
                                       ), 
                                       fluidRow(
                                         box(numericInput("Volcano_Width",
                                                          "width",
                                                          min = 1, max = 30, value = 7),
                                             numericInput("Volcano_Height",
                                                          "height",
                                                          min = 1, max = 30, value = 7),
                                             #selectizeInput("mybreaks","Mybreaks",choices = seq(-40, 40, 1), multiple = TRUE, size = 10),
                                             width = 8, collapsible = TRUE, collapsed = TRUE)),
                                       fluidRow(
                                         plotOutput("volcano", height = 600),
                                         downloadButton('downloadVolcano', 'Save volcano')
                                       )
                              )
                       ),
                       tabBox(title = "QC Plots", width = 12,
                              tabPanel(title = "MDS plot",
                                       fluidRow(width = 6,
                                         box(numericInput("MDStop","Top proteins",
                                                          500,min = 1,max = 1000), width = 6),
                                         box(selectInput("MDStype","Peptide or Protein",
                                                          choices = c("peptide","protein"),
                                                         selected = "protein"), width = 6),
                                         box(numericInput("MDS_Width",
                                                            "width",
                                                            min = 1, max = 30, value = 7), width = 6),
                                         box(numericInput("MDS_Height",
                                                            "height",
                                                            min = 1, max = 30, value = 7), width = 6)
                                         ),
                                         fluidRow(
                                           plotOutput("MDS"),
                                           downloadButton('downloadMDS', 'Save')
                                         )
                                       ),
                              tabPanel(title = "PCA plot",
                                         fluidRow(
                                           box(numericInput("PCAtop",
                                                            "Top proteins",
                                                            500,min = 1,max = 1000), width = 6),
                                           # box(selectizeInput("Indicate",
                                           #                    "Color and shape",
                                           #                    choices = c("condition", "replicate", "Condition", "Replicate"),
                                           #                    selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                             box(numericInput("PCA_width",
                                                              "width",
                                                              min = 1, max = 30, value = 7), width = 6),
                                             box(numericInput("PCA_height",
                                                              "height",
                                                              min = 1, max = 30, value = 7), width = 6)
                                           ),
                                           fluidRow(
                                             plotOutput("pca"),
                                             downloadButton('downloadPCA', 'Save')
                                           )
                                         
                              ),
                              tabPanel(title = "Pearson correlation",
                                       fluidRow(
                                         box(selectizeInput("Pearson_pal",
                                                            "color panel",
                                                            choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                            selected = c("PRGn"), multiple = FALSE), 
                                             checkboxInput("Pearson_pal_rev",
                                                           "pal rev",
                                                           value = FALSE), width = 4),
                                         box(uiOutput("cor_condition")),
                                         box(numericInput("Pearson_Width",
                                                          "width",
                                                          min = 1, max = 30, value = 7), 
                                             numericInput("Pearson_Height",
                                                          "height",
                                                          min = 1, max = 30, value = 7), width = 4),
                                         box(numericInput("Pearson_lower",
                                                          "lower",
                                                          min = -1, max = 1, value = -1), 
                                             numericInput("Pearson_upper",
                                                          "upper",
                                                          min = -1, max = 1, value = 1), width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("add_values_for_DEP_person",
                                                           "Add values",
                                                           value = FALSE), 
                                             numericInput("value_size_for_DEP_person",
                                                          "Value size",
                                                          min = 1, max = 30, value = 10),
                                             numericInput("value_digits_for_DEP_person",
                                                          "Value digits",
                                                          min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                       ),
                                       fluidRow(
                                         plotOutput("cor", height = 600),
                                         downloadButton('download_Pearson_correlation', 'Save')
                                       ),
                                       shinyBS::bsTooltip("Pearson_pal", "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("Pearson_pal_rev", "Whether or not to invert the color palette", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("Pearson_Width", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("Pearson_Height", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("Pearson_lower", "Set the lower limit of the color scale", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("Pearson_upper", "Set the upper limit of the color scale", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Sample CVs",
                                       fluidRow(
                                         box(numericInput("Sample_CVs_Width",
                                                          "width",
                                                          min = 1, max = 30, value = 7),
                                             width = 6),
                                         box(
                                           numericInput("Sample_CVs_Height",
                                                        "height",
                                                        min = 1, max = 30, value = 7), width = 6
                                         )
                                       ),
                                       fluidRow(
                                         plotOutput("CVs", height = 600),
                                         downloadButton('download_Sample_CVs', 'Save')
                                       ),
                                       shinyBS::bsTooltip("Sample_CVs_Width", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("Sample_CVs_Height", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                )
              )	
      ),
      tabItem(#** genelist tabItem----
              tabName = "genelist_tab",
              Genelist_UImodule("bbb") 
      ),
      ## Annotation tabItem ##----
      tabItem(
        tabName = "Annotation_tab",
        fluidRow(
          box(title = "Annotation Table",
              DT::dataTableOutput("annotation_Table"), width = 12)
        ),
        shinyBS::bsTooltip("annotation_Table", "The table of annotation result", "top", options = list(container = "body"))
        
      ),
      
      ## tabItem GO ##----------
      tabItem(tabName = "GO_tab",
              fluidRow(
                box(selectizeInput("go_ont",
                                   "Ontology",
                                   choices = c("BP", "CC", "MF", "ALL"),
                                   selected = c("BP"), multiple = FALSE), 
                    width = 2),
                box(numericInput("go_p",
                                 "P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("go_padj",
                                 "adj. P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("go_qvalue",
                                 "q value",
                                 min = 0, max = 1, value = 0.2),
                    width = 2),
                # box(radioButtons("go_color",
                #                  "colorBy",
                #                  c("pvalue", "p.adjust"),
                #                  selected = "p.adjust"), width = 4),
                # infoBoxOutput("significantBox_for_go", width = 4)
                shinyBS::bsTooltip("go_ont", "One of BP, MF, and CC subontologies, or ALL for all three", "top", options = list(container = "body")),
                shinyBS::bsTooltip("go_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("go_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("go_qvalue", "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("go_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       fluidRow(
                         box(checkboxInput("go_simplify", "removed redundancy of enriched GO terms", value = FALSE), width = 5),
                         # box(radioButtons("go_color",
                         #          "colorBy",
                         #          c("pvalue", "p.adjust"),
                         #          selected = "p.adjust"), width = 4),
                         infoBoxOutput("significantBox_for_go", width = 7)
                       ),
                       fluidRow(
                         box(title = "GO Top Table",
                             DT::dataTableOutput("go_Table"), width = 12)
                       ),
                       shinyBS::bsTooltip("go_simplify", "Whether simplify output by removing redundancy of enriched GO terms", "top", options = list(container = "body")),
                       shinyBS::bsTooltip("go_Table", "The table of the significant result", "top", options = list(container = "body"))
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Bar plot",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("go_wide_bar", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("go_high_bar", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("go_bar_if_Split_for_ont_ALL",
                                                           "If splited by ontology",
                                                           value = FALSE), width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("go_barplot", height = 500),
                                         downloadButton('download_go_barplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_bar", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_bar_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each BP, CC and MF", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Dot plot",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("go_wide_dot", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("go_high_dot", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("go_dot_if_Split_for_ont_ALL",
                                                           "If splited by ontology",
                                                           value = FALSE), width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("go_dotplot", height = 500),
                                         downloadButton('download_go_dotplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_dot", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_dot_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each BP, CC and MF", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Heatplot",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("go_wide_heat", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("go_high_heat", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("go_heatplot", height = 500),
                                         downloadButton('download_go_heatplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Cnetplot",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                             width = 4),
                                         box(numericInput("go_wide_cnet", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("go_high_cnet", "Height", min = 1, max = 50, value = 6),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("go_circular_cnet",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("go_cnetplot", height = 500),
                                         downloadButton('download_go_cnetplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Emaplot",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("go_wide_ema", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("go_high_ema", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("go_emaplot", height = 500),
                                         downloadButton('download_go_emaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Goplot",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_go", "ShowCategory", min = 1, max = 100, value = 10),
                                             width = 4),
                                         box(numericInput("go_wide_go", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("go_high_go", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("go_circular_go",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("go_goplot", height = 500),
                                         downloadButton('download_go_goplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_go", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_go", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_go", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_circular_go", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "GOgraph",
                                       fluidRow(
                                         box(numericInput("go_ShowCategory_GOgraph", "FirstSigNodes", min = 1, max = 100, value = 10),
                                             width = 4),
                                         box(numericInput("go_wide_GOgraph", "Width", min = 1, max = 50, value = 7),
                                             width = 4),
                                         box(numericInput("go_high_GOgraph", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(actionButton("plot_for_GOgraphplot", "Plot", width = "10%", icon = icon("caret-right"))),#step-forward
                                       fluidRow(
                                         plotOutput("go_GOgraphplot", height = 500),
                                         downloadButton('download_go_GOgraphplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("go_ShowCategory_GOgraph", "Number of significant nodes (retangle nodes in the graph)", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_wide_GOgraph", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("go_high_GOgraph", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                       
                )
              )
      ),
      ## tabItem kegg ##----------
      tabItem(tabName = "KEGG_tab",
              fluidRow(
                box(numericInput("kegg_p",
                                 "P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("kegg_padj",
                                 "adj. P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("kegg_qvalue",
                                 "q value",
                                 min = 0, max = 1, value = 0.2),
                    width = 2),
                infoBoxOutput("significantBox_for_kegg", width = 4),
                # box(radioButtons("kegg_color",
                #                  "colorBy",
                #                  c("pvalue", "p.adjust"),
                #                  selected = "p.adjust"), width = 2),
                shinyBS::bsTooltip("kegg_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("kegg_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("kegg_qvalue", "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("kegg_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       fluidRow(
                         box(title = "KEGG Top Table",
                             DT::dataTableOutput("kegg_Table"), width = 12)
                       ),
                       shinyBS::bsTooltip("kegg_Table", "The table of the significant result", "top", options = list(container = "body"))
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Bar plot",
                                       fluidRow(
                                         box(numericInput("kegg_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("kegg_wide_bar", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("kegg_high_bar", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("kegg_barplot", height = 500),
                                         downloadButton('download_kegg_barplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("kegg_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Dot plot",
                                       fluidRow(
                                         box(numericInput("kegg_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("kegg_wide_dot", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("kegg_high_dot", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("kegg_dotplot", height = 500),
                                         downloadButton('download_kegg_dotplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("kegg_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Heatplot",
                                       fluidRow(
                                         box(numericInput("kegg_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("kegg_wide_heat", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("kegg_high_heat", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("kegg_heatplot", height = 500),
                                         downloadButton('download_kegg_heatplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("kegg_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Cnetplot",
                                       fluidRow(
                                         box(numericInput("kegg_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                             width = 4),
                                         box(numericInput("kegg_wide_cnet", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("kegg_high_cnet", "Height", min = 1, max = 50, value = 6),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("kegg_circular_cnet",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("kegg_cnetplot", height = 500),
                                         downloadButton('download_kegg_cnetplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("kegg_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Emaplot",
                                       fluidRow(
                                         box(numericInput("kegg_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("kegg_wide_ema", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("kegg_high_ema", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("kegg_emaplot", height = 500),
                                         downloadButton('download_kegg_emaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("kegg_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("kegg_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                       
                )
              )
      ),
      ## tabItem reactome ##----------
      tabItem(tabName = "Reactome_tab",
              #         tags$head(tags$style(type="text/css", "
              #                  #loadmessage {
              #                  top: 0px; left: 0px;
              #                  width: 100%; padding: 5px 0px 5px 0px;
              #                  text-align: center; font-weight: bold;
              #                  font-size: 100%; color: #000000;
              #                  background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
              # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
              #                  tags$div("calculating...please wait...",id="loadmessage")),
              fluidRow(
                box(numericInput("reactome_p",
                                 "P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("reactome_padj",
                                 "adj. P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("reactome_qvalue",
                                 "q value",
                                 min = 0, max = 1, value = 0.2),
                    width = 2),
                infoBoxOutput("significantBox_for_reactome", width = 4),
                # box(radioButtons("reactome_color",
                #                  "colorBy",
                #                  c("pvalue", "p.adjust"),
                #                  selected = "p.adjust"), width = 2),
                shinyBS::bsTooltip("reactome_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("reactome_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("reactome_qvalue", "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("reactome_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       fluidRow(
                         box(title = "Reactome Top Table",
                             DT::dataTableOutput("reactome_Table"), width = 12),
                         shinyBS::bsTooltip("reactome_Table", "The table of the significant result", "top", options = list(container = "body"))
                       )
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Bar plot",
                                       fluidRow(
                                         box(numericInput("reactome_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("reactome_wide_bar", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("reactome_high_bar", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("reactome_barplot", height = 500),
                                         downloadButton('download_reactome_barplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("reactome_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Dot plot",
                                       fluidRow(
                                         box(numericInput("reactome_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("reactome_wide_dot", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("reactome_high_dot", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("reactome_dotplot", height = 500),
                                         downloadButton('download_reactome_dotplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("reactome_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Heatplot",
                                       fluidRow(
                                         box(numericInput("reactome_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("reactome_wide_heat", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("reactome_high_heat", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("reactome_heatplot", height = 500),
                                         downloadButton('download_reactome_heatplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("reactome_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Cnetplot",
                                       fluidRow(
                                         box(numericInput("reactome_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                             width = 4),
                                         box(numericInput("reactome_wide_cnet", "Width", min = 1, max = 50, value = 10),
                                             width = 4),
                                         box(numericInput("reactome_high_cnet", "Height", min = 1, max = 50, value = 6),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("reactome_circular_cnet",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("reactome_cnetplot", height = 500),
                                         downloadButton('download_reactome_cnetplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("reactome_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Emaplot",
                                       fluidRow(
                                         box(numericInput("reactome_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("reactome_wide_ema", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("reactome_high_ema", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("reactome_emaplot", height = 500),
                                         downloadButton('download_reactome_emaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("reactome_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("reactome_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                       
                )
              )
      ),
      
      ## tabitem gseGO##----
      tabItem(tabName = "gseGO_tab",
              fluidRow(
                box(selectizeInput("gsego_ont",
                                   "Ontology",
                                   choices = c("BP", "CC", "MF", "ALL"),
                                   selected = c("BP"), multiple = FALSE), 
                    width = 2),
                box(numericInput("gsego_p",
                                 "P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("gsego_padj",
                                 "adj. P value",
                                 min = 0, max = 1, value = 0.25),
                    width = 2),
                box(numericInput("gsego_NES",
                                 "NES",
                                 min = 0, max = 10, value = 1),
                    width = 2),
                box(selectizeInput("gsego_Phenotype",
                                   "Phenotype",
                                   choices = c("activated", "suppressed"),
                                   selected = c("activated", "suppressed"), multiple = TRUE), width = 2),
                # box(radioButtons("gsego_color",
                #                  "colorBy",
                #                  c("pvalue", "p.adjust"),
                #                  selected = "p.adjust"), width = 2),
                shinyBS::bsTooltip("gsego_ont", "One of BP, MF, and CC subontologies, or ALL for all three", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsego_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsego_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsego_NES", "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsego_Phenotype", "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsego_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       fluidRow(
                         box(checkboxInput("gsego_simplify", "removed redundancy of enriched GO terms", value = FALSE), width = 5),
                         # box(radioButtons("go_color",
                         #          "colorBy",
                         #          c("pvalue", "p.adjust"),
                         #          selected = "p.adjust"), width = 4),
                         infoBoxOutput("significantBox_for_gsego", width = 7)
                       ),
                       fluidRow(
                         box(title = "gseGO Top Table",
                             DT::dataTableOutput("gsego_Table"), width = 12)
                       ),
                       shinyBS::bsTooltip("gsego_simplify", "Whether simplify output by removing redundancy of enriched GO terms", "top", options = list(container = "body")),
                       shinyBS::bsTooltip("gsego_Table", "The table of the significant result", "top", options = list(container = "body"))
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Bar plot",
                                       fluidRow(
                                         box(numericInput("gsego_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("gsego_wide_bar", "Width", min = 1, max = 50, value = 12),
                                             width = 4),
                                         box(numericInput("gsego_high_bar", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("gsego_bar_if_Split_for_ont_ALL",
                                                           "If splited by ontology",
                                                           value = FALSE), width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsego_barplot", height = 500),
                                         downloadButton('download_gsego_barplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsego_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_high_bar", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_bar_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each phenotype that you selected based on each BP, CC and MF", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Dot plot",
                                       fluidRow(
                                         box(numericInput("gsego_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("gsego_wide_dot", "Width", min = 1, max = 50, value = 12),
                                             width = 4),
                                         box(numericInput("gsego_high_dot", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("gsego_dot_if_Split_for_ont_ALL",
                                                           "If splited by ontology",
                                                           value = FALSE), width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsego_dotplot", height = 500),
                                         downloadButton('download_gsego_dotplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsego_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_high_dot", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_dot_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each phenotype that you selected based on each BP, CC and MF", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Heatplot",
                                       fluidRow(
                                         box(numericInput("gsego_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("gsego_wide_heat", "Width", min = 1, max = 200, value = 80),
                                             width = 4),
                                         box(numericInput("gsego_high_heat", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsego_heatplot", height = 500),
                                         downloadButton('download_gsego_heatplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsego_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Cnetplot",
                                       fluidRow(
                                         box(numericInput("gsego_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                             width = 4),
                                         box(numericInput("gsego_wide_cnet", "Width", min = 1, max = 200, value = 30),
                                             width = 4),
                                         box(numericInput("gsego_high_cnet", "Height", min = 1, max = 200, value = 25),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("gsego_circular_cnet",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsego_cnetplot", height = 500),
                                         downloadButton('download_gsego_cnetplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsego_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Emaplot",
                                       fluidRow(
                                         box(numericInput("gsego_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("gsego_wide_ema", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("gsego_high_ema", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsego_emaplot", height = 500),
                                         downloadButton('download_gsego_emaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsego_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Gseaplot",
                                       fluidRow(box(uiOutput("gsego_term"), width = 12)),
                                       fluidRow(
                                         box(numericInput("gsego_wide_Gsea", "Width", min = 1, max = 50, value = 8),
                                             width = 4),
                                         box(numericInput("gsego_high_Gsea", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsego_Gseaplot", height = 500),
                                         downloadButton('download_gsego_Gseaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsego_term", "The term that you want to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_wide_Gsea", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsego_high_Gsea", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                       
                )
              )
      ),
      
      ## tabitem gseKEGG##------
      tabItem(tabName = "gseKEGG_tab",
              fluidRow(
                box(numericInput("gsekegg_p",
                                 "P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("gsekegg_padj",
                                 "adj. P value",
                                 min = 0, max = 1, value = 0.25),
                    width = 2),
                box(numericInput("gsekegg_NES",
                                 "NES",
                                 min = 0, max = 10, value = 1),
                    width = 2),
                box(selectizeInput("gsekegg_Phenotype",
                                   "Phenotype",
                                   choices = c("activated", "suppressed"),
                                   selected = c("activated", "suppressed"), multiple = TRUE), width = 4),
                # box(radioButtons("gsekegg_color",
                #                  "colorBy",
                #                  c("pvalue", "p.adjust"),
                #                  selected = "p.adjust"), width = 2),
                shinyBS::bsTooltip("gsekegg_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsekegg_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsekegg_NES", "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsekegg_Phenotype", "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsekegg_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       fluidRow(
                         infoBoxOutput("significantBox_for_gsekegg", width = 7)
                       ),
                       fluidRow(
                         box(title = "gseKEGG Top Table",
                             DT::dataTableOutput("gsekegg_Table"), width = 12)
                       ),
                       shinyBS::bsTooltip("gsekegg_Table", "The table of the significant result", "top", options = list(container = "body"))
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Bar plot",
                                       fluidRow(
                                         box(numericInput("gsekegg_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("gsekegg_wide_bar", "Width", min = 1, max = 50, value = 12),
                                             width = 4),
                                         box(numericInput("gsekegg_high_bar", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsekegg_barplot", height = 500),
                                         downloadButton('download_gsekegg_barplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsekegg_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Dot plot",
                                       fluidRow(
                                         box(numericInput("gsekegg_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("gsekegg_wide_dot", "Width", min = 1, max = 50, value = 12),
                                             width = 4),
                                         box(numericInput("gsekegg_high_dot", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsekegg_dotplot", height = 500),
                                         downloadButton('download_gsekegg_dotplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsekegg_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Heatplot",
                                       fluidRow(
                                         box(numericInput("gsekegg_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("gsekegg_wide_heat", "Width", min = 1, max = 200, value = 80),
                                             width = 4),
                                         box(numericInput("gsekegg_high_heat", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsekegg_heatplot", height = 500),
                                         downloadButton('download_gsekegg_heatplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsekegg_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Cnetplot",
                                       fluidRow(
                                         box(numericInput("gsekegg_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                             width = 4),
                                         box(numericInput("gsekegg_wide_cnet", "Width", min = 1, max = 200, value = 30),
                                             width = 4),
                                         box(numericInput("gsekegg_high_cnet", "Height", min = 1, max = 200, value = 25),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("gsekegg_circular_cnet",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsekegg_cnetplot", height = 500),
                                         downloadButton('download_gsekegg_cnetplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsekegg_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Emaplot",
                                       fluidRow(
                                         box(numericInput("gsekegg_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("gsekegg_wide_ema", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("gsekegg_high_ema", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsekegg_emaplot", height = 500),
                                         downloadButton('download_gsekegg_emaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsekegg_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Gseaplot",
                                       fluidRow(box(uiOutput("gsekegg_term"), width = 12)),
                                       fluidRow(
                                         box(numericInput("gsekegg_wide_Gsea", "Width", min = 1, max = 50, value = 8),
                                             width = 4),
                                         box(numericInput("gsekegg_high_Gsea", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsekegg_Gseaplot", height = 500),
                                         downloadButton('download_gsekegg_Gseaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsekegg_term", "The term that you want to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_wide_Gsea", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsekegg_high_Gsea", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                       
                )
              )),
      ## tabitem gsereactome ##------
      tabItem(tabName = "gseReactome_tab",
              fluidRow(
                box(numericInput("gsereactome_p",
                                 "P value",
                                 min = 0, max = 1, value = 0.05),
                    width = 2),
                box(numericInput("gsereactome_padj",
                                 "adj. P value",
                                 min = 0, max = 1, value = 0.25),
                    width = 2),
                box(numericInput("gsereactome_NES",
                                 "NES",
                                 min = 0, max = 10, value = 1),
                    width = 2),
                box(selectizeInput("gsereactome_Phenotype",
                                   "Phenotype",
                                   choices = c("activated", "suppressed"),
                                   selected = c("activated", "suppressed"), multiple = TRUE), width = 4),
                # box(radioButtons("gsereactome_color",
                #                  "colorBy",
                #                  c("pvalue", "p.adjust"),
                #                  selected = "p.adjust"), width = 2),
                shinyBS::bsTooltip("gsereactome_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsereactome_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsereactome_NES", "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsereactome_Phenotype", "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
                shinyBS::bsTooltip("gsereactome_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       fluidRow(
                         infoBoxOutput("significantBox_for_gsereactome", width = 7)
                       ),
                       fluidRow(
                         box(title = "gseReactome Top Table",
                             DT::dataTableOutput("gsereactome_Table"), width = 12)
                       ),
                       shinyBS::bsTooltip("gsereactome_Table", "The table of the significant result", "top", options = list(container = "body"))
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabPanel(title = "Bar plot",
                                       fluidRow(
                                         box(numericInput("gsereactome_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("gsereactome_wide_bar", "Width", min = 1, max = 50, value = 12),
                                             width = 4),
                                         box(numericInput("gsereactome_high_bar", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsereactome_barplot", height = 500),
                                         downloadButton('download_gsereactome_barplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsereactome_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Dot plot",
                                       fluidRow(
                                         box(numericInput("gsereactome_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                             width = 4),
                                         box(numericInput("gsereactome_wide_dot", "Width", min = 1, max = 50, value = 12),
                                             width = 4),
                                         box(numericInput("gsereactome_high_dot", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsereactome_dotplot", height = 500),
                                         downloadButton('download_gsereactome_dotplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsereactome_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Heatplot",
                                       fluidRow(
                                         box(numericInput("gsereactome_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("gsereactome_wide_heat", "Width", min = 1, max = 200, value = 80),
                                             width = 4),
                                         box(numericInput("gsereactome_high_heat", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsereactome_heatplot", height = 500),
                                         downloadButton('download_gsereactome_heatplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsereactome_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Cnetplot",
                                       fluidRow(
                                         box(numericInput("gsereactome_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                             width = 4),
                                         box(numericInput("gsereactome_wide_cnet", "Width", min = 1, max = 200, value = 30),
                                             width = 4),
                                         box(numericInput("gsereactome_high_cnet", "Height", min = 1, max = 200, value = 25),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         box(checkboxInput("gsereactome_circular_cnet",
                                                           "Circular",
                                                           value = TRUE),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsereactome_cnetplot", height = 500),
                                         downloadButton('download_gsereactome_cnetplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsereactome_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Emaplot",
                                       fluidRow(
                                         box(numericInput("gsereactome_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                             width = 4),
                                         box(numericInput("gsereactome_wide_ema", "Width", min = 1, max = 50, value = 11),
                                             width = 4),
                                         box(numericInput("gsereactome_high_ema", "Height", min = 1, max = 50, value = 10),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsereactome_emaplot", height = 500),
                                         downloadButton('download_gsereactome_emaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsereactome_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                              ),
                              tabPanel(title = "Gseaplot",
                                       fluidRow(box(uiOutput("gsereactome_term"), width = 12)),
                                       fluidRow(
                                         box(numericInput("gsereactome_wide_Gsea", "Width", min = 1, max = 50, value = 8),
                                             width = 4),
                                         box(numericInput("gsereactome_high_Gsea", "Height", min = 1, max = 50, value = 7),
                                             width = 4)
                                       ),
                                       fluidRow(
                                         plotOutput("gsereactome_Gseaplot", height = 500),
                                         downloadButton('download_gsereactome_Gseaplot', 'Save')
                                       ),
                                       shinyBS::bsTooltip("gsereactome_term", "The term that you want to show", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_wide_Gsea", "Width of the figure to export", "top", options = list(container = "body")),
                                       shinyBS::bsTooltip("gsereactome_high_Gsea", "Height of the figure to export", "top", options = list(container = "body"))
                              )
                       )
                       
                )
              )
      )
    )	
  )
)

## server ##--------
server <- function(input, output){
  options(shiny.maxRequestSize=200*1024^2)
  
  ## msqrob uioutput ##------------
  output$ID_colmun <- renderUI({
    selectizeInput("ID_colmun",
                   "Protein ID",
                   choices=c("none",colnames(file())),
                   selected = NULL)
  })
  
  output$Genename_colmun <- renderUI({
    selectizeInput("Genename_colmun",
                   "Genename",
                   choices=c("none",colnames(file())),
                   selected = NULL)
  })
  
  
  output$Control <- renderUI({
    selectizeInput("Control",
                   "Control Colmun",
                   choices=c(ctrl_ecols()),
                   selected = NULL)
  })
  
  output$downloadButton <- renderUI({
    downloadButton("downloadData", "Save", class = "downloadData")
  })
  
  ## annotation uioutput ##------
  output$genelist_tool_for_annotation <- renderUI({
    if(input$import_from_genelist_tool_for_annotation) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_annotation", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_annotation <- renderUI({
    if(!input$import_from_genelist_tool_for_annotation){
      textAreaInput(inputId = "text_for_annotation", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  ## msqrob reactive function ##------------------
  
  file <- reactive({
    file_location <- input$file
    if (is.null(file_location))
      return(NULL)
    read.table(file_location$datapath,header = T,sep = "\t")
  }) 
  
  ecols <- reactive({
    file <- file()
    get_ecols(file,input$intensity_colmun)
  })
  
  
  ctrl_ecols <- reactive({
    get_ctrl_ecols(ecols(),file())
  })
  
  
  corrected_file <-  reactive({
    get_corrected_file(input$Control,file(),ecols())
  })
  
  cond_vector <- reactive({
    condition_vector <- get_cond_vector(corrected_file(),ecols())
    # test_cond_vector <<- condition_vector
  })
  
  condition_vector <- reactive({
    condition_vector <- get_condition_vector(cond_vector())
    # test_condition_vector <<- condition_vector
  })
  
  condition_matched_frame <- reactive({
    condition_matched_frame <- get_condition_matched_frame(file = corrected_file(),ecols = ecols())
    test_condition_matched_frame <<- condition_matched_frame
  })
  
  # assign_result("condition_matched_frame",condition_matched_frame())
  
  condition <- reactive({
    get_condition_from_frame(condition_matched_frame())
  })
  
  
  ## reactive analysis ##-------
  pe <- reactive({
    pe <-  read2agg(corrected_file=corrected_file(),ecols=ecols(),
                    condition_matched_frame=condition_matched_frame(),ID_colmun=input$ID_colmun,
                    Genename_colmun=input$Genename_colmun,file_type=input$file_type)
    # pe_test <<- pe
  })
  
  pe_imp <- reactive({
    pe_imp <- my_imp(pe=pe(),filter_threshold=input$filter_threshold,cond_vector=cond_vector(),
                     condition_vector=condition_vector(),ID_colmun=input$ID_colmun,
                     imputation_type=input$imputation_type,file_type=input$file_type)
    # imp_test <<- pe_imp
  })
  
  pe_result <- reactive({
    pe_result <- my_test(pe=pe_imp(),condition_matched_frame=condition_matched_frame(),
                         file_type = input$file_type)
    # res_test <<-pe_result
  })
  
  Filtered_result <- reactive({
    Filtered_result <- Add_Filtered_result(object = pe_result(),
                                           LogFC = input$logFC,alpha = input$padj,
                                           addedname = "FilteredResult",
                                           condition_matched_frame = condition_matched_frame(),
                                           ID_colmun=input$ID_colmun)
    res_test <<-Filtered_result
  })
  
  ## GO uioutput##------
  output$genelist_tool_for_go <- renderUI({
    if(input$import_from_genelist_tool_for_go) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_go", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_go <- renderUI({
    if(!input$import_from_genelist_tool_for_go){
      textAreaInput(inputId = "text_for_go", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  ## KEGG uioutput##------
  output$genelist_tool_for_kegg <- renderUI({
    if(input$import_from_genelist_tool_for_kegg) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_kegg", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_kegg <- renderUI({
    if(!input$import_from_genelist_tool_for_kegg){
      textAreaInput(inputId = "text_for_kegg", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  ## reactome uioutput##------
  output$genelist_tool_for_reactome <- renderUI({
    if(input$import_from_genelist_tool_for_reactome) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_reactome", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_reactome <- renderUI({
    if(!input$import_from_genelist_tool_for_reactome){
      textAreaInput(inputId = "text_for_reactome", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  ## gsego uioutput##------
  output$genelist_tool_for_gsego <- renderUI({
    if(input$import_from_genelist_tool_for_gsego) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_gsego", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_gsego <- renderUI({
    if(!input$import_from_genelist_tool_for_gsego){
      textAreaInput(inputId = "text_for_gsego", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  ## gsekegg uioutput##------
  output$genelist_tool_for_gsekegg <- renderUI({
    if(input$import_from_genelist_tool_for_gsekegg) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_gsekegg", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_gsekegg <- renderUI({
    if(!input$import_from_genelist_tool_for_gsekegg){
      textAreaInput(inputId = "text_for_gsekegg", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  ## gsereactome uioutput##------
  output$genelist_tool_for_gsereactome <- renderUI({
    if(input$import_from_genelist_tool_for_gsereactome) {
      all_lists_of_genelist_tool <- get_all_lists()
      validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
      selectizeInput("genelist_for_gsereactome", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
    }
  })
  
  output$text_input_for_gsereactome <- renderUI({
    if(!input$import_from_genelist_tool_for_gsereactome){
      textAreaInput(inputId = "text_for_gsereactome", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
    }
  })
  
  
  ## observeEvent for msqrob  ---------------
  observeEvent(input$analyze,{
    
    
    
    output$cor_condition <- renderUI({
      selectizeInput("cor_condition",
                     "condition",
                     choices=c(condition()),
                     selected = NULL)
    })
    
    output$dis_condition <- renderUI({
      selectizeInput("dis_condition",
                     "condition",
                     choices=c(condition()),
                     selected = NULL)
    })
    
    output$select <- renderUI({
      res <- result_frame()
      cols <- grep("significant", colnames(res))
      names <- colnames(res)[cols]
      names <- gsub(".significant", "", names)
      selectizeInput("select",
                     "Select groups for comparing",
                     choices=names,
                     selected=names[1],
                     multiple = TRUE)
    })
    
    output$volcano_cond <- renderUI({
      selectizeInput("volcano_condition",
                     "condition",
                     choices=c(condition_matched_frame()$condition[2:nrow(condition_matched_frame())]),
                     selected = NULL)
    })
    
    ## for peptides file ##
    # if(input$file_type=="peptides"){
    #   pe <- readQFeatures(table = corrected_file(), fnames = 1, ecol = ecols(),name = "peptideRaw", sep = "\t")
    #   cat("1   ")
    #   
    #   
    #   col_vector <- reactive({
    #     get_col_vector(condition_matched_frame(),type = "peptide")
    #   })
    #   
    #   colData(pe)$condition <- as.factor(col_vector())
    #   cat("1   ")
    #   rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)##
    #   pe <- zeroIsNA(pe, "peptideRaw")
    #   #filter
    #   pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$Proteins%in% smallestUniqueGroups(rowData(pe[["peptideRaw"]])$Proteins), ]
    #   pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$Reverse != "+", ]
    #   pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$Potential.contaminant != "+", ]
    #   pe[["peptideRaw"]] <- pe[["peptideRaw"]][rowData(pe[["peptideRaw"]])$nNonZero >= 1, ]## 要改 ##
    #   
    #   pe_log <- pe
    #   pe_log <- logTransform_and_normalize(pe_log,method = "quantiles",type = "peptide")
    #   pe_log <- suppressWarnings(aggregateFeatures(pe_log,i = "peptideNorm",fcol = input$ID_colmun,na.rm = TRUE,name = "protein"))
    #   
    #   df_missval <- reactive({
    #     assay(pe_log[["protein"]])
    #   })
    #   
    #   cat("1   ")
    #   filtered <- my_filter_missval(assay(pe_log[["protein"]]),thr = input$filter_threshold,cond = cond_vector(),condition = condition_vector())
    #   colnames(filtered)[1] <- input$ID_colmun
    #   rowData(pe_log[["protein"]]) <- left_join(as.data.frame(rowData(pe_log[["protein"]])),
    #                                             filtered,by =input$ID_colmun)
    #   pe_log[["protein"]] <- pe_log[["protein"]][rowData(pe_log[["protein"]])$filter== "yes", ]
    #   cat("1   ")
    #   
    #   pe_imp <- pe_log
    #   if (input$imputation_type=="RF"){
    #     assay_protein <- assay(pe_imp[["protein"]])
    #     assay_protein <- missForest(assay_protein)
    #     # assay(test[["proteinNorm"]]) <- reactive({assay_protein$ximp})
    #     assay(pe_imp[["protein"]]) <- assay_protein$ximp
    #   }else if(input$imputation_type=="No imputation"){
    #     assay_protein <- assay(pe_imp[["protein"]])
    #     assay(pe_imp[["protein"]]) <- assay_protein
    #   }else{
    #     assay_protein <- assay(pe_imp[["protein"]])
    #     assay(pe_imp[["protein"]]) <- impute_matrix(assay_protein, method = input$imputation_type)
    #   }
    #   cat("1   ")
    #   
    #   pe_result <- pe_imp
    #   pe_result <- suppressWarnings(msqrob(object = pe_result,i = "protein",formula = ~condition))
    #   
    #   L <- my_makeContrast(condition_matched_frame())
    #   pe_result <- hypothesisTest(pe_result, i = "protein", L)
    #   
    #   test_pe <<- pe_result
    #   cat("1   ")
    #   
    #   assign_result(name = "primer",object = pe_result)
    #   
    # }
    # else if(input$file_type=="protein groups"){
    #   
    #   data_unique <- reactive({make_unique(corrected_file(), input$Genename_colmun, input$ID_colmun, delim = ";")})## 改 ##
    #   
    #   pe <-readQFeatures(table = data_unique(),fnames = 1,ecol = ecols(),name = "proteinRaw", sep="\t")
    #   
    #   cat("2   ")
    #   ## 要改 ##
    #   condition_matched_frame <- reactive({
    #     get_condition_matched_frame(file = data_unique(),ecols = ecols())
    #   })
    #   
    #   
    #   
    #   assign_result("condition_matched_frame",condition_matched_frame())
    #   
    #   col_vector <- reactive({
    #     get_col_vector(condition_matched_frame(),type = "protein")
    #   })
    #   
    #   colData(pe)$condition <- as.factor(col_vector())
    #   
    #   cat("2   ")
    #   pe <- zeroIsNA(pe, "proteinRaw")
    #   
    #   filtered <- reactive({
    #     my_filter_missval(assay(pe[["proteinRaw"]]),thr = input$filter_threshold,cond = cond_vector(),condition = condition_vector())
    #   })
    #   
    #   # test_filtered <<- filtered()
    #   
    #   cat("2   ")
    #   
    #   rowData(pe[["proteinRaw"]]) <- left_join(as.data.frame(rowData(pe[["proteinRaw"]])),
    #                                            filtered(),by ="Protein.IDs")
    #   try(pe[["proteinRaw"]] <- pe[["proteinRaw"]][rowData(pe[["proteinRaw"]])$Reverse != "+", ], silent = T)
    #   try(pe[["proteinRaw"]] <- pe[["proteinRaw"]][rowData(pe[["proteinRaw"]])$Potential.contaminant != "+", ], silent = T)
    #   pe[["proteinRaw"]] <- pe[["proteinRaw"]][rowData(pe[["proteinRaw"]])$filter== "yes", ]
    #   
    #   pe_log <- pe
    #   pe_log <- logTransform_and_normalize(pe_log,method = "quantiles",type = "protein")
    #   log_test <<- pe_log
    #   df_missval <- reactive({
    #     assay(pe_log[["protein"]])
    #   })
    #   
    #   ## imputation ##
    #   if (input$imputation_type=="RF"){
    #     assay_protein <- missForest(df_missval())
    #     # assay(test[["protein"]]) <- reactive({assay_protein$ximp})
    #     assay(pe_log[["protein"]]) <- assay_protein$ximp
    #   }else if(input$imputation_type=="No imputation"){
    #     assay(pe_log[["protein"]]) <- df_missval()
    #   }else{
    #     assay(pe_log[["protein"]]) <- impute_matrix(df_missval(), method = input$imputation_type)
    #   }
    #   
    #  
    #   pe_result <- pe_log
    #   pe_result <- suppressWarnings(msqrob(object = pe_log,i = "protein",formula = ~condition))
    #   L <- my_makeContrast(condition_matched_frame())
    #   
    #   
    #   pe_result <- hypothesisTest(pe_result, i = "protein", L)
    #   test_pro <<- pe_result
    #   # pe <- reactive({Add_Filtered_result(pe,type = "protein groups",LogFC = input$logFC,alpha = input$padj,addedname = "FilteredResult")})
    #   cat("2   ")
    #   assign_result(name = "primer",object = pe_result)
    #   
    # }
    
    ## msqrob reactive function ##------------
    msqrob_result <- reactive({
      # Add_Filtered_result(object = pe_result(),
      #                     LogFC = input$logFC,alpha = input$padj,
      #                     addedname = "FilteredResult",
      #                     condition_matched_frame = condition_matched_frame(),
      #                     ID_colmun=input$ID_colmun)
      msqrob_result <-Filtered_result()
    })
    
    result_frame <- reactive({
      result_from_object(object=msqrob_result(),condition_matched_frame())
    })
    # assign_result(name = "frame",object = result_frame())
    
    # res <- reactive({
    #   get("pe",envir = get_msqrob_result_Env())
    # })
    
    ## input function ##----------
    cat("3 ")
    MDS_input <- reactive({
      plot_MDS(msqrob_result(),type=input$MDStype,top=input$MDStop)
    })
    
    cat("4 ")
    cor_input <- reactive({
      plot_cor(object = msqrob_result(),condition_matched_frame = condition_matched_frame(),condition = input$cor_condition)
    })
    
    cat("5 ")
    # dist_input <- reactive({
    #   plot_dist(object = msqrob_result(),condition_matched_frame = condition_matched_frame(),condition = input$dis_condition)
    # })
    
    cat("6 ")
    CVs_input <- reactive({
      plot_cvs(msqrob_result())
    })
    
    cat("7 ")
    volcano_input <- reactive({
      plot_volcano(row_data = result_frame(),condition_matched_frame = condition_matched_frame(),
                   condition = input$volcano_condition,plot = TRUE)
    })
    
    cat("8 ")
    pca_input <- reactive({
      plot_pca(object = msqrob_result(),n=input$PCAtop,
               condition_matched_frame = condition_matched_frame(),plot = TRUE)
    })
    
    cat("9 ")
    heatmap_input <- reactive({
      plot_heatmap(object = msqrob_result(),row_data = result_frame(),
                   type = input$heatmaptype,clustering_distance = input$heatmapcluster,
                   select=input$select,condition_matched_frame=condition_matched_frame(),
                   col_limit = input$col_limit,color = input$colorbar)
    })
    
    cat("10 ")
    table_input <- reactive({
      frame <- result_frame()
      table_select(frame=result_frame(),select=input$select)
    })
    
    # output function ##-----------
    cat("11 ")
    output$table <- DT::renderDataTable({
      table_input()
    }, options = list(pageLength = 25, scrollX = T),
    selection = list(selected = c(1)))
    
    cat("12 ")
    output$MDS <- renderPlot({
      MDS_input()
    })
    
    cat("13 ")
    output$cor <- renderPlot({
      cor_input()
    })
    
    cat("14 ")
    # output$dist <- renderPlot({
    #   dist_input()
    # })
    
    cat("15 ")
    output$CVs <- renderPlot({
      CVs_input()
    })
    
    cat("16 ")
    output$volcano <- renderPlot({
      volcano_input()
    })
    
    cat("17 ")
    output$pca <- renderPlot({
      pca_input()
    })
    
    cat("18 ")
    output$heatmap <- renderPlot({
      heatmap_input()
    })
    
    # stop("done")
    cat("19 ")
    
    
    ## download function ##--------
    output$downloadData <- downloadHandler(
      filename = function() { paste("data",".txt", sep = "") },
      content = function(file) {
        write.table(result_frame(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )
    
    output$downloadHeatmap <- downloadHandler(
      filename = function(){
        paste0("Heatmap",".pdf")
      },
      content = function(file) {
        pdf(file, width = input$Heatmap_Width, height = input$Heatmap_Height)
        print(heatmap_input())
        dev.off()
      }
    )
    
    output$downloadVolcano <- downloadHandler(
      filename = function(){
        paste0("Volcano",".pdf")
      },
      content = function(file) {
        pdf(file, width = input$Volcano_Width, height = input$Volcano_Height)
        print(volcano_input())
        dev.off()
      }
    )
    
    output$downloadMDS <- downloadHandler(
      filename = function(){
        paste0("MDS",".pdf")
      },
      content = function(file) {
        pdf(file, width = input$MDS_Width, height = input$MDS_Height)
        print(MDS_input())
        dev.off()
      }
    )
    
    output$downloadPCA <- downloadHandler(
      filename = function(){
        paste0("PCA",".pdf")
      },
      content = function(file) {
        pdf(file, width = input$PCA_width, height = input$PCA_height)
        print(pca_input())
        dev.off()
      }
    )
    
    output$download_Pearson_correlation <- downloadHandler(
      filename = function(){
        paste0("Pearson_correlation",".pdf")
      },
      content = function(file) {
        pdf(file, width = input$Pearson_Width, height = input$Pearson_Height)
        print(cor_input())
        dev.off()
      }
    )
    
    output$download_Sample_CVs <- downloadHandler(
      filename = function(){
        paste0("Sample_CVs",".pdf")
      },
      content = function(file) {
        pdf(file, width = input$Sample_CVs_Width, height = input$Sample_CVs_Height)
        print(CVs_input())
        dev.off()
      }
    )
    
    
  })
  
  
  
 
  ## observeEvent for GO ##----------------
  
  # output$genelist_for_annotation <- renderUI({
  #   if(input$import_from_genelist_annotation) {
  #     all_lists_of_genelist_tool <- get_all_lists()
  #     validate(need(all_lists_of_genelist_tool != "", message = "There are no gene lists that can be obtained from gene list tool. Please do DEP-LFQ or DEG-RNAseq analysis or import gene list from gene list tool options"))
  #     selectizeInput("genelist_for_annotation", "Choose gene list", choices = get_all_lists(), selected = NULL, width = '100%')
  #   }
  # })
  
  
  # output$input_text <- renderUI({
  #   if(!input$import_from_genelist_annotation){
  #     textAreaInput(inputId = "input_text", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
  #   }
  # })
  observeEvent(
    input$analyze_for_go,{
      ### Interactive UI functions ### 
      output$downloadTable_go <- renderUI({
        selectizeInput("dataset_for_go",
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })
      
      output$downloadButton_go <- renderUI({
        downloadButton('downloadgo', 'Save table', class = "downloadgo")
      })
      
      if(!input$import_from_genelist_tool_for_go) {
        genelist <- reactive({ strsplit(input$text_for_go,'\n')[[1]]})
        gene_name <-reactive({unique(unlist(strsplit(genelist(),";")[]))})
      }else{
        genelist <- reactive({ get(input$genelist_for_go,envir = get_all_genelist_Env()) })
        # genelist <- reactive({ get(input$genelist_for_go,envir = get_msqrob_genelist_Env()) })
        gene_name <- reactive({unique(unlist(strsplit(genelist()$symbol,";")[]))})
      }
      
      # genelist <- reactive({ strsplit(input$text_input_for_go,'\n')[[1]] })
      # gene_name <-reactive(unlist(strsplit(genelist(),";")))
      
      gene_df <- reactive({
        gene_df <- data.frame(name=gene_name())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        if(!input$df_with_lg2fc){
          colnames(gene_df) = "name"
          gene_df <- as.data.frame(gene_df)
          test <- gene_df
        } else {
          colnames(gene_df) = c("name","fc")
          gene_df <- as.data.frame(gene_df)
          gene_df$fc = as.numeric(as.character(gene_df$fc))
          gene_df$name = as.character(gene_df$name)
          test <- gene_df
        }
        gene_df
      })
      # }
      cat("1   ")
      
      
      # if(!input$import_from_for_go) {
      gene_df_check <- try(gene_df())
      if(!class(gene_df_check) == "try-error") {
        gene_df_check$name = rm_digit_end(gene_df_check$name)
        if(ncol(gene_df_check) == 1) {
          gene_df <- data.frame(name = gene_df_check[!duplicated(gene_df_check$name),], stringsAsFactors = F)
        }
        if(ncol(gene_df_check) == 2) {
          gene_df <- gene_df_check[!duplicated(gene_df_check$name),]
        }
      }
      # }
      cat("1   ")
      # if(input$import_from_for_go) {
      #   gene_df$name = rm_digit_end(gene_df$name)
      #   if(ncol(gene_df) == 1) {
      #     gene_df <- data.frame(name = gene_df[!duplicated(gene_df$name),], stringsAsFactors = F)
      #   }
      #   if(ncol(gene_df) == 2) {
      #     gene_df <- gene_df[!duplicated(gene_df$name),]
      #   }
      # }
      
      pkg_for_go <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_go]  
      # reat <<- reactive({
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Please wait ...", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      reat <- try(goAnalysis(df = gene_df, df_with_lg2fc = input$df_with_lg2fc, organism = input$organism_for_go, species_df = annoSpecies_df), silent = TRUE)
      # })
      cat("1   ")
      # test1 <<- reat()
      res <- reactive({
        try(giveGO_res_and_table(reat = reat, ont = input$go_ont, pCutoff = input$go_p, p.adj.cutoff = input$go_padj, q.cutoff = input$go_qvalue, simplify = input$go_simplify), silent = TRUE)
      })
      
      test2 <- res()
      
      output$go_Table <- DT::renderDataTable({
        # if(!input$import_from_for_go) {
        #   shiny::validate(
        #     need(!class(gene_df_check) == "try-error", message = "No genes meet your requirements, and can not do the GO analysis")
        #   )      
        # } 
        
        # if(input$import_from_for_go) {
        #   shiny::validate(
        #     need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the GO analysis")
        #   )      
        # }    
        
        shiny::validate(need(require(pkg_for_go, character.only = TRUE), message = paste0("The package ", pkg_for_go, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_go, "')")))
        
        validate( need(reat, "No genes can be mapped to ENTREZID, and can not do the GO analysis"))   
        DT::datatable(res()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
      })
      
      # output$table <- DT::renderDataTable({
      #   table()
      # }, options = list(pageLength = 25, scrollX = T),
      # selection = list(selected = c(1)))
      
      
      output$significantBox_for_go <- renderInfoBox({
        # if(!input$import_from_for_go) {
        #   shiny::validate(
        #     need(!class(gene_df_check) == "try-error", message = "No genes meet your requirements, and can not do the GO analysis")
        #   )      
        # } 
        
        # if(input$import_from_for_go) {
        #   shiny::validate(
        #     need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the GO analysis")
        #   )      
        # }
        
        shiny::validate(need(require(pkg_for_go, character.only = TRUE), message = paste0("The package ", pkg_for_go, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_go, "')")))
        
        validate( need(reat, "No genes can be mapped to ENTREZID, and can not do the GO analysis"))
        
        
        num_total_go <- res()$all_table %>%
          nrow()
        num_signif_go <- res()$sig_table %>%
          nrow()
        frac_go <- num_signif_go / num_total_go
        cat("1   ")
        if(frac_go == 0) {
          info_box_go <- infoBox("Significant terms",
                                 paste0(num_signif_go,
                                        " out of ",
                                        num_total_go),
                                 "No terms enriched",
                                 icon = icon("thumbs-down", lib = "glyphicon"),
                                 color = "red",
                                 width = 4)
        }
        if(!frac_go == 0) {
          info_box_go <-     infoBox("Significant terms",
                                     paste0(num_signif_go,
                                            " out of ",
                                            num_total_go),
                                     paste0(signif(frac_go * 100, digits = 3),
                                            "% of terms enriched"),
                                     icon = icon("thumbs-up", lib = "glyphicon"),
                                     color = "green",
                                     width = 4)
        }
        info_box_go
      })
      
      go_barplot_input <- reactive({
        # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_barplot(res = res(), ShowCategory = input$go_ShowCategory_bar, color = input$go_color, ont = input$go_ont, Split = input$go_bar_if_Split_for_ont_ALL)
      })
      
      go_dotplot_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_dotplot(res = res(), ShowCategory = input$go_ShowCategory_dot, color = input$go_color, ont = input$go_ont, Split = input$go_dot_if_Split_for_ont_ALL)
      })
      
      # go_dotplot_opt_input <- reactive({
      #   validate(
      #     need(input$go_ont != "ALL", "Please go to the panel : Dot plot ")
      #   )
      #   progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      #   my_dotplot_opt(res = res(), color = input$go_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$go_ShowCategory_dot_opt)
      # })
      
      go_heatplot_input <- reactive({
        # progress <- shiny::Progress$new()
        #   progress$set(message = "Plotting", value = 0.66)
        #   # Close the progress when this reactive exits (even if there's an error)
        #   on.exit(progress$close())
        withProgress(message = 'Plotting', value = 0.66, {
          my_heatplot(res = res(), ShowCategory = input$go_ShowCategory_heat, df_with_lg2fc = input$df_with_lg2fc, ont = input$go_ont)        
        })
      })
      
      go_cnetplot_input <- reactive({
        # progress <- shiny::Progress$new()
        #   progress$set(message = "Plotting", value = 0.66)
        #   # Close the progress when this reactive exits (even if there's an error)
        #   on.exit(progress$close())
        withProgress(message = 'Plotting', value = 0.66, {
          my_cnetplot(res = res(), ShowCategory = input$go_ShowCategory_cnet, circular = input$go_circular_cnet, colorEdge = TRUE, df_with_lg2fc = input$df_with_lg2fc, ont = input$go_ont)        
        })
      })
      
      go_emaplot_input <- reactive({
        # progress <- shiny::Progress$new()
        #   progress$set(message = "Plotting", value = 0.66)
        #   # Close the progress when this reactive exits (even if there's an error)
        #   on.exit(progress$close())
        withProgress(message = 'Plotting', value = 0.66, {
          my_emaplot(res = res(), ShowCategory = input$go_ShowCategory_ema, color = input$go_color, layout = "kk", ont = input$go_ont)        
        })
      })
      
      go_goplot_input <- reactive({
        validate(
          need(input$go_ont != "ALL", "Ontology ALL: can not plot for goplot")
        )
        # progress <- shiny::Progress$new()
        #   progress$set(message = "Plotting", value = 0.66)
        #   # Close the progress when this reactive exits (even if there's an error)
        #   on.exit(progress$close())
        withProgress(message = 'Plotting', value = 0.66, {
          my_goplot(res = res(), ShowCategory = input$go_ShowCategory_go, color = input$go_color, ont = input$go_ont, Layout = "kk",circular = input$go_circular_go)        
        })
      })
      
      # go_GOgraphplot_input <- reactive({
      #   validate(
      #   need(input$go_ont != "ALL", "Ontology ALL: can not plot for plotGOgraph")
      # )
      #   # progress <- shiny::Progress$new()
      #   #   progress$set(message = "Plotting", value = 0.66)
      #   #   # Close the progress when this reactive exits (even if there's an error)
      #   #   on.exit(progress$close())
      #   # withProgress(message = 'Plotting', value = 0.66, {
      #   my_plotGOgraph(res = res(), firstSigNodes = input$go_ShowCategory_GOgraph, ont = input$go_ont)        
      #     # })
      #   })
      
      output$go_barplot <- renderPlot({
        go_barplot_input()
      })
      
      output$go_dotplot <- renderPlot({
        go_dotplot_input()
      })
      
      # output$go_dotplot_opt <- renderPlot({
      #   go_dotplot_opt_input()
      # })
      
      output$go_heatplot <- renderPlot({
        go_heatplot_input()
      })
      
      output$go_cnetplot <- renderPlot({
        go_cnetplot_input()
      })
      
      output$go_emaplot <- renderPlot({
        go_emaplot_input()
      })
      
      output$go_goplot <- renderPlot({
        go_goplot_input()
      })
      
      #  output$go_GOgraphplot <- renderPlot({
      #   go_GOgraphplot_input()
      # })
      
      
      ### Download objects and functions ### ------------------------------------
      datasetInput_for_go <- reactive({
        table_for_go = res()
        if(input$go_ont == "ALL") {
          if(input$go_bar_if_Split_for_ont_ALL | input$go_dot_if_Split_for_ont_ALL){
            switch(input$dataset_for_go,
                   "full_results" = res()$all_table,
                   "significant_results" = res()$sig_table)
            
          } else {
            switch(input$dataset_for_go,
                   "full_results" = res()$all_table %>% dplyr::arrange(p.adjust),
                   "significant_results" = res()$sig_table %>% dplyr::arrange(p.adjust))
            
          }
          
        } else {
          switch(input$dataset_for_go,
                 "full_results" = res()$all_table,
                 "significant_results" = res()$sig_table)
        }
        # switch(input$dataset_for_go,
        #        "full_results" = res()$all_table,
        #        "significant_results" = res()$sig_table)
      })
      
      
      output$downloadgo <- downloadHandler(
        filename = function() { paste(input$dataset_for_go, ".txt", sep = "") },
        content = function(file) {
          write.table(datasetInput_for_go(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t") }
      )
      
      output$download_go_barplot <- downloadHandler(
        filename = 'barplot.pdf',
        content = function(file) {
          pdf(file, width = input$go_wide_bar, height = input$go_high_bar)
          print(go_barplot_input())
          dev.off()
        }
      )
      
      output$download_go_dotplot <- downloadHandler(
        filename = 'dotplot.pdf',
        content = function(file) {
          pdf(file, width = input$go_wide_dot, height = input$go_high_dot)
          print(go_dotplot_input())
          dev.off()
        }
      )
      
      # output$download_go_dotplot_opt <- downloadHandler(
      #   filename = 'dotplot_opt.pdf',
      #   content = function(file) {
      #     pdf(file, width = input$go_wide_dot_opt, height = input$go_high_dot_opt)
      #     print(go_dotplot_opt_input())
      #     dev.off()
      #   }
      # )
      
      output$download_go_heatplot <- downloadHandler(
        filename = 'heatplot.pdf',
        content = function(file) {
          pdf(file, width = input$go_wide_heat, height = input$go_high_heat)
          print(go_heatplot_input())
          dev.off()
        }
      )
      
      output$download_go_cnetplot <- downloadHandler(
        filename = 'cnetplot.pdf',
        content = function(file) {
          pdf(file, width = input$go_wide_cnet, height = input$go_high_cnet)
          print(go_cnetplot_input())
          dev.off()
        }
      )
      
      output$download_go_emaplot <- downloadHandler(
        filename = 'emaplot.pdf',
        content = function(file) {
          pdf(file, width = input$go_wide_ema, height = input$go_high_ema)
          print(go_emaplot_input())
          dev.off()
        }
      )
      
      output$download_go_goplot <- downloadHandler(
        filename = 'goplot.pdf',
        content = function(file) {
          pdf(file, width = input$go_wide_go, height = input$go_high_go)
          print(go_goplot_input())
          dev.off()
        }
      )
      
      # output$download_go_GOgraphplot <- downloadHandler(
      #   filename = 'GOgraphplot.pdf',
      #   content = function(file) {
      #     pdf(file, width = input$go_wide_GOgraph, height = input$go_high_GOgraph)
      #     print(go_GOgraphplot_input())
      #     dev.off()
      #   }
      # )
      
      observeEvent(input$plot_for_GOgraphplot, {
        go_GOgraphplot_input <- reactive({
          validate(
            need(input$go_ont != "ALL", "Ontology ALL: can not plot for plotGOgraph")
          )
          # progress <- shiny::Progress$new()
          #   progress$set(message = "Plotting", value = 0.66)
          #   # Close the progress when this reactive exits (even if there's an error)
          #   on.exit(progress$close())
          withProgress(message = 'Plotting', value = 0.66, {
            my_plotGOgraph(res = res(), firstSigNodes = input$go_ShowCategory_GOgraph, ont = input$go_ont)        
          })
        })
        
        output$go_GOgraphplot <- renderPlot({
          go_GOgraphplot_input()
        })
        
        output$download_go_GOgraphplot <- downloadHandler(
          filename = 'GOgraphplot.pdf',
          content = function(file) {
            pdf(file, width = input$go_wide_GOgraph, height = input$go_high_GOgraph)
            print(go_GOgraphplot_input())
            dev.off()
          }
        )
      })
    }
  )
  
  ## observeEvent for annotation ##---------
  observeEvent(
    input$analyze_for_annotation,{
      withProgress(message = 'Please wait', value = 0.66, {
        output$downloadannotation_for_output <- renderUI({
          #align left;width: 10em; vertical-align: center
          downloadButton("downloadannotation", "Save table", class = "downloadannotation")#btn btn-success
        })
        
        # genelist <- strsplit(input$text_input_for_annotation,'\n')[[1]]
        # gene_name <-reactive({unlist(strsplit(genelist,";")[])})
        # gene_df <<- data.frame(name=gene_name())
        cat("yeah")
        
        if(!input$import_from_genelist_tool_for_annotation) {
          genelist <- strsplit(input$text_for_annotation,'\n')[[1]]
          gene_name <-unique(unlist(strsplit(genelist,";")[]))
          gene_df <- data.frame(name=gene_name)
        }else{
          genelist <- get(input$genelist_for_annotation,envir = get_all_genelist_Env())
          gene_name <- unique(unlist(strsplit(genelist$symbol,";")[]))
          gene_df <- data.frame(name=gene_name)
        }
        
        gene_df$name = rm_digit_end(gene_df$name)
        gene_df = data.frame(name = gene_df[!duplicated(gene_df$name),], stringsAsFactors = F)
        
        cat("yeah")
        
        if(!nrow(gene_df) == 0) {
          gene_df <- gene_df %>% dplyr::filter(!is.na(name))
        }
        
        if(nrow(gene_df) == 0) {
          output$annotation_Table <- DT::renderDataTable({
            shiny::validate(need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the annotation analysis"))
            DT::datatable(data.frame(a = "I am for check"), filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
          })
        }
        shiny::validate(need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the annotation analysis"))
        
        gene_df$name <-  sapply(gene_df$name,the1stname)
        
        # pkg_for_annotation <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_annotation]
        pkg_for_annotation <- "org.Hs.eg.db"
        if(!require(pkg_for_annotation, character.only = TRUE)) {
          output$annotation_Table <- DT::renderDataTable({
            shiny::validate(need(require(pkg_for_annotation, character.only = TRUE), message = paste0("The package ", pkg_for_annotation, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_annotation, "')")))
            DT::datatable(data.frame(a = "I am for check"), filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
          })
        }
        
        shiny::validate(
          need(
            require(pkg_for_annotation, character.only = TRUE),
            paste0("The package ", pkg_for_annotation, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_annotation, "')")
          )
        )    
        
        cat("yeah")
        # if(input$organism_for_annotation == "Human"){
        #   library(org.Hs.eg.db)
        #   orgDB <- org.Hs.eg.db
        #   # kegg_organism <- "hsa"
        # }else if(input$organism_for_annotation == "Mouse"){
        #   library(org.Mm.eg.db)
        #   orgDB <- org.Mm.eg.db
        #   # kegg_organism <- "mmu"
        # }
        # ids <- gene_df$ENTREZID <-  mapIds(x = orgDB,
        #                              keys = gene_df$name,
        #                              keytype = "SYMBOL",
        #                              column = "ENTREZID")
        check_ids <- my_to_entrezid(orgDB = get(pkg_for_annotation), gene = as.character(gene_df$name))$id
        # shiny::validate(need(class(check_ids) != "logical", message = "No genes can be mapped to ENTREZID, and can not do the  annotation analysis")) 
        if(class(check_ids) == "logical") {
          output$annotation_Table <- DT::renderDataTable({
            shiny::validate(need(class(check_ids) != "logical", message = "No genes can be mapped to ENTREZID, and can not do the  annotation analysis"))
            DT::datatable(data.frame(a = "I am for check"), filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
          })
        }  
        shiny::validate(need(class(check_ids) != "logical", message = "No genes can be mapped to ENTREZID, and can not do the  annotation analysis"))                             
        ids <- gene_df$ENTREZID <- check_ids
        
        gene_df$ENTREZID <- as.character(gene_df$ENTREZID)
        
        ids <- na.omit(unlist(ids))
        
        ids <- as.vector(ids)
        # if(input$organism_for_annotation == "human"){
        #   library(org.Hs.eg.db)
        #   orgDB <- org.Hs.eg.db
        #   kegg_organism <- "hsa"
        # }else if(input$organism_for_annotation == "mouse"){
        #   library(org.Mm.eg.db)
        #   orgDB <- org.Mm.eg.db
        #   kegg_organism <- "mmu"
        # }
        bgo <- Geneannotate(ids,genedb=get(pkg_for_annotation))
        IDS <- unique(bgo$ENTREZID)
        
        cat("yeah")
        
        library(parallel)
        cpus = detectCores(logical = F)
        cl <- makeCluster(cpus)
        annots <- parSapply(cl,IDS,FUN  = mergego,bgo=bgo)
        stopCluster(cl)
        
        rownames(annots)<-c("GO.CC.name","GO.BP.name","GO.MF.name", "PFAM", "Gene.Description", "KEGG.name", "reactome.name")
        annots2 <- annots %>% t() %>%as.data.frame() %>%mutate(ENTREZID=colnames(annots))
        mergepg <- merge(gene_df,annots2,by.x="ENTREZID",by.y="ENTREZID")
        mergepg = mergepg[ , c(2, 7, 1, 3:6, 8:9)]
        
        
        output$annotation_Table <- DT::renderDataTable(
          mergepg, filter = 'top', options = list( autoWidth = F,scrollX = TRUE
          )
        )
        
        output$downloadannotation <- downloadHandler(
          filename = function() { paste("gene annotation", ".txt", sep = "") },
          content = function(file) {
            write.table(mergepg,
                        file,
                        col.names = TRUE,
                        row.names = FALSE,
                        sep ="\t") }
        ) 
      })
    }
  )
  
  ## observeEvent for KEGG ##--------
  observeEvent(
    input$analyze_for_kegg,{
      ### Interactive UI functions ### ------------------------------------------
      output$downloadTable_kegg <- renderUI({
        selectizeInput("dataset_for_kegg",
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })
      
      output$downloadButton_kegg <- renderUI({
        downloadButton('downloadkegg', 'Save table', class = "downloadkegg")
      })
      
      if(!input$import_from_genelist_tool_for_kegg) {
        genelist_kegg <- reactive({ strsplit(input$text_for_kegg,'\n')[[1]]})
        gene_name_kegg <-reactive({unique(unlist(strsplit(genelist_kegg(),";")[]))})
      }else{
        genelist_kegg <- reactive({ get(input$genelist_for_kegg,envir = get_all_genelist_Env()) })
        gene_name_kegg <- reactive({unique(unlist(strsplit(genelist_kegg()$symbol,";")[]))})
      }
      
      # testgene <<- gene_name_kegg()
      # genelist_kegg <- reactive({ strsplit(input$text_input_for_kegg,'\n')[[1]] })
      # gene_name_kegg <-reactive(unlist(strsplit(genelist_kegg(),";")))
      
      gene_df_kegg <- reactive({
        gene_df <- data.frame(name=toupper(gene_name_kegg()))
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        if(!input$df_with_lg2fc_for_kegg){
          colnames(gene_df) = "name"
          gene_df = as.data.frame(gene_df)
        } else {
          colnames(gene_df) = c("name","fc")
          gene_df = as.data.frame(gene_df)
          gene_df$fc = as.numeric(as.character(gene_df$fc))
          gene_df$name = as.character(gene_df$name)
          # test <<- gene_df
        }
        gene_df    
      })    	
      # }
      # testgene <<- gene_df_kegg()
      cat("2   ")
      
      
      # if(!input$import_from_for_kegg) {
      gene_df_kegg_check <- try(gene_df_kegg())
      if(!class(gene_df_kegg_check) == "try-error") {
        gene_df_kegg_check$name = rm_digit_end(gene_df_kegg_check$name)
        if(ncol(gene_df_kegg_check) == 1) {
          gene_df_kegg <- data.frame(name = gene_df_kegg_check[!duplicated(gene_df_kegg_check$name),], stringsAsFactors = F)
        }
        if(ncol(gene_df_kegg_check) == 2) {
          gene_df_kegg <- gene_df_kegg_check[!duplicated(gene_df_kegg_check$name),]
        }
        
      }
      
      # testgene <<- gene_df_kegg
      # }
      cat("2   ")
      # if(input$import_from_for_kegg) {
      #   gene_df_kegg$name = rm_digit_end(gene_df_kegg$name)
      #   if(ncol(gene_df_kegg) == 1) {
      #     gene_df_kegg <- data.frame(name = gene_df_kegg[!duplicated(gene_df_kegg$name),], stringsAsFactors = F)
      #   }
      #   if(ncol(gene_df_kegg) == 2) {
      #     gene_df_kegg <- gene_df_kegg[!duplicated(gene_df_kegg$name),]
      #   }
      # }
      
      pkg_for_kegg <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_kegg]  
      
      # reat_kegg <<- reactive({
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Please wait ...", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      reat_kegg <- try(keggAnalysis(df = gene_df_kegg, organism = input$organism_for_kegg, df_with_lg2fc = input$df_with_lg2fc_for_kegg, species_df = annoSpecies_df), silent = TRUE)      
      # })
      test1<-reat_kegg
      res_kegg <- reactive({
        # if(input$import_from_for_kegg) {
        #   shiny::validate(
        #   need(nrow(gene_df_kegg) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
        # )      
        # }
        #    # validate( need(reat_kegg, "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
        #    shiny::validate( need(class(reat_kegg) != "try-error", message = "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
        try(givekegg_reat_res_and_table(reat = reat_kegg, pCutoff = input$kegg_p, p.adj.cutoff = input$kegg_padj, q.cutoff = input$kegg_qvalue), silent = TRUE)
      })
      
      # test2 <<- res_kegg()
      
      cat("2   ")
      # output$kegg_Table <- DT::renderDataTable(
      #   res_kegg()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE
      #   )
      # )
      output$kegg_Table <- DT::renderDataTable({
        # if(!input$import_from_for_kegg) {
        #   shiny::validate(
        #     need(!class(gene_df_kegg_check) == "try-error", message = "No genes meet your requirements, and can not do the KEGG analysis")
        #   )      
        # }
        
        # if(input$import_from_for_kegg) {
        #   shiny::validate(
        #     need(nrow(gene_df_kegg) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
        #   )      
        # }
        
        shiny::validate(need(require(pkg_for_kegg, character.only = TRUE), message = paste0("The package ", pkg_for_kegg, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_kegg, "')")))
        
        validate( need(reat_kegg, "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
        # shiny::validate( need(class(reat_kegg) != "try-error", message = "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
        DT::datatable(res_kegg()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
      })
      
      output$significantBox_for_kegg <- renderInfoBox({
        # if(!input$import_from_for_kegg) {
        #   shiny::validate(
        #     need(!class(gene_df_kegg_check) == "try-error", message = "No genes meet your requirements, and can not do the KEGG analysis")
        #   )      
        # }
        
        # if(input$import_from_for_kegg) {
        #   shiny::validate(
        #     need(nrow(gene_df_kegg) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
        #   )      
        # }
        cat("2   ")
        shiny::validate(need(require(pkg_for_kegg, character.only = TRUE), message = paste0("The package ", pkg_for_kegg, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_kegg, "')")))
        
        validate( need(res_kegg, "No genes can be mapped to ENTREZID, and can not do the reactome analysis"))
        
        num_total_kegg <- res_kegg()$all_table %>%
          nrow()
        num_signif_kegg <- res_kegg()$sig_table %>%
          nrow()
        frac_kegg <- num_signif_kegg / num_total_kegg
        cat("2   ")
        if(frac_kegg == 0) {
          info_box_kegg <- infoBox("Significant terms",
                                   paste0(num_signif_kegg,
                                          " out of ",
                                          num_total_kegg),
                                   "No terms enriched",
                                   icon = icon("thumbs-down", lib = "glyphicon"),
                                   color = "red",
                                   width = 4)
        }
        if(!frac_kegg == 0) {
          info_box_kegg <-     infoBox("Significant terms",
                                       paste0(num_signif_kegg,
                                              " out of ",
                                              num_total_kegg),
                                       paste0(signif(frac_kegg * 100, digits = 3),
                                              "% of terms enriched"),
                                       icon = icon("thumbs-up", lib = "glyphicon"),
                                       color = "green",
                                       width = 4)
        }
        info_box_kegg
      })
      
      kegg_barplot_input <- reactive({
        # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        barplot(res_kegg()$sig_res, showCategory = input$kegg_ShowCategory_bar, color = input$kegg_color)
      })
      
      kegg_dotplot_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        dotplot(res_kegg()$sig_res, showCategory = input$kegg_ShowCategory_dot, color = input$kegg_color)
      })
      
      # kegg_dotplot_opt_input <- reactive({
      #   progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      #   my_dotplot_opt(res = res_kegg(), color = input$kegg_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$kegg_ShowCategory_dot_opt)
      # })
      
      kegg_heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          heatplot_for_react_kegg(res = res_kegg(), ShowCategory = input$kegg_ShowCategory_heat, df_with_lg2fc = input$df_with_lg2fc_for_kegg)        
        })
      })
      
      kegg_cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          cnetplot_for_react_kegg(res = res_kegg(), ShowCategory = input$kegg_ShowCategory_cnet, circular = input$kegg_circular_cnet, colorEdge = TRUE, df_with_lg2fc = input$df_with_lg2fc_for_kegg)        
        })
      })
      
      kegg_emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          emaplot_for_react_kegg(res = res_kegg(), ShowCategory = input$kegg_ShowCategory_ema, color = input$kegg_color, layout = "kk")        
        })
      })
      
      output$kegg_barplot <- renderPlot({
        kegg_barplot_input()
      })
      
      output$kegg_dotplot <- renderPlot({
        kegg_dotplot_input()
      })
      
      # output$kegg_dotplot_opt <- renderPlot({
      #   kegg_dotplot_opt_input()
      # })
      
      output$kegg_heatplot <- renderPlot({
        kegg_heatplot_input()
      })
      
      output$kegg_cnetplot <- renderPlot({
        kegg_cnetplot_input()
      })
      
      output$kegg_emaplot <- renderPlot({
        kegg_emaplot_input()
      })
      
      ### Download objects and functions ### ------------------------------------
      datasetInput_for_kegg <- reactive({
        table_for_kegg = res_kegg()
        switch(input$dataset_for_kegg,
               "full_results" = res_kegg()$all_table,
               "significant_results" = res_kegg()$sig_table)
      })
      
      output$downloadkegg <- downloadHandler(
        filename = function() { paste(input$dataset_for_kegg, ".txt", sep = "") },
        content = function(file) {
          write.table(datasetInput_for_kegg(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t") }
      )
      
      output$download_kegg_barplot <- downloadHandler(
        filename = 'barplot.pdf',
        content = function(file) {
          pdf(file, width = input$kegg_wide_bar, height = input$kegg_high_bar)
          print(kegg_barplot_input())
          dev.off()
        }
      )
      
      output$download_kegg_dotplot <- downloadHandler(
        filename = 'dotplot.pdf',
        content = function(file) {
          pdf(file, width = input$kegg_wide_dot, height = input$kegg_high_dot)
          print(kegg_dotplot_input())
          dev.off()
        }
      )
      
      # output$download_kegg_dotplot_opt <- downloadHandler(
      #   filename = 'dotplot_opt.pdf',
      #   content = function(file) {
      #     pdf(file, width = input$kegg_wide_dot_opt, height = input$kegg_high_dot_opt)
      #     print(kegg_dotplot_opt_input())
      #     dev.off()
      #   }
      # )
      
      output$download_kegg_heatplot <- downloadHandler(
        filename = 'heatplot.pdf',
        content = function(file) {
          pdf(file, width = input$kegg_wide_heat, height = input$kegg_high_heat)
          print(kegg_heatplot_input())
          dev.off()
        }
      )
      
      output$download_kegg_cnetplot <- downloadHandler(
        filename = 'cnetplot.pdf',
        content = function(file) {
          pdf(file, width = input$kegg_wide_cnet, height = input$kegg_high_cnet)
          print(kegg_cnetplot_input())
          dev.off()
        }
      )
      
      output$download_kegg_emaplot <- downloadHandler(
        filename = 'emaplot.pdf',
        content = function(file) {
          pdf(file, width = input$kegg_wide_ema, height = input$kegg_high_ema)
          print(kegg_emaplot_input())
          dev.off()
        }
      )
    }
    
    
  )
  
  ## observeEvent for reactome ##-------
  observeEvent(
      input$analyze_for_reactome,{
      ### Interactive UI functions ### ------------------------------------------
      output$downloadTable_reactome <- renderUI({
        selectizeInput("dataset_for_reactome",
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })
      
      output$downloadButton_reactome <- renderUI({
        downloadButton('downloadreactome', 'Save table', class = "downloadreactome")
      })
      
      if(!input$import_from_genelist_tool_for_reactome) {
        genelist_reactome <- reactive({ strsplit(input$text_for_reactome,'\n')[[1]]})
        gene_name_reactome <-reactive({unique(unlist(strsplit(genelist_reactome(),";")[]))})
      }else{
        genelist_reactome <- reactive({ get(input$genelist_for_reactome,envir = get_all_genelist_Env()) })
        gene_name_reactome <- reactive({unique(unlist(strsplit(genelist_reactome()$symbol,";")[]))})
      }
      
      # genelist_reactome <- reactive({ strsplit(input$text_input_for_reactome,'\n')[[1]] })
      # gene_name_reactome <-reactive(unlist(strsplit(genelist_reactome(),";")))  
      
      gene_df_reactome <- reactive({
        gene_df <- data.frame(name=toupper(gene_name_reactome()))
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        if(!input$df_with_lg2fc_for_reactome){
          colnames(gene_df) = "name"
          gene_df = as.data.frame(gene_df)
        } else {
          colnames(gene_df) = c("name","fc")
          gene_df = as.data.frame(gene_df)
          gene_df$fc = as.numeric(as.character(gene_df$fc))
          gene_df$name = as.character(gene_df$name)
          # test <<- gene_df
        }
        gene_df  
      })      
      # }
      cat("3   ")
      # if(input$import_from_for_reactome) {
      #   if(input$import_for_reactome == "UPregu for DEP-LFQ" | input$import_for_reactome == "DOWNregu for DEP-LFQ" | input$import_for_reactome == "UPDOWN for DEP-LFQ") {
      #     if(input$import_contrast_for_reactome == "Any significant") {
      #       if(input$import_for_reactome == "UPDOWN for DEP-LFQ") {
      #         gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(significant) %>% dplyr::select(name)
      #       } else {
      #         if(input$import_for_reactome == "UPregu for DEP-LFQ") {
      #           index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
      #           df <- get_results(dep()) %>% dplyr::filter(significant)              
      #           cols_diff = df %>% dplyr::select(index)
      #           cols_diff_reject = cols_diff > 0
      #           cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
      #           index1 = apply(cols_diff_reject, 1, all)
      #           gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
      #         } else {
      #           if(input$import_for_reactome == "DOWNregu for DEP-LFQ") {
      #             index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
      #             df <- get_results(dep()) %>% dplyr::filter(significant)
      #             cols_diff = df %>% dplyr::select(index)
      #             cols_diff_reject = cols_diff < 0
      #             cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
      #             index1 = apply(cols_diff_reject, 1, all)
      #             gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
      #           }
      #         }
      #       }
      #     } else {
      #       if(input$import_for_reactome == "UPDOWN for DEP-LFQ") {
      #         if(!input$df_with_lg2fc_for_reactome) {
      #           gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name)
      #         } else {
      #           gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "ratio", sep = "_"))
      #         }
      #       } else {
      #         if(input$import_for_reactome == "UPregu for DEP-LFQ") {
      #           if(!input$df_with_lg2fc_for_reactome) {
      #             gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) > 0) %>% dplyr::select(name)
      #           } else {
      #             gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "ratio", sep = "_"))
      #           }                  
      #         } else {
      #           if(input$import_for_reactome == "DOWNregu for DEP-LFQ") {
      #             if(!input$df_with_lg2fc_for_reactome) {
      #               gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) < 0) %>% dplyr::select(name)
      #             } else {
      #               gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "ratio", sep = "_"))
      #             }
      #           }
      #         }
      #       }
      #     }
      #   } else {
      #     if(input$import_for_reactome == "UPregu for DEG-RNAseq" | input$import_for_reactome == "DOWNregu for DEG-RNAseq" | input$import_for_reactome == "UPDOWN for DEG-RNAseq") {
      #       if(input$import_contrast_for_reactome == "Any significant") {
      #         if(input$import_for_reactome == "UPDOWN for DEG-RNAseq") {
      #           gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant) %>% dplyr::select(name)
      #         } else {
      #           if(input$import_for_reactome == "UPregu for DEG-RNAseq") {
      #             index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
      #             df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
      #             cols_diff = df %>% dplyr::select(index)
      #             cols_diff_reject = cols_diff > 0
      #             cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
      #             index1 = apply(cols_diff_reject, 1, all)
      #             gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
      #           } else {
      #             if(input$import_for_reactome == "DOWNregu for DEG-RNAseq") {
      #               index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
      #               df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
      #               cols_diff = df %>% dplyr::select(index)
      #               cols_diff_reject = cols_diff < 0
      #               cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
      #               index1 = apply(cols_diff_reject, 1, all)
      #               gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
      #             }
      #           }
      #         }
      #       } else {
      #         if(input$import_for_reactome == "UPDOWN for DEG-RNAseq") {
      #           if(!input$df_with_lg2fc_for_reactome) {
      #             gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name)
      #           } else {
      #             gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_"))
      #           }                  
      #         } else {
      #           if(input$import_for_reactome == "UPregu for DEG-RNAseq") {
      #             if(!input$df_with_lg2fc_for_reactome) {
      #               gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name)
      #             } else {
      #               gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_"))
      #             }                      
      #           } else {
      #             if(input$import_for_reactome == "DOWNregu for DEG-RNAseq") {
      #               if(!input$df_with_lg2fc_for_reactome) {
      #                 gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name)
      #               } else {
      #                 gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_"))
      #               }                          
      #             }
      #           }
      #         }
      #       }
      #     }
      #     
      #   }
      # }
      
      # if(!input$import_from_for_reactome) {
      gene_df_reactome_check <- try(gene_df_reactome())
      if(!class(gene_df_reactome_check) == "try-error") {
        gene_df_reactome_check$name = rm_digit_end(gene_df_reactome_check$name)
        if(ncol(gene_df_reactome_check) == 1) {
          gene_df_reactome <- data.frame(name = gene_df_reactome_check[!duplicated(gene_df_reactome_check$name),], stringsAsFactors = F)
        }
        if(ncol(gene_df_reactome_check) == 2) {
          gene_df_reactome <- gene_df_reactome_check[!duplicated(gene_df_reactome_check$name),]
        }
        
      }
      # }
      cat("3    ")
      # if(input$import_from_for_reactome) {
      #   gene_df_reactome$name = rm_digit_end(gene_df_reactome$name)
      #   if(ncol(gene_df_reactome) == 1) {
      #     gene_df_reactome <- data.frame(name = gene_df_reactome[!duplicated(gene_df_reactome$name),], stringsAsFactors = F)
      #   }
      #   if(ncol(gene_df_reactome) == 2) {
      #     gene_df_reactome <- gene_df_reactome[!duplicated(gene_df_reactome$name),]
      #   }
      # }
      
      pkg_for_reactome <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_reactome]      
      
      # reat_reactome <<- reactive({
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Please wait ...", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      reat_reactome <- try(reactAnalysis(df = gene_df_reactome, organism = input$organism_for_reactome, df_with_lg2fc = input$df_with_lg2fc_for_reactome, species_df = annoSpecies_df_for_reactome), silent = TRUE)
      # })
      
      # test1 <<- reat_reactome()
      res_reactome <- reactive({
        try(givekegg_reat_res_and_table(reat = reat_reactome, pCutoff = input$reactome_p, p.adj.cutoff = input$reactome_padj, q.cutoff = input$reactome_qvalue), silent = TRUE)
      })
      
      test2 <- res_reactome()
      cat("3   ")
      output$reactome_Table <- DT::renderDataTable({
        # if(!input$import_from_for_reactome) {
        #   shiny::validate(
        #     need(!class(gene_df_reactome_check) == "try-error", message = "No genes meet your requirements, and can not do the reactome analysis")
        #   )      
        # }
        
        # if(input$import_from_for_reactome) {
        #   shiny::validate(
        #     need(nrow(gene_df_reactome) != 0, message = "No genes meet your requirements, and can not do the reactome analysis")
        #   )      
        # }
        
        shiny::validate(need(require(pkg_for_reactome, character.only = TRUE), message = paste0("The package ", pkg_for_reactome, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_reactome, "')")))
        
        validate( need(res_reactome(), "No genes can be mapped to ENTREZID, and can not do the reactome analysis"))
        DT::datatable(res_reactome()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
      })
      
      output$significantBox_for_reactome <- renderInfoBox({
        # if(!input$import_from_for_reactome) {
        #   shiny::validate(
        #     need(!class(gene_df_reactome_check) == "try-error", message = "No genes meet your requirements, and can not do the reactome analysis")
        #   )      
        # }
        
        # if(input$import_from_for_reactome) {
        #   shiny::validate(
        #     need(nrow(gene_df_reactome) != 0, message = "No genes meet your requirements, and can not do the reactome analysis")
        #   )      
        # }
        
        shiny::validate(need(require(pkg_for_reactome, character.only = TRUE), message = paste0("The package ", pkg_for_reactome, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_reactome, "')")))
        
        validate( need(res_reactome, "No genes can be mapped to ENTREZID, and can not do the reactome analysis"))
        
        num_total_reactome <- res_reactome()$all_table %>%
          nrow()
        num_signif_reactome <- res_reactome()$sig_table %>%
          nrow()
        frac_reactome <- num_signif_reactome / num_total_reactome
        cat("3   ")
        if(frac_reactome == 0) {
          info_box_reactome <- infoBox("Significant terms",
                                       paste0(num_signif_reactome,
                                              " out of ",
                                              num_total_reactome),
                                       "No terms enriched",
                                       icon = icon("thumbs-down", lib = "glyphicon"),
                                       color = "red",
                                       width = 4)
        }
        if(!frac_reactome == 0) {
          info_box_reactome <-     infoBox("Significant terms",
                                           paste0(num_signif_reactome,
                                                  " out of ",
                                                  num_total_reactome),
                                           paste0(signif(frac_reactome * 100, digits = 3),
                                                  "% of terms enriched"),
                                           icon = icon("thumbs-up", lib = "glyphicon"),
                                           color = "green",
                                           width = 4)
        }
        info_box_reactome
      })
      
      reactome_barplot_input <- reactive({
        # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        barplot(res_reactome()$sig_res, showCategory = input$reactome_ShowCategory_bar, color = input$reactome_color)
      })
      
      reactome_dotplot_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        dotplot(res_reactome()$sig_res, showCategory = input$reactome_ShowCategory_dot, color = input$reactome_color)
      })
      
      # reactome_dotplot_opt_input <- reactive({
      #   progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      #   my_dotplot_opt(res = res_reactome(), color = input$reactome_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$reactome_ShowCategory_dot_opt)
      # })
      
      reactome_heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          heatplot_for_react_kegg(res = res_reactome(), ShowCategory = input$reactome_ShowCategory_heat, df_with_lg2fc = input$df_with_lg2fc_for_reactome)        
        })
      })
      
      reactome_cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          cnetplot_for_react_kegg(res = res_reactome(), ShowCategory = input$reactome_ShowCategory_cnet, circular = input$reactome_circular_cnet, colorEdge = TRUE, df_with_lg2fc = input$df_with_lg2fc_for_reactome)        
        })
      })
      
      reactome_emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          emaplot_for_react_kegg(res = res_reactome(), ShowCategory = input$reactome_ShowCategory_ema, color = input$reactome_color, layout = "kk")        
        })
      })
      
      output$reactome_barplot <- renderPlot({
        reactome_barplot_input()
      })
      
      output$reactome_dotplot <- renderPlot({
        reactome_dotplot_input()
      })
      
      output$reactome_dotplot_opt <- renderPlot({
        reactome_dotplot_opt_input()
      })
      
      output$reactome_heatplot <- renderPlot({
        reactome_heatplot_input()
      })
      
      output$reactome_cnetplot <- renderPlot({
        reactome_cnetplot_input()
      })
      
      output$reactome_emaplot <- renderPlot({
        reactome_emaplot_input()
      })
      
      ### Download objects and functions ### ------------------------------------
      datasetInput_for_reactome <- reactive({
        table_for_reactome = res_reactome()
        switch(input$dataset_for_reactome,
               "full_results" = res_reactome()$all_table,
               "significant_results" = res_reactome()$sig_table)
      })
      
      output$downloadreactome <- downloadHandler(
        filename = function() { paste(input$dataset_for_reactome, ".txt", sep = "") },
        content = function(file) {
          write.table(datasetInput_for_reactome(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t") }
      )
      
      output$download_reactome_barplot <- downloadHandler(
        filename = 'barplot.pdf',
        content = function(file) {
          pdf(file, width = input$reactome_wide_bar, height = input$reactome_high_bar)
          print(reactome_barplot_input())
          dev.off()
        }
      )
      
      output$download_reactome_dotplot <- downloadHandler(
        filename = 'dotplot.pdf',
        content = function(file) {
          pdf(file, width = input$reactome_wide_dot, height = input$reactome_high_dot)
          print(reactome_dotplot_input())
          dev.off()
        }
      )
      
      # output$download_reactome_dotplot_opt <- downloadHandler(
      #   filename = 'dotplot_opt.pdf',
      #   content = function(file) {
      #     pdf(file, width = input$reactome_wide_dot_opt, height = input$reactome_high_dot_opt)
      #     print(reactome_dotplot_opt_input())
      #     dev.off()
      #   }
      # )
      
      output$download_reactome_heatplot <- downloadHandler(
        filename = 'heatplot.pdf',
        content = function(file) {
          pdf(file, width = input$reactome_wide_heat, height = input$reactome_high_heat)
          print(reactome_heatplot_input())
          dev.off()
        }
      )
      
      output$download_reactome_cnetplot <- downloadHandler(
        filename = 'cnetplot.pdf',
        content = function(file) {
          pdf(file, width = input$reactome_wide_cnet, height = input$reactome_high_cnet)
          print(reactome_cnetplot_input())
          dev.off()
        }
      )
      
      output$download_reactome_emaplot <- downloadHandler(
        filename = 'emaplot.pdf',
        content = function(file) {
          pdf(file, width = input$reactome_wide_ema, height = input$reactome_high_ema)
          print(reactome_emaplot_input())
          dev.off()
        }
      )
    }
    
    
  )
  
  ## observeEvent for GSEAGO ##--------
  observeEvent(
    input$analyze_for_gsego,{
      ### Interactive UI functions ### 
      output$downloadTable_gsego <- renderUI({
        selectizeInput("dataset_for_gsego",
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })
      
      output$downloadButton_gsego <- renderUI({
        downloadButton('downloadgsego', 'Save table', class = "downloadgsego")
      })
      
      if(!input$import_from_genelist_tool_for_gsego) {
        genelist <- reactive({ strsplit(input$text_for_gsego,'\n')[[1]]})
        gene_name <-reactive({unique(unlist(strsplit(genelist(),";")[]))})
      }else{
        genelist <- reactive({ get(input$genelist_for_gsego,envir = get_all_genelist_Env()) })
        gene_name <- reactive({unique(unlist(strsplit(genelist()$symbol,";")[]))})
      }
      
      # genelist <- reactive({ strsplit(input$text_input_for_gsego,'\n')[[1]] })
      # gene_name <-reactive(unlist(strsplit(genelist(),";")))  
      
      gene_df_gsego <- reactive({
        gene_df <- data.frame(name=gene_name())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        colnames(gene_df) = c("name","fc")
        gene_df = as.data.frame(gene_df)
        gene_df$fc = as.numeric(as.character(gene_df$fc))
        gene_df$name = as.character(gene_df$name)
        gene_df
      })      
      # }  
      
      # if(input$import_from_for_gsego) {
      #   if(input$import_for_gsego == "DEP-LFQ") {
      #     gene_df_gsego <<- get_results(dep()) %>% dplyr::select(name, paste(input$import_contrast_for_gsego, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsego, "ratio", sep = "_"))
      #   } else {
      #     if(input$import_for_gsego == "DEG-RNAseq") {
      #       gene_df_gsego <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::select(name, paste(input$import_contrast_for_gsego, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsego, "log2FoldChange", sep = "_"))
      #     }
      #   }
      # }
      
      # if(!input$import_from_for_gsego) {
      gene_df_gsego_check <- try(gene_df_gsego())
      # }
      
      # if(!input$import_from_for_gsego) {
      if(!class(gene_df_gsego_check) == "try-error") {
        gene_df_gsego_check$name = rm_digit_end(gene_df_gsego_check$name)
        gene_df_gsego <- gene_df_gsego_check
      }
      # }
      
      # if(input$import_from_for_gsego) {
      #   gene_df_gsego$name = rm_digit_end(gene_df_gsego$name)
      # }
      
      pkg_for_gsego <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_gsego]
      
      test <- gene_df_gsego
      
      # reat_gsego <<- reactive({
      #Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Please wait ...", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      reat_gsego <- try(gsegoAnalysis(df = gene_df_gsego, organism = input$organism_for_gsego, species_df = annoSpecies_df), silent = TRUE)
      # })
      
      # test1 <<- reat_gsego()
      res_gsego <- reactive({
        #   validate(
        #   need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
        #   need(class(reat_gsego()) != "try-error", "Your order ranked geneList can not do gsea analysis, Please go to GO panel if desired")
        # )
        try(givegseGO_res_and_table(reat = reat_gsego, ont = input$gsego_ont, pCutoff = input$gsego_p, p.adj.cutoff = input$gsego_padj, NES.cutoff = input$gsego_NES, simplify = input$gsego_simplify, Phenotype = input$gsego_Phenotype), silent = TRUE)
      })
      
      test2 <- - res_gsego()
      
      output$gsego_Table <- DT::renderDataTable({
        # if(input$import_from_for_gsego) {
        #   shiny::validate(
        #     need(nrow(gene_df_gsego) != 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
        #   )      
        # }
        
        # if(!input$import_from_for_gsego) {
        shiny::validate(
          need(class(gene_df_gsego_check) != "try-error", message = "No genes meet your requirements, and can not do the gseGO analysis")
        )      
        # }
        
        shiny::validate(need(require(pkg_for_gsego, character.only = TRUE), message = paste0("The package ", pkg_for_gsego, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_gsego, "')")))
        
        validate(
          need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
          need(reat_gsego, "Sorry, I can not do the gseGO analysis")
        )
        
        validate(
          need(nrow(res_gsego()$all_table) != 0, "Sorry, I can not do the gseGO analysis")
        )      
        DT::datatable(res_gsego()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
      })#selection = list(selected = c(1))
      
      
      output$significantBox_for_gsego <- renderInfoBox({
        # if(input$import_from_for_gsego) {
        #   shiny::validate(
        #     need(nrow(gene_df_gsego) != 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
        #   )      
        # }
        
        # if(!input$import_from_for_gsego) {
        shiny::validate(
          need(class(gene_df_gsego_check) != "try-error", message = "No genes meet your requirements, and can not do the gseGO analysis")
        )      
        # }
        
        shiny::validate(need(require(pkg_for_gsego, character.only = TRUE), message = paste0("The package ", pkg_for_gsego, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_gsego, "')")))
        
        validate(
          need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
          need(reat_gsego, "Sorry, I can not do the gseGO analysis")
        )
        
        validate(
          need(nrow(res_gsego()$all_table) != 0, "Sorry, I can not do the gseGO analysis")
        )      
        
        num_total_gsego <- res_gsego()$all_table %>%
          nrow()
        num_signif_gsego <- res_gsego()$sig_table %>%
          nrow()
        frac_gsego <- num_signif_gsego / num_total_gsego
        
        if(frac_gsego == 0) {
          info_box_gsego <- infoBox("Significant terms",
                                    paste0(num_signif_gsego,
                                           " out of ",
                                           num_total_gsego),
                                    "No terms enriched",
                                    icon = icon("thumbs-down", lib = "glyphicon"),
                                    color = "red",
                                    width = 4)
        }
        if(!frac_gsego == 0) {
          info_box_gsego <-     infoBox("Significant terms",
                                        paste0(num_signif_gsego,
                                               " out of ",
                                               num_total_gsego),
                                        paste0(signif(frac_gsego * 100, digits = 3),
                                               "% of terms enriched"),
                                        icon = icon("thumbs-up", lib = "glyphicon"),
                                        color = "green",
                                        width = 4)
        }
        info_box_gsego
      })
      
      output$gsego_term <- renderUI({
        if (!nrow(res_gsego()$sig_table) == 0) {
          selectizeInput("gsego_term",
                         "Term",
                         choices = res_gsego()$sig_table$Description, selected = res_gsego()$sig_table$Description[1], multiple = FALSE)
        }
      })    
      
      gsego_barplot_input <- reactive({
        # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        gse_barplot(res = res_gsego(), ShowCategory = input$gsego_ShowCategory_bar, color = input$gsego_color, ont = input$gsego_ont, Split = input$gsego_bar_if_Split_for_ont_ALL)
      })
      
      gsego_dotplot_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        gse_dotplot(res = res_gsego(), ShowCategory = input$gsego_ShowCategory_dot, color = input$gsego_color, ont = input$gsego_ont, Split = input$gsego_dot_if_Split_for_ont_ALL)
      })
      
      gsego_heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          try(Heatplot(res_gsego()$sig_res, showCategory = input$gsego_ShowCategory_heat, foldChange = res_gsego()$de), silent = T)  
        })
      })
      
      gsego_cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          try(Cnetplot(x = res_gsego()$sig_res, showCategory = input$gsego_ShowCategory_cnet, foldChange = res_gsego()$de, circular = input$gsego_circular_cnet, colorEdge = TRUE), silent = T)
        })
      })
      
      gsego_emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          try(Emapplot(res_gsego()$sig_res, showCategory = input$gsego_ShowCategory_ema, color = input$gsego_color, layout = "kk"), silent = T)
        })
      })
      
      gsego_Gseaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          my_gseaplot2(res_gsego()$sig_res, geneSetID = match(input$gsego_term, res_gsego()$sig_res$Description), title = input$gsego_term,# geneSetID = input$gsego_Table_rows_selected, title = res_gsego()$sig_table[input$gsego_Table_rows_selected,3]
                       color = "green",
                       base_size = 11,
                       rel_heights = c(1.5, 0.5, 1),
                       subplots = 1:3,
                       pvalue_table = TRUE,
                       ES_geom = "line")        
        })
      })
      
      output$gsego_barplot <- renderPlot({
        gsego_barplot_input()
      })
      
      output$gsego_dotplot <- renderPlot({
        gsego_dotplot_input()
      })
      
      output$gsego_heatplot <- renderPlot({
        gsego_heatplot_input()
      })
      
      output$gsego_cnetplot <- renderPlot({
        gsego_cnetplot_input()
      })
      
      output$gsego_emaplot <- renderPlot({
        gsego_emaplot_input()
      })
      
      output$gsego_Gseaplot <- renderPlot({
        gsego_Gseaplot_input()
      })
      
      ### Download objects and functions ### ------------------------------------
      datasetInput_for_gsego <- reactive({
        table_for_gsego = res_gsego()
        if(input$gsego_ont == "ALL") {
          if(input$gsego_bar_if_Split_for_ont_ALL | input$gsego_dot_if_Split_for_ont_ALL){
            switch(input$dataset_for_gsego,
                   "full_results" = res_gsego()$all_table %>% dplyr::arrange(ONTOLOGY),
                   "significant_results" = res_gsego()$sig_table %>% dplyr::arrange(ONTOLOGY))
            
          } else {
            switch(input$dataset_for_gsego,
                   "full_results" = res_gsego()$all_table,
                   "significant_results" = res_gsego()$sig_table)
            
          }
          
        } else {
          switch(input$dataset_for_gsego,
                 "full_results" = res_gsego()$all_table,
                 "significant_results" = res_gsego()$sig_table)
        }
        # switch(input$dataset_for_gsego,
        #        "full_results" = res_gsego()$all_table,
        #        "significant_results" = res_gsego()$sig_table)
      })
      
      output$downloadgsego <- downloadHandler(
        filename = function() { paste(input$dataset_for_gsego, ".txt", sep = "") },
        content = function(file) {
          write.table(datasetInput_for_gsego(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t") }
      )
      
      output$download_gsego_barplot <- downloadHandler(
        filename = 'barplot.pdf',
        content = function(file) {
          pdf(file, width = input$gsego_wide_bar, height = input$gsego_high_bar)
          print(gsego_barplot_input())
          dev.off()
        }
      )
      
      output$download_gsego_dotplot <- downloadHandler(
        filename = 'dotplot.pdf',
        content = function(file) {
          pdf(file, width = input$gsego_wide_dot, height = input$gsego_high_dot)
          print(gsego_dotplot_input())
          dev.off()
        }
      )
      
      output$download_gsego_heatplot <- downloadHandler(
        filename = 'heatplot.pdf',
        content = function(file) {
          pdf(file, width = input$gsego_wide_heat, height = input$gsego_high_heat)
          print(gsego_heatplot_input())
          dev.off()
        }
      )
      
      output$download_gsego_cnetplot <- downloadHandler(
        filename = 'cnetplot.pdf',
        content = function(file) {
          pdf(file, width = input$gsego_wide_cnet, height = input$gsego_high_cnet)
          print(gsego_cnetplot_input())
          dev.off()
        }
      )
      
      output$download_gsego_emaplot <- downloadHandler(
        filename = 'emaplot.pdf',
        content = function(file) {
          pdf(file, width = input$gsego_wide_ema, height = input$gsego_high_ema)
          print(gsego_emaplot_input())
          dev.off()
        }
      )
      
      output$download_gsego_Gseaplot <- downloadHandler(
        filename = 'Gseaplot.pdf',
        content = function(file) {
          pdf(file, width = input$gsego_wide_Gsea, height = input$gsego_high_Gsea)
          print(gsego_Gseaplot_input())
          dev.off()
        }
      )
    }
  )
  
 
  ## observeEvent for gsekegg ##------
  observeEvent(
    input$analyze_for_gsekegg,{
    ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_gsekegg <- renderUI({
      selectizeInput("dataset_for_gsekegg",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                     ))
    })
    
    output$downloadButton_gsekegg <- renderUI({
      downloadButton('downloadgsekegg', 'Save table', class = "downloadgsekegg")
    })
    
    if(!input$import_from_genelist_tool_for_gsekegg) {
      genelist <- reactive({ strsplit(input$text_for_gsekegg,'\n')[[1]]})
      gene_name <-reactive({unique(unlist(strsplit(genelist(),";")[]))})
    }else{
      genelist <- reactive({ get(input$genelist_for_gsekegg,envir = get_all_genelist_Env()) })
      gene_name <- reactive({unique(unlist(strsplit(genelist()$symbol,";")[]))})
    }
    
    # genelist <- reactive({ strsplit(input$text_input_for_gsekegg,'\n')[[1]] })
    # gene_name <-reactive(unlist(strsplit(genelist(),";")))  
    
    gene_df_gsekegg <- reactive({
      gene_df <- data.frame(name=gene_name())
      gene_df$name <-  sapply(gene_df$name,the1stname)
      gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
      rownames(gene_df) = NULL
      colnames(gene_df) = c("name","fc")
      gene_df = as.data.frame(gene_df)
      gene_df$fc = as.numeric(as.character(gene_df$fc))
      gene_df$name = as.character(gene_df$name)
      gene_df
    })      
    # }
    
    # if(input$import_from_for_gsekegg) {
    #   if(input$import_for_gsekegg == "DEP-LFQ") {
    #     gene_df_gsekegg <<- get_results(dep()) %>% dplyr::select(name, paste(input$import_contrast_for_gsekegg, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsekegg, "ratio", sep = "_"))
    #   } else {
    #     if(input$import_for_gsekegg == "DEG-RNAseq") {
    #       gene_df_gsekegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::select(name, paste(input$import_contrast_for_gsekegg, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsekegg, "log2FoldChange", sep = "_"))
    #     }
    #   }
    # }
    
    # if(!input$import_from_for_gsekegg) {
    gene_df_gsekegg_check <- try(gene_df_gsekegg())
    # }
    
    # if(!input$import_from_for_gsekegg) {
    if(!class(gene_df_gsekegg_check) == "try-error") {
      gene_df_gsekegg_check$name = rm_digit_end(gene_df_gsekegg_check$name)
      gene_df_gsekegg <- gene_df_gsekegg_check
    }
    # }
    
    # if(input$import_from_for_gsekegg) {
    #   gene_df_gsekegg$name = rm_digit_end(gene_df_gsekegg$name)
    # }
    
    pkg_for_gsekegg <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_gsekegg]
    
    # reat_gsekegg <<- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Please wait ...", value = 0.66)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    reat_gsekegg <- try(gsekeggAnalysis(df = gene_df_gsekegg, organism = input$organism_for_gsekegg, species_df = annoSpecies_df), silent = TRUE)
    # })
    
    res_gsekegg <- reactive({
      try(give_gsekegg_gsereat_res_and_table(reat = reat_gsekegg, pCutoff = input$gsekegg_p, p.adj.cutoff = input$gsekegg_padj, NES.cutoff = input$gsekegg_NES, Phenotype = input$gsekegg_Phenotype), silent = TRUE)
    })
    
    test2 <<- res_gsekegg()
    
    output$gsekegg_Table <- DT::renderDataTable({
      # if(input$import_from_for_gsekegg) {
      #   shiny::validate(
      #     need(nrow(gene_df_gsekegg) != 0, message = "No genes meet your requirements, and can not do the gsekegg analysis")
      #   )      
      # }
      
      # if(!input$import_from_for_gsekegg) {
      shiny::validate(
        need(class(gene_df_gsekegg_check) != "try-error", message = "No genes meet your requirements, and can not do the gsekegg analysis")
      )      
      # }
      
      shiny::validate(need(require(pkg_for_gsekegg, character.only = TRUE), message = paste0("The package ", pkg_for_gsekegg, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_gsekegg, "')")))
      
      validate(
        need(length(input$gsekegg_Phenotype) != 0, "Please choose at least one Phenotype"),
        need(reat_gsekegg, "Sorry, I can not do the gsekegg analysis")
      )
      validate(
        need(nrow(res_gsekegg()$all_table) != 0, "Sorry, I can not do the gsekegg analysis")
      )                
      DT::datatable(res_gsekegg()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
    })#selection = list(selected = c(1))
    
    output$significantBox_for_gsekegg <- renderInfoBox({
      # if(!input$import_from_for_gsekegg) {
      shiny::validate(
        need(class(gene_df_gsekegg_check) != "try-error", message = "No genes meet your requirements, and can not do the gsekegg analysis")
      )      
      # }
      
      # if(input$import_from_for_gsekegg) {
      #   shiny::validate(
      #     need(nrow(gene_df_gsekegg) != 0, message = "No genes meet your requirements, and can not do the gsekegg analysis")
      #   )      
      # }
      
      shiny::validate(need(require(pkg_for_gsekegg, character.only = TRUE), message = paste0("The package ", pkg_for_gsekegg, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_gsekegg, "')")))
      
      validate(
        need(length(input$gsekegg_Phenotype) != 0, "Please choose at least one Phenotype"),
        need(reat_gsekegg, "Sorry, I can not do the gsekegg analysis")
      ) 
      validate(
        need(nrow(res_gsekegg()$all_table) != 0, "Sorry, I can not do the gsekegg analysis")
      )             
      num_total_gsekegg <- res_gsekegg()$all_table %>%
        nrow()
      num_signif_gsekegg <- res_gsekegg()$sig_table %>%
        nrow()
      frac_gsekegg <- num_signif_gsekegg / num_total_gsekegg
      
      if(frac_gsekegg == 0) {
        info_box_gsekegg <- infoBox("Significant terms",
                                    paste0(num_signif_gsekegg,
                                           " out of ",
                                           num_total_gsekegg),
                                    "No terms enriched",
                                    icon = icon("thumbs-down", lib = "glyphicon"),
                                    color = "red",
                                    width = 4)
      }
      if(!frac_gsekegg == 0) {
        info_box_gsekegg <-     infoBox("Significant terms",
                                        paste0(num_signif_gsekegg,
                                               " out of ",
                                               num_total_gsekegg),
                                        paste0(signif(frac_gsekegg * 100, digits = 3),
                                               "% of terms enriched"),
                                        icon = icon("thumbs-up", lib = "glyphicon"),
                                        color = "green",
                                        width = 4)
      }
      info_box_gsekegg
    })
    
    output$gsekegg_term <- renderUI({
      if (!nrow(res_gsekegg()$sig_table) == 0) {
        selectizeInput("gsekegg_term",
                       "Term",
                       choices = res_gsekegg()$sig_table$Description, selected = res_gsekegg()$sig_table$Description[1], multiple = FALSE)
      }
    })    
    
    gsekegg_barplot_input <- reactive({
      # # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Plotting", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      enrichplot:::barplot.enrichResult(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_bar, color = input$gsekegg_color, x = "NES", split="phenotype") + labs(y = "NES")
    })
    
    gsekegg_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
      progress$set(message = "Plotting", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      enrichplot:::dotplot(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_dot, color = input$gsekegg_color, x = "NES", split="phenotype")
    })
    
    gsekegg_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        try(Heatplot(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_heat, foldChange = res_gsekegg()$de), silent = T)  
      })
    })
    
    gsekegg_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        try(Cnetplot(x = res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_cnet, foldChange = res_gsekegg()$de, circular = input$gsekegg_circular_cnet, colorEdge = TRUE), silent = T)
      })
    })
    
    gsekegg_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        try(Emapplot(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_ema, color = input$gsekegg_color, layout = "kk"), silent = T)
      })
    })
    
    gsekegg_Gseaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        my_gseaplot2(res_gsekegg()$sig_res, geneSetID = match(input$gsekegg_term, res_gsekegg()$sig_res$Description), title = input$gsekegg_term,# geneSetID = input$gsekegg_Table_rows_selected, title = res_gsekegg()$sig_table[input$gsekegg_Table_rows_selected,3]
                     color = "green",
                     base_size = 11,
                     rel_heights = c(1.5, 0.5, 1),
                     subplots = 1:3,
                     pvalue_table = TRUE,
                     ES_geom = "line")        
      })
    })
    
    output$gsekegg_barplot <- renderPlot({
      gsekegg_barplot_input()
    })
    
    output$gsekegg_dotplot <- renderPlot({
      gsekegg_dotplot_input()
    })
    
    output$gsekegg_heatplot <- renderPlot({
      gsekegg_heatplot_input()
    })
    
    output$gsekegg_cnetplot <- renderPlot({
      gsekegg_cnetplot_input()
    })
    
    output$gsekegg_emaplot <- renderPlot({
      gsekegg_emaplot_input()
    })
    
    output$gsekegg_Gseaplot <- renderPlot({
      gsekegg_Gseaplot_input()
    })
    
    ### Download objects and functions ### ------------------------------------
    datasetInput_for_gsekegg <- reactive({
      table_for_gsekegg = res_gsekegg()
      switch(input$dataset_for_gsekegg,
             "full_results" = res_gsekegg()$all_table,
             "significant_results" = res_gsekegg()$sig_table)
    })
    
    output$downloadgsekegg <- downloadHandler(
      filename = function() { paste(input$dataset_for_gsekegg, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_gsekegg(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )
    
    output$download_gsekegg_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_bar, height = input$gsekegg_high_bar)
        print(gsekegg_barplot_input())
        dev.off()
      }
    )
    
    output$download_gsekegg_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_dot, height = input$gsekegg_high_dot)
        print(gsekegg_dotplot_input())
        dev.off()
      }
    )
    
    output$download_gsekegg_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_heat, height = input$gsekegg_high_heat)
        print(gsekegg_heatplot_input())
        dev.off()
      }
    )
    
    output$download_gsekegg_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_cnet, height = input$gsekegg_high_cnet)
        print(gsekegg_cnetplot_input())
        dev.off()
      }
    )
    
    output$download_gsekegg_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_ema, height = input$gsekegg_high_ema)
        print(gsekegg_emaplot_input())
        dev.off()
      }
    )
    
    output$download_gsekegg_Gseaplot <- downloadHandler(
      filename = 'Gseaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_Gsea, height = input$gsekegg_high_Gsea)
        print(gsekegg_Gseaplot_input())
        dev.off()
      }
    )
  })
  
  ## observeEvent for gsereactome ##------
  observeEvent(
    input$analyze_for_gsereactome,{
    ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_gsereactome <- renderUI({
      selectizeInput("dataset_for_gsereactome",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                     ))
    })
    
    output$downloadButton_gsereactome <- renderUI({
      downloadButton('downloadgsereactome', 'Save table', class = "downloadgsereactome")
    })
    
    if(!input$import_from_genelist_tool_for_gsereactome) {
      genelist <- reactive({ strsplit(input$text_for_gsereactome,'\n')[[1]]})
      gene_name <-reactive({unique(unlist(strsplit(genelist(),";")[]))})
    }else{
      genelist <- reactive({ get(input$genelist_for_gsereactome,envir = get_all_genelist_Env()) })
      gene_name <- reactive({unique(unlist(strsplit(genelist()$symbol,";")[]))})
    }
    # genelist <- reactive({ strsplit(input$text_input_for_gsereactome,'\n')[[1]] })
    # gene_name <-reactive(unlist(strsplit(genelist(),";")))  
    
    gene_df_gsereactome <- reactive({
      gene_df <- data.frame(name=gene_name())
      gene_df$name <-  sapply(gene_df$name,the1stname)
      gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
      rownames(gene_df) = NULL
      colnames(gene_df) = c("name","fc")
      gene_df = as.data.frame(gene_df)
      gene_df$fc = as.numeric(as.character(gene_df$fc))
      gene_df$name = as.character(gene_df$name)
      gene_df
    })      
    # }
    
    # if(input$import_from_for_gsereactome) {
    #   if(input$import_for_gsereactome == "DEP-LFQ") {
    #     gene_df_gsereactome <<- get_results(dep()) %>% dplyr::select(name, paste(input$import_contrast_for_gsereactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsereactome, "ratio", sep = "_"))
    #   } else {
    #     if(input$import_for_gsereactome == "DEG-RNAseq") {
    #       gene_df_gsereactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::select(name, paste(input$import_contrast_for_gsereactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsereactome, "log2FoldChange", sep = "_"))
    #     }
    #   }
    # }
    
    # if(!input$import_from_for_gsereactome) {
    gene_df_gsereactome_check <- try(gene_df_gsereactome())
    # }
    
    # if(!input$import_from_for_gsereactome) {
    if(!class(gene_df_gsereactome_check) == "try-error") {
      gene_df_gsereactome_check$name = rm_digit_end(gene_df_gsereactome_check$name)
      gene_df_gsereactome <- gene_df_gsereactome_check
    }
    # }
    
    # if(input$import_from_for_gsereactome) {
    #   gene_df_gsereactome$name = rm_digit_end(gene_df_gsereactome$name)
    # }
    
    pkg_for_gsereactome <- annoSpecies_df$pkg[annoSpecies_df$species == input$organism_for_gsereactome]              
    
    # reat_gsereactome <<- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Please wait ...", value = 0.66)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    reat_gsereactome <- try(gsereactAnalysis(df = gene_df_gsereactome, organism = input$organism_for_gsereactome, species_df = annoSpecies_df_for_reactome), silent = TRUE)
    # })
    
    res_gsereactome <- reactive({
      try(give_gsekegg_gsereat_res_and_table(reat = reat_gsereactome, pCutoff = input$gsereactome_p, p.adj.cutoff = input$gsereactome_padj, NES.cutoff = input$gsereactome_NES, Phenotype = input$gsereactome_Phenotype), silent = TRUE)
    })
    
    # test2 <<- res_gsereactome()
    
    output$gsereactome_Table <- DT::renderDataTable({
      # if(input$import_from_for_gsereactome) {
      #   shiny::validate(
      #     need(nrow(gene_df_gsereactome) != 0, message = "No genes meet your requirements, and can not do the gseReactome analysis")
      #   )      
      # }
      
      # if(!input$import_from_for_gsereactome) {
      shiny::validate(
        need(class(gene_df_gsereactome_check) != "try-error", message = "No genes meet your requirements, and can not do the gseReactome analysis")
      )      
      # }
      
      shiny::validate(need(require(pkg_for_gsereactome, character.only = TRUE), message = paste0("The package ", pkg_for_gsereactome, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_gsereactome, "')")))
      
      validate(
        need(length(input$gsereactome_Phenotype) != 0, "Please choose at least one Phenotype"),
        need(reat_gsereactome, "Sorry, I can not do the gsereactome analysis")
      ) 
      validate(
        need(nrow(res_gsereactome()$all_table) != 0, "Sorry, I can not do the gsereactome analysis")
      )             
      DT::datatable(res_gsereactome()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
    })#selection = list(selected = c(1))
    
    output$significantBox_for_gsereactome <- renderInfoBox({
      # if(input$import_from_for_gsereactome) {
      #   shiny::validate(
      #     need(nrow(gene_df_gsereactome) != 0, message = "No genes meet your requirements, and can not do the gseReactome analysis")
      #   )      
      # }
      
      # if(!input$import_from_for_gsereactome) {
      shiny::validate(
        need(class(gene_df_gsereactome_check) != "try-error", message = "No genes meet your requirements, and can not do the gseReactome analysis")
      )      
      # }   
      
      shiny::validate(need(require(pkg_for_gsereactome, character.only = TRUE), message = paste0("The package ", pkg_for_gsereactome, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_gsereactome, "')")))
      
      validate(
        need(length(input$gsereactome_Phenotype) != 0, "Please choose at least one Phenotype"),
        need(reat_gsereactome, "Sorry, I can not do the gsereactome analysis")
      )
      validate(
        need(nrow(res_gsereactome()$all_table) != 0, "Sorry, I can not do the gsereactome analysis")
      )
      num_total_gsereactome <- res_gsereactome()$all_table %>%
        nrow()
      num_signif_gsereactome <- res_gsereactome()$sig_table %>%
        nrow()
      frac_gsereactome <- num_signif_gsereactome / num_total_gsereactome
      
      if(frac_gsereactome == 0) {
        info_box_gsereactome <- infoBox("Significant terms",
                                        paste0(num_signif_gsereactome,
                                               " out of ",
                                               num_total_gsereactome),
                                        "No terms enriched",
                                        icon = icon("thumbs-down", lib = "glyphicon"),
                                        color = "red",
                                        width = 4)
      }
      if(!frac_gsereactome == 0) {
        info_box_gsereactome <-     infoBox("Significant terms",
                                            paste0(num_signif_gsereactome,
                                                   " out of ",
                                                   num_total_gsereactome),
                                            paste0(signif(frac_gsereactome * 100, digits = 3),
                                                   "% of terms enriched"),
                                            icon = icon("thumbs-up", lib = "glyphicon"),
                                            color = "green",
                                            width = 4)
      }
      info_box_gsereactome
    })
    
    output$gsereactome_term <- renderUI({
      if (!nrow(res_gsereactome()$sig_table) == 0) {
        selectizeInput("gsereactome_term",
                       "Term",
                       choices = res_gsereactome()$sig_table$Description, selected = res_gsereactome()$sig_table$Description[1], multiple = FALSE)
      }
    })    
    
    gsereactome_barplot_input <- reactive({
      # # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Plotting", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      enrichplot:::barplot.enrichResult(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_bar, color = input$gsereactome_color, x = "NES", split="phenotype") + labs(y = "NES")
    })
    
    gsereactome_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
      progress$set(message = "Plotting", value = 0.66)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      enrichplot:::dotplot(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_dot, color = input$gsereactome_color, x = "NES", split="phenotype")
    })
    
    gsereactome_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        try(Heatplot(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_heat, foldChange = res_gsereactome()$de), silent = T)  
      })
    })
    
    gsereactome_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        try(Cnetplot(x = res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_cnet, foldChange = res_gsereactome()$de, circular = input$gsereactome_circular_cnet, colorEdge = TRUE), silent = T)
      })
    })
    
    gsereactome_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        try(Emapplot(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_ema, color = input$gsereactome_color, layout = "kk"), silent = T)
      })
    })
    
    gsereactome_Gseaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        my_gseaplot2(res_gsereactome()$sig_res, geneSetID = match(input$gsereactome_term, res_gsereactome()$sig_res$Description), title = input$gsereactome_term,# geneSetID = input$gsereactome_Table_rows_selected, title = res_gsereactome()$sig_table[input$gsereactome_Table_rows_selected,3]
                     color = "green",
                     base_size = 11,
                     rel_heights = c(1.5, 0.5, 1),
                     subplots = 1:3,
                     pvalue_table = TRUE,
                     ES_geom = "line")        
      })
    })
    
    output$gsereactome_barplot <- renderPlot({
      gsereactome_barplot_input()
    })
    
    output$gsereactome_dotplot <- renderPlot({
      gsereactome_dotplot_input()
    })
    
    output$gsereactome_heatplot <- renderPlot({
      gsereactome_heatplot_input()
    })
    
    output$gsereactome_cnetplot <- renderPlot({
      gsereactome_cnetplot_input()
    })
    
    output$gsereactome_emaplot <- renderPlot({
      gsereactome_emaplot_input()
    })
    
    output$gsereactome_Gseaplot <- renderPlot({
      gsereactome_Gseaplot_input()
    })
    
    ### Download objects and functions ### ------------------------------------
    datasetInput_for_gsereactome <- reactive({
      table_for_gsereactome = res_gsereactome()
      switch(input$dataset_for_gsereactome,
             "full_results" = res_gsereactome()$all_table,
             "significant_results" = res_gsereactome()$sig_table)
    })
    
    output$downloadgsereactome <- downloadHandler(
      filename = function() { paste(input$dataset_for_gsereactome, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_gsereactome(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )
    
    output$download_gsereactome_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_bar, height = input$gsereactome_high_bar)
        print(gsereactome_barplot_input())
        dev.off()
      }
    )
    
    output$download_gsereactome_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_dot, height = input$gsereactome_high_dot)
        print(gsereactome_dotplot_input())
        dev.off()
      }
    )
    
    output$download_gsereactome_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_heat, height = input$gsereactome_high_heat)
        print(gsereactome_heatplot_input())
        dev.off()
      }
    )
    
    output$download_gsereactome_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_cnet, height = input$gsereactome_high_cnet)
        print(gsereactome_cnetplot_input())
        dev.off()
      }
    )
    
    output$download_gsereactome_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_ema, height = input$gsereactome_high_ema)
        print(gsereactome_emaplot_input())
        dev.off()
      }
    )
    
    output$download_gsereactome_Gseaplot <- downloadHandler(
      filename = 'Gseaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_Gsea, height = input$gsereactome_high_Gsea)
        print(gsereactome_Gseaplot_input())
        dev.off()
      }
    )
  })
  
  ## for genelist ##----------
  Genelist_server_module("bbb",alpha=reactive({input$padj}),lfc=reactive({input$logFC}),object = Filtered_result,
                         condition_matched_frame=condition_matched_frame)
}



# server <-   function(input, output){}

## run app##-------
shinyApp(ui=ui, server=server)

