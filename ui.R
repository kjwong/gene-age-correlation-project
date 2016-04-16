## Name: Kenny Wong
## Collaborators: Arjun Krishnan
## Murphy Lab

if (!require("shiny")) install.packages('shiny')
if (!require("DT")) install.packages('DT')
if (!require("d3heatmap")) install.packages('d3heatmap')
if (!require("plotly")) install.packages('plotly')
if (!require("shinyBS")) install.packages('shinyBS')
if (!require("networkD3")) install.packages('networkD3')

library(shiny)
library(DT)
library(d3heatmap)
library(shinyBS)
library(plotly)
library(networkD3)

shinyUI(fluidPage(
  
  h1("Web interface for analyzing age-associated gene expression in humans"),
  hr(),
  h4("This app lets you visualize and analyze the genes most correlated with age in a specified age range."),
  hr(),
  
  sidebarLayout(
    sidebarPanel(
      h4("Select an expression data file."),
      bsTooltip('file1',"Row names: Entrez Gene ID, Column names: Sample ID",
                placement="right",trigger="hover"),
      fileInput('file1',label=NULL,accept = c('.pcl','.txt')),
      radioButtons('sep1', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), inline=FALSE, '\t'),
      checkboxInput('header1', 'Header', TRUE),
      tags$hr(),
      h4("Select a file with samples and their ages."),
      
      fileInput('file2',label=NULL,accept = c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        '.csv',
        '.tsv'
      )),
      bsTooltip('file2',"Row names: Sample ID, Column names: age (required), other (optional)",
                placement="right",
                trigger="hover"),
      radioButtons('sep2', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), inline=FALSE, ','),
      checkboxInput('header2', 'Header', TRUE),    
      hr(),
      h4("Select a default dataset."),
      uiOutput("defaultinfo"),
      selectInput('default',label=NULL,choices=c("blood.f","GSE58137"),
                  selected = "blood.f"),
      h6("blood.f: Whole blood samples from the Affymetrix expression platform (GPL570)."),
      h6("GSE58137: Whole blood samples from the Illumina Beadchip platform (GPL6947, GPL10558)."),
      hr(),
      helpText("If no files are selected, the default dataset will be used."),
     
      conditionalPanel(
        condition = "output.filters",
        hr()
      ),
      uiOutput("filters"),
      conditionalPanel(
        condition = "output.filters",
        helpText("If left blank, all samples will be used."),
        hr(),
        actionButton("run_samples", "Run analysis")
      ),
      hr(),
      h6("Created by the Murphy Lab at Princeton University")
    ),
    
    mainPanel(
      uiOutput("nosamp"),
      conditionalPanel(
        condition = "input.run_samples",
        h3("Sample count by age"),
        helpText("The age distribution of your filtered samples."),
        plotlyOutput("plot"),
        uiOutput("plot1_caption"),
        uiOutput("plot1_caption2"),
        conditionalPanel(
          condition = "output.plot",
          hr(),
          fluidRow(
            column(4,
                   uiOutput("slider_plot")
            ),
            column(4,
                   selectInput("corr",label="Choose correlation method:",choices=c("spearman","pearson"),selected="spearman")
            ),
            column(4,
                   numericInput("numruns",label="Choose number of runs [1,100]:",value=10,min=1,max=100),
                   bsTooltip('numruns',"Increasing the number of runs will yield more accurate correlation scores but will take longer.",
                             placement="bottom",trigger="hover")
            )
            
          ),
          
          actionButton("runpcl","Run correlation"),
          hr()
        )
      ),
      
      uiOutput("nosamp2"),
      conditionalPanel(
        condition = "output.plot2",
        h3("Predictive gene count"),
        helpText("Predictive genes are those most correlated with age within your selected age range.")
      ),
      conditionalPanel(
        condition = "input.runpcl",
        plotlyOutput("plot2"),
        uiOutput("plot2_caption"),
        uiOutput("plot2_caption2"),
        hr(),
        conditionalPanel(
          condition = "output.plot2",
          uiOutput("slider_plot2")
        )
      ),
      
      conditionalPanel(
        condition = "output.plot2",
        fluidRow(
          column(3,
                 actionButton("tablepcl","Select genes"))
        ),
        hr()
      ),
      conditionalPanel(
        condition = "input.tablepcl",
        h3("Age-associated genes")
      ),
      conditionalPanel(
        condition = "input.tablepcl",
        tabsetPanel(
          id = 'dataset',
          tabPanel('Positive', h4("Gene list"), helpText("Genes positively correlated with age."),DT::dataTableOutput('ptable'),
                   conditionalPanel(condition="output.ptable",
                                    downloadButton("ptable_dl","Download"),
                                    hr(),
                                    h4("Expression value heatmap and clustering"),
                                    helpText("Create a heat map for positively predictive genes."),
                                    fluidRow(
                                      column(4,
                                             uiOutput("heatposage")),
                                      column(4,
                                             numericInput("heatposnum",label="Choose a number of genes to display:",value=100,min=3,max=NA),
                                             bsTooltip('heatposnum',"minimum = 3",
                                                       placement="bottom",trigger="hover")),
                                      column(4,
                                             textInput("userposgenes",label="Choose specific genes to display:",value=""),
                                             bsTooltip('userposgenes',"Enter 3 or more gene names separated by space.",
                                                       placement="bottom",trigger="hover"),
                                             helpText("If left blank, all genes are used."))
                                    ),
                                    actionButton("runposheat","Generate heat map"),
                                    hr(),
                                    conditionalPanel(
                                      condition = "input.runposheat",
                                      d3heatmapOutput("posheat"),
                                      conditionalPanel(
                                        condition = "output.posheat",
                                        helpText("Click and drag to zoom in. Click to zoom out. Hover to see values.")
                                      ),
                                      hr()
                                    ),
                                    h4("Gene ontology analysis"),
                                    helpText("Perform GO term analysis for positively predictive genes."),
                                    fluidRow(
                                      column(4,
                                             bsTooltip('postype',"BP = Biological Process, MF = Molecular Function, CC = Cellular Component",
                                                       placement = "top",trigger="hover"),
                                             selectInput(inputId = "postype", label = "Choose a gene ontology:", choices = c("BP", "MF", "CC"),selected="BP")
                                      ),
                                      column(4,
                                          
                                             selectInput(inputId = "posstat", label = "Choose a statistical test:", choices = c("fisher", "ks")),
                                             bsTooltip('posstat',"fisher = Fisher exact test, ks = Kolmogorov-Smirnov test",
                                                       placement = "top",trigger="hover")
                                      )
                                    ),
                                    actionButton("runposgt", "Generate GO terms"),
                                    hr(),
                                    conditionalPanel(
                                      condition = "input.runposgt",
                                      DT::dataTableOutput('pos_goterms'),
                                      conditionalPanel(condition="output.pos_goterms",
                                                       downloadButton("pgo_dl","Download"),
                                                       hr(),
                                                       forceNetworkOutput("pos_go_graph"),
                                                       helpText("Significant GO terms (blue) are connected by their common ancestors (orange)."),
                                                       uiOutput("posn")))
                   )
          ),
          
          tabPanel('Negative', h4("Gene list"), helpText("Genes negatively correlated with age."),DT::dataTableOutput('ntable'),
                   conditionalPanel(condition="output.ntable",
                                    downloadButton("ntable_dl","Download"),
                                    hr(),
                                    h4("Expression value heatmap and clustering"),
                                    helpText("Create a heat map for negatively predictive genes."),
                                    fluidRow(
                                      column(4,
                                             uiOutput("heatnegage")),
                                      column(4,
                                             numericInput("heatnegnum",label="Choose a number of genes to display:",value=100,min=3,max=NA),
                                             bsTooltip('heatnegnum',"minimum = 3",
                                                       placement="bottom",trigger="hover")),
                                      column(4,
                                             textInput("userneggenes",label="Choose genes:",value=""),
                                             bsTooltip('userneggenes',"Enter 3 or more gene names separated by space.",
                                                       placement="bottom",trigger="hover"),
                                             helpText("If left blank, all genes are used."))
                                    ),
                                    actionButton("runnegheat","Generate heat map"),
                                    hr(),
                                    conditionalPanel(
                                      condition = "input.runnegheat",
                                      d3heatmapOutput("negheat"),
                                      conditionalPanel(
                                        condition = "output.negheat",
                                        helpText("Click and drag to zoom in. Click to zoom out. Hover to see values.")
                                      ),
                                      hr()
                                    ),
                                    h4("Gene ontology analysis"),
                                    helpText("Perform GO term analysis for negatively predictive genes."),
                                    fluidRow(
                                      column(4,
                                             bsTooltip('negtype',"BP = Biological Process, MF = Molecular Function, CC = Cellular Component",
                                                       placement = "top",trigger="hover"),
                                             selectInput(inputId = "negtype", label = "Choose a gene ontology:", choices = c("BP", "MF", "CC"),selected="BP")
                                      ),
                                      column(4,
                                             bsTooltip('negstat',"fisher = Fisher exact test, ks = Kolmogorov-Smirnov test",
                                                       placement = "top",trigger="hover"),
                                             selectInput(inputId = "negstat", label = "Choose a test for p-value:", choices = c("fisher", "ks"))
                                      )
                                    ),
                                    actionButton("runneggt", "Generate GO terms"),
                                    hr(),
                                    conditionalPanel(
                                      condition = "input.runneggt",
                                      DT::dataTableOutput('neg_goterms'),
                                      conditionalPanel(condition="output.neg_goterms",
                                                       downloadButton("ngo_dl","Download"),
                                                       hr(),
                                                       forceNetworkOutput("neg_go_graph"),
                                                       helpText("Significant GO terms (blue) are connected by their common ancestors (orange)."),
                                                       uiOutput("negn"),
                                                       hr()
                                      )
                                    )
                   )
          )
        )
      )
    )
  )
))