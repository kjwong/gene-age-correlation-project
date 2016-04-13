## Name: Kenny Wong
## Collaborators: Arjun Krishnan
## Murphy Lab

if (!require("shiny")) install.packages('shiny')
if (!require("DT")) install.packages('DT')
if (!require("d3heatmap")) install.packages('d3heatmap')
if (!require("plotly")) install.packages('plotly')
if (!require("shinyBS")) install.packages('shinyBS')

library(shiny)
library(DT)
library(d3heatmap)
library(shinyBS)
library(plotly)

shinyUI(fluidPage(
  
  h1("A GUI to visualize and analyze age-associated genes in humans"),
  hr(),
  h4("This app lets you visualize and analyze the genes most correlated with age in a specified age range."),
  h4("Upload (1) expression data and (2) a table of the samples and their respective ages, and hit run!"),
  hr(),

  sidebarLayout(
    sidebarPanel(
      h4("Select a sample expression file."),
      bsTooltip('file1',"Row names: Entrez Gene ID, Column names: Sample ID",
                placement="right",trigger="hover"),

      fileInput('file1',label=NULL,accept = c('.pcl')),
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
      uiOutput("upload2"),
      conditionalPanel(
        condition = "input.upload2",
        tags$hr()
      ),
      uiOutput("filters"),
      conditionalPanel(
        condition = "input.upload2",
        helpText("If left blank, all samples will be used."),
        hr(),
        actionButton("run_samples", "Run samples")
      ),
      hr(),
      helpText("Content by Murphy Lab at Princeton University")
    ),
    
    mainPanel(
      uiOutput("nosamp"),
      conditionalPanel(
        condition = "input.run_samples",
        h3("Sample count by age:"),
        helpText("The age distribution of your filtered samples."),
        plotlyOutput("plot"),
        uiOutput("plot_caption"),
        hr(),
        fluidRow(
          column(6,
                  uiOutput("slider_plot")
            ),
          column(6,
                 numericInput("numruns",label="Choose number of bootstrap runs [1,50]:",value=10,min=1,max=50),
                 bsTooltip('numruns',"Increasing the number of runs will yield more accurate correlation scores but will take longer.",
                           placement="bottom",trigger="hover")
          )
            
        ),
        
        actionButton("runpcl","Run gene expression-age correlation"),
        hr()
      ),
      
      uiOutput("nosamp2"),
      conditionalPanel(
        condition = "output.plot2",
        h3("Predictive gene count:"),
        helpText("Predictive genes are those most correlated with age within your selected age range.")
      ),
      conditionalPanel(
        condition = "input.runpcl",
        plotlyOutput("plot2"),
        uiOutput("plot2_caption"),
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
        h3("Age-associated genes:")
      ),
      conditionalPanel(
        condition = "input.tablepcl",
        h4("Gene lists"),
        tabsetPanel(
          id = 'dataset',
          tabPanel('Pos. correlation', DT::dataTableOutput('ptable'),
                   conditionalPanel(condition="output.ptable",
                                    downloadButton("ptable_dl","Download table"),
                                    hr(),
                                    helpText("Create a heat map for positively predictive genes."),
                                    fluidRow(
                                      column(6,
                                        uiOutput("heatposage")),
                                      column(6,
                                        textInput("userposgenes",label="Choose genes for heat map:",value=""),
                                        bsTooltip('userposgenes',"Enter 3 or more gene names separated by space.",
                                                  placement="bottom",trigger="hover"),
                                        helpText("If left blank, genes with the top 50 scores will be used."))
                                    ),
                                    actionButton("runposheat","Generate heat map"),
                                    hr(),
                                    helpText("Perform GO term analysis for positively predictive genes."),
                                    selectInput(inputId = "posstat", label = "Choose a test for p-value:", choices = c("fisher", "ks", "t")),
                                    actionButton("runposgt", "Generate GO terms")
                                    )
                   ),
                                    
          tabPanel('Neg. correlation', DT::dataTableOutput('ntable'),
                   conditionalPanel(condition="output.ntable",
                                    downloadButton("ntable_dl","Download table"),
                                    hr(),
                                    helpText("Create a heat map for negatively predictive genes."),
                                    fluidRow(
                                      column(6,
                                             
                                             uiOutput("heatnegage")),
                                      column(6,
                                             textInput("userneggenes",label="Choose genes for heat map:",value=""),
                                             bsTooltip('userneggenes',"Enter 3 or more gene names separated by space.",
                                                       placement="bottom",trigger="hover"),
                                             helpText("If left blank, genes with the top scores will be used."))
                                    ),
                                    actionButton("runnegheat","Generate heat map"),
                                    hr(),
                                    helpText("Perform GO term analysis for negatively predictive genes."),
                                    selectInput(inputId = "negstat", label = "Choose a test for p-value:", choices = c("fisher", "ks", "t")),
                                    actionButton("runneggt", "Generate GO terms")
                                    )
                   )
        ),
        hr()
      ),
      conditionalPanel(
        condition = "input.runposheat || input.runnegheat",
        h4("Expression value heatmap and clustering"),
        tabsetPanel(
          id = 'heatmap',
          tabPanel('Pos. correlation', d3heatmapOutput("posheat")),
          tabPanel('Neg. correlation', d3heatmapOutput("negheat"))
        ),
        hr()
      ),
      conditionalPanel(
        condition = "input.runposgt || input.runneggt",
        h4("Gene ontology analysis"),
        tabsetPanel(
          id = 'goterms',
          tabPanel('Pos. correlation',DT::dataTableOutput('pos_goterms'),
                   conditionalPanel(condition="output.pos_goterms",
                                    downloadButton("pgo_dl","Download table"),
                                    hr(),
                                    helpText("Induced subgraph of the most significant GO terms."),
                                    plotOutput("pos_go_graph"),
                                    sliderInput("posnodes",label="Number of significant nodes:",value=5,min=1,max=10))),
          tabPanel('Neg. correlation',DT::dataTableOutput('neg_goterms'),
                   conditionalPanel(condition="output.neg_goterms",
                                    downloadButton("ngo_dl","Download table"),
                                    hr(),
                                    helpText("Induced subgraph of the most significant GO terms."),
                                    plotOutput("neg_go_graph"),
                                    sliderInput("negnodes",label="Number of significant nodes:",value=5,min=1,max=10)))  
        ),
        hr()
      )
    )
  )
))