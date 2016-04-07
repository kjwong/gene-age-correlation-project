## Name: Kenny Wong
## Collaborators: Arjun Krishnan
## Murphy Lab

library(shiny)
library(DT)
library(d3heatmap)
library(shinyBS)
library(plotly)

shinyUI(fluidPage(
  
  h1("Web interface for gene-age correlation in humans"),
  hr(),
  h4("This app lets you visualize and analyze the genes most correlated with age in a specified age range."),
  h4("Upload (1) expression data and (2) a table of the samples and their respective ages, and hit run!"),
  hr(),

  sidebarLayout(
    sidebarPanel(
      h4("Select a GSM sample expression file."),
      bsTooltip('file1',"Row names: Entrez Gene ID, Column names: GSM accession number",
                placement="right",trigger="hover"),

      fileInput('file1',label=NULL,accept = c('.pcl')),
      radioButtons('sep1', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), inline=FALSE, '\t'),
      checkboxInput('header1', 'Header', TRUE),
      tags$hr(),
      h4("Select a file with GSM samples and their ages."),

      fileInput('file2',label=NULL,accept = c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        '.csv',
        '.tsv'
      )),
      bsTooltip('file2',"Row names: GSM accession number, Column names: age (required), other (optional)",
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
        h3("GSM sample count by age:"),
        plotlyOutput("plot"),
        fluidRow(
          column(6,
            uiOutput("slider_plot")),
          column(4,
                 hr(),
            textOutput("plot_caption"))
        ),
        helpText("Use the slider to select samples of an age range."),
        actionButton("runpcl","Run gene-age correlation"),
        hr()
      ),
      
      uiOutput("nosamp2"),
      conditionalPanel(
        condition = "output.plot2",
        h3("Predictive gene count:"),
        helpText("Predictive genes are those most correlated with age (Spearman rho of the highest magnitude).")
      ),
      conditionalPanel(
        condition = "input.runpcl",
        plotlyOutput("plot2"),
        conditionalPanel(
          condition = "output.plot2",
          fluidRow(
            column(6,
                   uiOutput("slider_plot2")),
            column(4,
                   hr(),
                   textOutput("plot2_caption"))
          )
        )
      ),

      conditionalPanel(
        condition = "output.plot2",
        helpText("Use the slider to set a correlation score threshold."),
        fluidRow(
          column(3,
            actionButton("tablepcl","See selected genes"))
        ),
        hr()
      ),
      conditionalPanel(
        condition = "input.tablepcl",
        h3("Your selected genes:")
      ),
      conditionalPanel(
        condition = "input.tablepcl",
        h4("Gene lists"),
        tabsetPanel(
          id = 'dataset',
          tabPanel('Pos. correlation', DT::dataTableOutput('ptable'),
                   conditionalPanel(condition="output.ptable",downloadButton("ptable_dl","Download table"))),
          tabPanel('Neg. correlation', DT::dataTableOutput('ntable'),
                   conditionalPanel(condition="output.ntable",downloadButton("ntable_dl","Download table"))),
          tabPanel('Pos. correlation (exp. values)', DT::dataTableOutput('ppcl'),
                   conditionalPanel(condition="output.ppcl",downloadButton("ppcl_dl","Download table"))),
          tabPanel('Neg. correlation (exp. values)', DT::dataTableOutput('npcl'),
                   conditionalPanel(condition="output.npcl",downloadButton("npcl_dl","Download table")))
        ),
        hr(),
        h4("Expression value heatmap and clustering"),
        tabsetPanel(
          id = 'heatmap',
          tabPanel('Pos. correlation', d3heatmapOutput("posheat")),
          tabPanel('Neg. correlation', d3heatmapOutput("negheat"))
        )
      ),
      conditionalPanel(
        condition = "output.ptable",
        selectInput(inputId = "stat", label = "Choose a statistical test (for p-value)", choices = c("fisher", "ks", "t")),
        actionButton("rungt", "Run GO term analysis"),
        helpText("See the most enriched GO terms in your selected genes."),
        hr()
      ),
      conditionalPanel(
        condition = "input.rungt",
        h3("GO term analysis:"),
        tabsetPanel(
          id = 'goterms',
          tabPanel('Pos. correlation',DT::dataTableOutput('pos_goterms'),
                   conditionalPanel(condition="output.pos_goterms",
                                    helpText("The 'classic' column denotes the p-value."),
                                    downloadButton("pgo_dl","Download table"))),
          tabPanel('Neg. correlation',DT::dataTableOutput('neg_goterms'),
                   conditionalPanel(condition="output.neg_goterms",
                                    helpText("The 'classic' column denotes the p-value."),
                                    downloadButton("ngo_dl","Download table")))  
        ),
        hr()
      )
    )
  )
))