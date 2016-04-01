library(shiny)
library(DT)
library(d3heatmap)
library(ggplot2)

shinyUI(fluidPage(
  
  titlePanel("Web interface for gene-age correlation in GSM samples"),
  hr(),
  h4("This app lets you visualize which genes are most correlated with age in a specified age range."),
  h4("Upload (1) GSM sample expression data and (2) a table of GSM samples and their ages, and hit run!"),
  hr(),

  sidebarLayout(
    sidebarPanel(
      helpText("Select a GSM sample expression file."),
      fileInput('file1', 'Choose file to upload',accept = c('.pcl')),
      radioButtons('sep1', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), inline=TRUE, '\t'),
      checkboxInput('header1', 'Header', TRUE),
      tags$hr(),
      helpText("Select a file with GSM samples and their ages."),
      fileInput('file2', 'Choose file to upload',accept = c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        '.csv',
        '.tsv'
      )),
      radioButtons('sep2', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), inline=TRUE, ','),
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
      conditionalPanel(
        condition = "input.run_samples",
        h3("GSM sample count by age:"),
        plotOutput("plot"),
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
      
      conditionalPanel(
        condition = "output.plot2",
        h3("Predictive gene count (corr. scores of magnitude > 1):")
      ),
      conditionalPanel(
        condition = "input.runpcl",
        plotOutput("plot2"),
        fluidRow(
          column(6,
                 uiOutput("slider_plot2")),
          column(4,
                 hr(),
                 textOutput("plot2_caption"))
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
        tabsetPanel(
          id = 'dataset',
          tabPanel('Pos. scores', DT::dataTableOutput('ptable'),
                   conditionalPanel(condition="output.ptable",downloadButton("ptable_dl","Download"))),
          tabPanel('Neg. scores', DT::dataTableOutput('ntable'),
                   conditionalPanel(condition="output.ntable",downloadButton("ntable_dl","Download"))),
          tabPanel('Pos. exp. values by age', DT::dataTableOutput('ppcl'),
                   conditionalPanel(condition="output.ppcl",downloadButton("ppcl_dl","Download"))),
          tabPanel('Neg. exp. values by age', DT::dataTableOutput('npcl'),
                   conditionalPanel(condition="output.npcl",downloadButton("npcl_dl","Download")))
        ),
        hr(),
        tabsetPanel(
          id = 'heatmap',
          tabPanel('Pos. gene by age heatmap', d3heatmapOutput("posheat")),
          tabPanel('Neg. gene by age heatmap', d3heatmapOutput("negheat"))
        )
      ),
      conditionalPanel(
        condition = "output.ptable",
            actionButton("rungt", "Run GO term analysis"),
            helpText("See the most enriched GO terms in your selected genes."),
          
        
        hr()
      ),
      conditionalPanel(
        condition = "input.rungt",
        h3("GO term analysis"),
        tabsetPanel(
          id = 'goterms',
          tabPanel('Pos. gene GO terms',DT::dataTableOutput('pos_goterms'),
                   conditionalPanel(condition="output.pos_goterms",downloadButton("pgo_dl","Download"))),
          tabPanel('Neg. gene GO terms',DT::dataTableOutput('neg_goterms'),
                   conditionalPanel(condition="output.neg_goterms",downloadButton("ngo_dl","Download")))  
        ),
        hr()
      )
    )
  )
))