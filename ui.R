library(shiny)
library(DT)
library(d3heatmap)
library(ggplot2)

shinyUI(fluidPage(
  
  titlePanel("Web interface for gene-age correlation in GSM samples"),
  hr(),
#   fluidRow(
#     column(3,
#            helpText("Select a GSM expression file. (.pcl)"),
#            fileInput('file1', 'Choose file to upload',
#                      accept = c('.pcl')
#            ),
#            tags$hr(),
#            checkboxInput('header1', 'Header', TRUE),
#            radioButtons('sep1', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), '\t'),
#            actionButton("upload1", "Enter"),
#            tags$hr()
#     ),
#     column(3,
#            helpText("Select a file with gene IDs, symbols, and names. (.txt)"),
#            fileInput('file2', 'Choose file to upload',
#                      accept = c('.txt')
#            ),
#            tags$hr(),
#            checkboxInput('header2', 'Header', TRUE),
#            radioButtons('sep2', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), '\t'),
#            actionButton("upload2", "Enter"),
#            tags$hr()
#     ),
#     column(3,
#            helpText("Select a file with sample age annotations. (.txt)"),
#            fileInput('file3', 'Choose file to upload',
#                      accept = c('.txt')
#            ),
#            tags$hr(),
#            checkboxInput('header3', 'Header', TRUE),    
#            radioButtons('sep3', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), '\t'),
#            actionButton("upload3", "Enter"),
#            tags$hr()
#     ), 
#     column(3,
#            helpText("Select a file with your selected samples. (.txt)"),
#            fileInput('file4', 'Choose file to upload',
#                      accept = c('.txt')
#            ),
#            tags$hr(),
#            checkboxInput('header4', 'Header', TRUE),    
#            radioButtons('sep4', 'Separator', c(Comma=',', Semicolon=';',Tab='\t'), '\t'),
#            actionButton("upload4", "Enter"),
#            tags$hr()
#     )       
#   ),
  sidebarLayout(
    sidebarPanel(
      # checkbox button for sex
      checkboxGroupInput("sex", label = h4("Sex:"), 
                         choices = list("Male" = 1, "Female" = 2), selected = (2)),
      checkboxGroupInput("tissue", label = h4("Tissue:"), choices = list("Blood" = 1, "Other" = 2), selected = 1),
      checkboxGroupInput("status", label = h4("Status:"), choices = list("Healthy" = 1, "Other" = 2), selected = 1),
      hr(),
      actionButton("run_filter", "Run filter"),
      # change later to min and max values filtered by tissue type and sex
      conditionalPanel('input.dataset === "Pos. scores"'),
      conditionalPanel('input.dataset === "Neg. scores"'),
      conditionalPanel('input.dataset === "Pos. exp. values by age"'),
      conditionalPanel('input.dataset === "Neg. exp. values by age"')
    ),
    
    mainPanel(
      textOutput("test"),
      h3("GSM sample count by age:"),
      plotOutput("plot"),
      fluidRow(
        column(6,
          uiOutput("slider_plot")),
        column(4,
          textOutput("plot_caption"))
      ),
      helpText("Use the slider to select samples of an age range."),
      actionButton("runpcl","Run gene-age correlation"),
      hr(),
      conditionalPanel(
        condition = "output.plot2",
        h3("Predictive gene count (corr. scores of magnitude > 1):")
      ),
      plotOutput("plot2"),
      fluidRow(
        column(6,
          uiOutput("slider_plot2")),
        column(4,
          textOutput("plot2_caption"))
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
        fluidRow(
          column(4,
            actionButton("rungt", "Run GO term analysis"),
            helpText("See the most enriched GO terms in your selected genes.")
          )
        ),
        hr()
      ),
      conditionalPanel(
        condition = "input.rungt",
        h3("GO term analysis"),
        tabsetPanel(
          id = 'goterms',
          tabPanel('Pos. gene GO terms',DT::dataTableOutput('pos_goterms')),
          tabPanel('Neg. gene GO terms',DT::dataTableOutput('neg_goterms'))  
        ),
        hr()
      ),
      helpText("Content by Murphy Lab at Princeton University")
    )
  )
))