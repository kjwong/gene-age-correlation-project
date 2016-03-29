library(shiny)
library(DT)
library(d3heatmap)
library(ggplot2)

shinyUI(fluidPage(
  
  titlePanel("Gene-age correlation in GSM samples"),
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
      h3("# of GSM samples by age:"),
      plotOutput("plot"),
      fluidRow(
        column(6,
          uiOutput("slider_plot")),
        column(4,
          textOutput("plot_caption"))
      ),
      actionButton("runpcl","Run gene-age correlation"),
      hr(),
      h3("# of 'predictive' genes (scores of magnitude > 1):"),
      plotOutput("plot2"),
      fluidRow(
        column(6,
          uiOutput("slider_plot2")),
        column(4,
          textOutput("plot2_caption"))
      ),
      fluidRow(
        column(3,
          actionButton("tablepcl","See selected genes")),
        column(3,
          actionButton("rungt", "See GO term enrichment"))
      ),
      helpText("Note: GO term enrichment may take a while."),
      hr(),
      uiOutput("tablescap"),
      tabsetPanel(
        id = 'dataset',
        tabPanel('Pos. scores', DT::dataTableOutput('ptable')),
        tabPanel('Neg. scores', DT::dataTableOutput('ntable')),
        tabPanel('Pos. exp. values by age', DT::dataTableOutput('ppcl')),
        tabPanel('Neg. exp. values by age', DT::dataTableOutput('npcl'))
      ),
      hr(),
      tabsetPanel(
        id = 'heatmap',
        tabPanel('Pos. gene-age heatmap', d3heatmapOutput("posheat")),
        tabPanel('Neg. gene-age heatmap', d3heatmapOutput("negheat"))
      ),
      hr(),
      tabsetPanel(
        id = 'goterms',
        tabPanel('Pos. gene GO terms',DT::dataTableOutput('pos_goterms')),
        tabPanel('Neg. gene GO terms',DT::dataTableOutput('neg_goterms'))  
      )
    )
  )
))