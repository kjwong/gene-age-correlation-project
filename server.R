## Name: Kenny Wong
## Collaborators: Arjun Krishnan
## Murphy Lab


if (!require("topGO")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("topGO")
}
if (!require("org.Hs.eg.db")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite('org.Hs.eg.db')
}
if (!require("GO.db")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GO.db")
}
if (!require("GOstats")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GOstats")
}
if (!require("Rgraphviz")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("Rgraphviz")
} 
if (!require("S4Vectors")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("S4Vectors")
} 

if (!require("networkD3")) install.packages("networkD3")
if (!require("igraph")) install.packages('igraph')
if (!require("data.table")) install.packages('data.table')
if (!require("parallel")) install.packages('parallel')
if (!require("shiny")) install.packages('shiny')
if (!require("DT")) install.packages('DT')
if (!require("d3heatmap")) install.packages('d3heatmap')
if (!require("ggplot2")) install.packages('ggplot2')
if (!require("plotly")) install.packages('plotly')
if (!require("stats")) install.packages('stats')
if (!require("graph")) install.packages('graph')
if (!require("stats4")) install.packages('stats4')

library("data.table")
library(parallel)
library(shiny)
library(DT)
library(d3heatmap)
library(ggplot2)
library(org.Hs.eg.db)
library(topGO)
library(plotly)
library(stats)
library(Rgraphviz)
library(S4Vectors)
library(stats4)
library(graph)
library(networkD3)
library(Rgraphviz)

options(shiny.maxRequestSize = 1000*1024^2) # for large inputs
load("default.RData")

shinyServer(function(input, output) {
  
  # INPUT AND FILTERS
  
  # reading the entire pcl 
  gsm_pcl <- reactive({
    inFile <- input$file1
    withProgress(message = "Loading expression data",value=0.1,{
      if (is.null(inFile)) {
        if (input$default == "blood.f") gsm_pcl <- blood.f
        else if (input$default == "GSE58137") gsm_pcl <- GSE58137.entrez
      }
      else gsm_pcl <- data.frame(fread(inFile$datapath, header = input$header1,sep=input$sep1), row.names = 1)
      incProgress(0.9)
    })
    gsm_pcl
  })

  # reading-in sample annotations -INPUT
  # this has the ages and other parameters/filters
  gsm_input <- reactive ({
      inFile <- input$file2
      if (is.null(inFile)) {
        if (input$default == "blood.f") gsm_input <- sample_f
        else if (input$default == "GSE58137") gsm_input <- GSE58137.samples
        colnames(gsm_input) <- tolower(colnames(gsm_input))
      }
      else { 
        gsm_input <- read.table(inFile$datapath, header=input$header2, sep = input$sep2, row.names=1)
        colnames(gsm_input) <- tolower(colnames(gsm_input))
      }
    gsm_input$age <- as.integer(round(gsm_input$age,0))
    gsm_input
  })

  # dynamic filters
  output$filters <- renderUI({
    input$default
    df <- gsm_input()
    df <- df[, !(colnames(df) %in% c("age"))]
    numcol <- length(colnames(df)) 
    lapply(1:numcol,function(i){
      selectInput(inputId=paste0("filter",i),label=paste("Filter samples by",colnames(df)[i]),as.character(unique(df[,i])),multiple=TRUE)
    })
  })

  # list of samples filtered down
  st_samples <- eventReactive(input$run_samples,{
    df <- gsm_input()
    df <- df[!(is.na(df$age) | df$age==""), ]
    df <- df[, !(colnames(df) %in% c("age"))]
    numcol <- length(colnames(df)) 
    for (i in seq(1:numcol)) {
      filter <- input[[paste0("filter",i)]]
    
      if (length(filter) != 0) { 
        df <- df[which(df[,i] %in% filter),]
      }
    }
    rownames(df)
  })
  
  output$nosamp <- renderUI({
    gsm_age <- gsm_age()
    if (nrow(gsm_age) == 0) (helpText("There are 0 samples that fit your specified filters."))
    else return()
  })
  
  # list of gene names
  gene_sym <- reactive({
    read.delim("human_gene-info_ncbi.txt", header=T, row.names=1, sep="\t")
  })

  # list of sample ages
  gsm_age <- eventReactive(input$run_samples,{
    gsm_input <- gsm_input()
    st_samples <- st_samples()
    
    if (length(st_samples) == 0) return(data.frame())
    gsm_age <- data.frame(gsm_input["age"])
    gsm_age <- data.frame(gsm_age[which(rownames(gsm_age) %in% st_samples),],rownames=FALSE)
    rownames(gsm_age) <- intersect(rownames(gsm_input),st_samples)
    gsm_age
  })

  # sample age range
  lwr <- reactive({min(gsm_age()[,1])})
  upr <- reactive({max(gsm_age()[,1])})
  arng <- reactive({
    arng <- lwr():upr()
    count = 1
    for (i in lwr():upr()) {
      if (!(i %in% st_gsm_age()[,1])) {
        arng <- arng[-count]
      }
      else count = count + 1
    }
    arng
  })
  # input age range
  inlwr <- reactive({input$age_range[1]})
  inupr <- reactive({input$age_range[2]})
  inarng <- reactive({
    arng <- inlwr():inupr()
    count = 1
    for (i in inlwr():inupr()) {
      if (!(i %in% st_gsm_age()[,1])) {
        arng <- arng[-count]
      }
      else count = count + 1
    }
    arng
  })
  # pos heat map age range
  heatposlwr <- reactive({input$heatposage[1]})
  heatposupr <- reactive({input$heatposage[2]})
  heatposarng <- reactive({
    arng <- heatposlwr():heatposupr()
    count = 1
    for (i in heatposlwr():heatposupr()) {
      if (!(i %in% st_gsm_age()[,1])) {
        arng <- arng[-count]
      }
      else count = count + 1
    }
    arng
  })
  
  # neg heat map age range
  heatneglwr <- reactive({input$heatnegage[1]})
  heatnegupr <- reactive({input$heatnegage[2]})
  heatnegarng <- reactive({
    arng <- heatneglwr():heatnegupr()
    count = 1
    for (i in heatneglwr():heatnegupr()) {
      if (!(i %in% st_gsm_age()[,1])) {
        arng <- arng[-count]
      }
      else count = count + 1
    }
    arng
  })

  # captions for plot1
  plotcap <- reactive({
    gsm_age <- gsm_age()
    paste("No. of samples selected:\n",dim(gsm_age[which((gsm_age[,1]>=inlwr()) & (gsm_age[,1]<=inupr())),])[1])
  })
  output$plot1_caption <- renderUI({
    h5(plotcap())
  })
  plotcap1 <- reactive({
    gsm_age <- gsm_age()
    num <- dim(gsm_age[which((gsm_age[,1]>=inlwr()) & (gsm_age[,1]<=inupr())),])[1]
    paste0("Percentage of samples selected: ",round(num/dim(gsm_age)[1] * 100,3),"%")
  })
  output$plot1_caption2 <- renderUI({
    h5(plotcap1())
  })
  
  # plot 1 slider
  output$slider_plot <- renderUI({
    gsm_age <- gsm_age()
    quad <- (max(gsm_age[,1]) - min(gsm_age[,1])) / 4
    if (round(min(gsm_age[,1]) + quad,-1) >= min(gsm_age[,1])) min <- round(min(gsm_age[,1]) + quad,-1)
    else min <- min(gsm_age[,1]) + quad
    if (round(max(gsm_age[,1] - 2 * quad),-1) <= max(gsm_age[,1])) max <- round(max(gsm_age[,1] - 2 * quad),-1)
    else max <- max(gsm_age[,1] - 2 * quad)
    sliderInput("age_range", label = "Choose an age range:", min = min(gsm_age[,1]), 
                max = max(gsm_age[,1]), value = c(min,max))
                
  })
  plot1 <- reactive({
    .e <- environment()
    gsm_age <- gsm_age()
    colnames(gsm_age) <- c("age")
    pdf(NULL)
    p <- ggplot(data=gsm_age,aes(x=age),environment=.e) +
      geom_histogram(binwidth=2,fill="blue",col="gray",alpha=0.4) +
      labs(x="age",y="") +
      theme(axis.text=element_text(size=9),
            axis.title=element_text(size=9))
    p <- ggplotly(p)
  })
  # plot 1
  output$plot <- renderPlotly({
    p <- plot1()

#     p <- p %>%
#       add_trace(x = c(inlwr(), inlwr()), y= c(0, 100000), mode = "lines", name = "Lower", text = "", line=list(
#         dash = "dashdot", color = "gray",alpha=0.5)) %>%
#       add_trace(x = c(inupr(),inupr()), y= c(0, 100000), mode = "lines", name = "Upper", text="",line=list(
#         dash = "dashdot", color = "gray",alpha=0.5))

    layout(p, hovermode = "closest", showlegend=FALSE)
    
  })

  
  
#
# running gene-age correlation
#
  


  # intersecting col names with sample ids
  st_make <- eventReactive(input$runpcl,{
    gsm_age <- gsm_age()
    gsm_pcl <- gsm_pcl()
    intersect(rownames(gsm_age)[which((gsm_age[,1]>=lwr()) & (gsm_age[,1]<=upr()))],
              colnames(gsm_pcl))
  })

  output$nosamp2 <- renderUI({
    if (length(st_make()) != 0) return()
    else (helpText("There are 0 samples in your expression data that fit your specified filters."))
  })

  # creating pcl for samples
  st_gsm_pcl <- reactive({
    gsm_pcl <- gsm_pcl()
    st_make <- st_make()
    gsm_pcl[,st_make]
  })

  # gsm age for samples
  st_gsm_age <- reactive({
    gsm_age <- gsm_age()
    st_make <- st_make()
    gsm_age[st_make,]
  })

  # creating pcl of genes per age with median of gene-expression values per gene
  st_age_pcl <- function(arng){
    st_gsm_age <- st_gsm_age()
    st_gsm_pcl <- st_gsm_pcl()
    st_age_pcl <- array(NaN, c(nrow(st_gsm_pcl()), length(arng)))
    for(j in 1:length(arng)) {
      age_gsm <- rownames(st_gsm_age)[which(st_gsm_age[,1]==arng[j])]
      if(length(age_gsm)==1) {
        st_age_pcl[,j] <- st_gsm_pcl[,age_gsm]
      } else {
        st_age_pcl[,j] <- apply(st_gsm_pcl[,age_gsm], 1, median)
      }
    }
    rownames(st_age_pcl) <- rownames(st_gsm_pcl)
    st_age_pcl
  }

  # scaling genes (rows) to mean 0 and variance 1
  st_age_pcl_rowz <- function(arng){
    withProgress(message="Processing heat map", detail="Creating expression matrix of genes by age",value = 0.1, {
      st_age_pcl <- st_age_pcl(arng)
      st_gsm_pcl <- st_gsm_pcl()
      incProgress(0.6,detail="Scaling rows to mean 0 and variance 1")
      idx <- apply(st_age_pcl,1,var)
      st_age_pcl_rowz <- t(apply(st_age_pcl, 1, scale))
      rownames(st_age_pcl_rowz) <- rownames(st_gsm_pcl)
      incProgress(0.3,detail="Removing null values")
      st_age_pcl_rowz <- st_age_pcl_rowz[-which(idx==0),]
      colnames(st_age_pcl_rowz) <- arng
      st_age_pcl_rowz
    })
  }

  # boot_rho is the output of the bootstrap runs
  # It contains <nboot> rows and <#genes> columns
  boot_rho <- reactive({
    clust <- makeCluster(10) # Initiate cluster
    arng <- inarng()
    nboot <- input$numruns
    st_gsm_age <- st_gsm_age()
    st_gsm_pcl <- st_gsm_pcl()
    boot_rho <- array(NaN, c(nboot, nrow(st_gsm_pcl)))
    withProgress(message = 'Calculating correlation scores', value = 0, {
      for(n in 1:nboot) {
        incProgress(1/nboot,detail=paste("Run",n))
        cat(n, "...\n")
        bag_gsm <- {}
        # for each age...
        for(j in 1:length(arng)) {
          age_gsm <- rownames(st_gsm_age)[which(st_gsm_age[,1]==arng[j])] # GSM names of age
          
          if(length(age_gsm)==1) {
            bag_gsm <- c(bag_gsm, age_gsm) # if only 1 sample for the age, just add it to the bag
          } else {
            sel_gsm <- sample(age_gsm, replace=T) # sel_gsm will differ every run
            bag_gsm <- c(bag_gsm, unique(sel_gsm))
          }
        }
        boot_pcl <- st_gsm_pcl[,bag_gsm] # exp. values
        boot_age <- st_gsm_age[bag_gsm,][,1] # ages
        clusterExport(clust,varlist = c("boot_pcl","boot_age"),envir=environment()) # cluster environment
        rho <- parApply(clust, boot_pcl, 1, cor, y=boot_age, method=input$corr) 
        boot_rho[n,] <- rho
      }
    })
    stopCluster(clust)
    colnames(boot_rho) <- rownames(st_gsm_pcl) # gene names
    boot_rho
  })

  # processing rho into comparable correlation scores (fisherz)
  bxs_boot_fisherz <- reactive({
      st_gsm_pcl <- st_gsm_pcl()
      boot_fisherz <- t(apply(boot_rho(), 1, function(x) { scale(0.5*log((1+x)/(1-x))) })) # Fisher z-transformation
      bxs_boot_fisherz <- apply(boot_fisherz, 2, quantile, probs=c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0), na.rm=T)
      bxs_boot_fisherz <- data.frame(bxs_boot_fisherz)
      colnames(bxs_boot_fisherz) <- rownames(st_gsm_pcl)
      bxs_boot_fisherz
  })


#
# predictive genes
#


  # magnitudes of scores of predictive genes
  abs_scores <- eventReactive(input$runpcl,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    all_predg <- all_predg()
    data.frame(abs(bxs_boot_fisherz[4,all_predg]))
  })

  # scores of predictive genes
  raw_scores <- eventReactive(input$runpcl,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    all_predg <- all_predg()
    data.frame(bxs_boot_fisherz[4,all_predg])
  })
  
  
  # number of genes
  num_genes <- reactive({
    abs_scores <- abs_scores()
    df<-abs_scores[which(as.numeric(abs_scores[1,])>input$score_mag)]
    ncol(df)
  })
  
  # histogram plot2 of magnitudes
  output$plot2 <- renderPlotly({
    .e <- environment()
    df <- cbind(data.frame(all_predg()), t(raw_scores()))
    score <- df[,2]
    pdf(NULL)
    p <- ggplot(data=df, aes(x=score), environment = .e) + 
      geom_histogram(binwidth=0.1,fill="red",col="gray",alpha=0.4) +
      labs(x=paste("correlation coefficient"), y="") +
      theme(axis.text=element_text(size=9),
            axis.title=element_text(size=9))
    p <- ggplotly(p)
    layout(p, hovermode="closest",showlegend=FALSE)
  })
  
  # slider for plot2
  output$slider_plot2 <- renderUI({
    lwrbound = signif(min(abs_scores()),3)
    uprbound = signif(max(abs_scores()),3)
    sliderInput("score_mag", label = "Choose a score magnitude cutoff:", ticks=FALSE,min = lwrbound, 
                max = uprbound,value= lwrbound)
  })
  # captions for plot2
  plotcap2 <- reactive({
    paste("No. of genes selected:\n",num_genes())
  })
  plotcap22 <- reactive({
    paste0("Percentage of genes selected: ",round(num_genes()/dim(gsm_pcl())[1] * 100,3),"%")
  })
  output$plot2_caption <- renderUI({
    h5(plotcap2())
  })
  output$plot2_caption2 <- renderUI({
    h5(plotcap22())
  })
  

  # positive predictors
  pos_predg <- reactive({
    bxs_boot_fisherz <- bxs_boot_fisherz()
    pos_predg <- colnames(bxs_boot_fisherz)[which(bxs_boot_fisherz[4,]>=input$score_mag)]
    pos_predg
  })

  # negative predictors
  neg_predg <- reactive({
    bxs_boot_fisherz <- bxs_boot_fisherz()
    neg_predg <- colnames(bxs_boot_fisherz)[which(bxs_boot_fisherz[4,]<0 & abs(bxs_boot_fisherz[4,])>=input$score_mag)]
  })

  # all predictive genes
  all_predg <- eventReactive(input$runpcl, {
    bxs_boot_fisherz <- bxs_boot_fisherz()
    withProgress(message="Identifying predictive genes", value = 0.1, {
      all_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]>=quantile(bxs_boot_fisherz[4,],0.95,na.rm=TRUE)) | 
                                                      (bxs_boot_fisherz[4,]<0 & bxs_boot_fisherz[4,]<=quantile(bxs_boot_fisherz[4,],0.05,na.rm=TRUE)))]
      incProgress(0.9)
    })
    all_predg
  })

  # positive pcl table with exp values
  pos_pcl <- eventReactive(input$runposheat,{
      bxs_boot_fisherz <- bxs_boot_fisherz()
      st_age_pcl_rowz <- st_age_pcl_rowz(heatposarng())
      pos_predg <- pos_predg()
      pos_predg <- intersect(pos_predg, rownames(st_age_pcl_rowz))
      pos_predg_pcl_rowz <- cbind(gene_sym()[pos_predg,],
                                  t(bxs_boot_fisherz[4,pos_predg]), 
                                  apply(st_age_pcl_rowz[pos_predg,], 2, sprintf, fmt="%.6f"))
      colnames(pos_predg_pcl_rowz) <- c("symbol", "name","score", heatposarng())
    pos_predg_pcl_rowz[order(-pos_predg_pcl_rowz$score,pos_predg_pcl_rowz$name),]
  })
  # positive score table w/o exp values
  pos_score <- eventReactive(input$tablepcl,{
    withProgress(message="Generating positive score table",value=0.1,{
      bxs_boot_fisherz <- bxs_boot_fisherz()
      pos_predg <- pos_predg()
      pos_predg <- intersect(pos_predg, rownames(st_gsm_pcl()))
      pos_predg_pcl_rowz <- cbind(gene_sym()[pos_predg,],
                                  round(t(bxs_boot_fisherz[4,pos_predg]),3))
      colnames(pos_predg_pcl_rowz) <- c("symbol", "name","score")
      incProgress(0.9)
      res <- pos_predg_pcl_rowz[order(-pos_predg_pcl_rowz$score,pos_predg_pcl_rowz$name),]
    })
  })
  # negative pcl table with exp values
  neg_pcl <- eventReactive(input$runnegheat,{
      bxs_boot_fisherz <- bxs_boot_fisherz()
      st_age_pcl_rowz <- st_age_pcl_rowz(heatnegarng())
      neg_predg <- neg_predg()
      neg_predg <- intersect(neg_predg, rownames(st_age_pcl_rowz))
      neg_predg_pcl_rowz <- cbind(gene_sym()[neg_predg,],
                                  t(bxs_boot_fisherz[4,neg_predg]), 
                                  apply(st_age_pcl_rowz[neg_predg,], 2, sprintf, fmt="%.6f"))
      colnames(neg_predg_pcl_rowz) <- c("symbol", "name","score", heatnegarng())
    neg_predg_pcl_rowz[order(neg_predg_pcl_rowz$score,neg_predg_pcl_rowz$name),]
  })
  # negative score table w/o exp values
  neg_score <- eventReactive(input$tablepcl,{
    withProgress(message="Generating negative score table",value=0.1,{
      bxs_boot_fisherz <- bxs_boot_fisherz()
      neg_predg <- neg_predg()
      neg_predg <- intersect(neg_predg, rownames(st_gsm_pcl()))
      neg_predg_pcl_rowz <- cbind(gene_sym()[neg_predg,],
                                  round(t(bxs_boot_fisherz[4,neg_predg]),3))
      colnames(neg_predg_pcl_rowz) <- c("symbol", "name","score")
      incProgress(0.9)
      res <- neg_predg_pcl_rowz[order(neg_predg_pcl_rowz$score,neg_predg_pcl_rowz$name),]
    })
  })

  


#
# outputs
#

  # heat maps
  posheatmap <- eventReactive(input$runposheat,{
    pos_p <- pos_pcl()
    pos_p <- pos_p[!(is.na(pos_p[,1]) | pos_p[,1]==""), ]
    rownames(pos_p) <- pos_p[,1]
    pos_p <- pos_p[,-c(1,2,3)] # remove columns other than exp. values
    # convert factor to numeric data frame
    indx <- sapply(pos_p, is.factor)
    pos_p[indx] <- lapply(pos_p[indx], function(x) as.numeric(as.character(x)))
    length <- input$heatposnum                       # if user inputs number of genes
    if (dim(pos_p)[1] < length) length <- dim(pos_p)[1]           
    usergenes <- toupper(unlist(strsplit(input$userposgenes, " "))) # if user inputs gene names
    if (length(usergenes) == 0) pos_op <- pos_p[1:length,]
    else pos_op <- pos_p[usergenes,]
    d3heatmap(pos_op,Colv = FALSE,xaxis_font_size=7,yaxis_font_size=7)
  })
 
  negheatmap <- eventReactive(input$runnegheat,{
    neg_p <- neg_pcl()
    neg_p <- neg_p[!(is.na(neg_p[,1]) | neg_p[,1]==""), ]
    rownames(neg_p) <- neg_p[,1]
    neg_p <- neg_p[,-c(1,2,3)]
    # convert factor to numeric data frame
    indx <- sapply(neg_p, is.factor)
    neg_p[indx] <- lapply(neg_p[indx], function(x) as.numeric(as.character(x)))
    length <- input$heatnegnum
    if (dim(neg_p)[1] < length) length <- dim(neg_p)[1] 
    usergenes <- toupper(unlist(strsplit(input$userneggenes, " ")))
    if (length(usergenes) == 0) neg_op <- neg_p[1:length,]
    else neg_op <- neg_p[usergenes,]
    d3heatmap(neg_op,Colv = FALSE,xaxis_font_size=7,yaxis_font_size=7)
  })

  # predictive gene tables
  output$ptable <- DT::renderDataTable({
    DT::datatable(pos_score(), rownames=FALSE)
  })
  
  output$ntable <- DT::renderDataTable({
    DT::datatable(neg_score(), rownames=FALSE)
  })
  
  # heat map age sliders
  output$heatposage <- renderUI({
    sliderInput("heatposage",label="Choose an age range:",value=c(input$age_range[1],input$age_range[2]),min=min(arng()),max=max(arng()))
  })
  output$heatnegage <- renderUI({
    sliderInput("heatnegage",label="Choose an age range:",value=c(input$age_range[1],input$age_range[2]),min=min(arng()),max=max(arng()))
  })
  
  # rendering heat maps
  output$posheat <- renderD3heatmap({
    posheatmap()
  })
  
  output$negheat <- renderD3heatmap({
    negheatmap()
  })

  # GO term p values
  pos_go_results <- reactive({
    sampleGOdata <- pos_go()
    result <- runTest(sampleGOdata, algorithm = "weight01", statistic = input$posstat)
  })

  neg_go_results <- reactive({
    sampleGOdata <- neg_go()
    result <- runTest(sampleGOdata, algorithm = "weight01", statistic = input$negstat)
  })

  # GO term tables
  output$pos_goterms <- DT::renderDataTable({
    allRes <- pos_go_table()
    DT::datatable(allRes, rownames=FALSE)
  })

  output$neg_goterms <- DT::renderDataTable({
    allRes <- neg_go_table()
    DT::datatable(allRes, rownames=FALSE)
  })

  # mapping for GO term hierarchy
  posmap <- eventReactive(input$runposgt,{
    if (input$postype == "BP") map <- GOBPPARENTS
    else if (input$postype == "MF") map <- GOMFPARENTS
    else if (input$postype == "CC") map <- GOCCPARENTS
    map
  })
  negmap <- eventReactive(input$runneggt,{
    if (input$negtype == "BP") map <- GOBPPARENTS
    else if (input$negtype == "MF") map <- GOMFPARENTS
    else if (input$negtype == "CC") map <- GOCCPARENTS
    map
  })
  # creating positive GO term data and tables
  pos_go <- eventReactive(input$runposgt,{
    all_predg <- pos_predg()
    withProgress(message = "Generating positive GO terms\n", detail = "Compiling tables", value = 0.05, {
      gsm_pcl <- gsm_pcl()
      incProgress(0.05, detail = "Compiling tables")
      ann <- annFUN.org(input$postype,feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
      incProgress(0.1, detail = "Building GO DAG topology and annotations")
      geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg))
      names(geneList) <- rownames(gsm_pcl)
      selector <- function(x) {return (x==1)}
      sampleGOdata <- new("topGOdata",ontology=input$postype,allGenes=geneList,
                          geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
      incProgress(0.8)
    })
    sampleGOdata
  })
  
  pos_go_table <- reactive({
    sampleGOdata <- pos_go()
    withProgress(message = "Generating positive GO terms\n", detail = paste("Running", toupper(input$posstat),"test"), value = 0.1,{
      result <- pos_go_results()
      incProgress(0.7, detail = "Aggregating results")
      allRes <- GenTable(sampleGOdata, pValue = result,
                         orderBy = "pValue",topNodes=1000)
      allRes[,"Fold Enrichment"] <- round(allRes[,"Significant"] / allRes[,"Expected"],1)
      incProgress(0.2)
      allRes <- allRes[c("GO.ID","Term","Annotated","Significant","Fold Enrichment","pValue")]
      p <- allRes[,"pValue"]
      allRes[,"FDR"] <- round(p.adjust(p, method = "BH", n = length(p)),5)
      colnames(allRes) <- c("GO.ID", "GO Term","No. of genes","Overlap","Fold Enrichment","P-value","FDR")
      allRes <- allRes[allRes$P-value < 0.05,]
      allRes <- allRes[order(allRes$FDR),]
    })
    allRes
  })
  
  # creating negative GO term data and tables
  neg_go <- eventReactive(input$runneggt,{
    all_predg <- neg_predg()
    withProgress(message = "Generating negative GO terms\n", detail = "Compiling tables", value = 0.05, {
      gsm_pcl <- gsm_pcl()
      incProgress(0.05, detail = "Compiling tables")
      ann <- annFUN.org(input$negtype,feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
      incProgress(0.1, detail = "Building GO DAG topology and annotations")
      geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg)) #0 if in list and 1 if not
      names(geneList) <- rownames(gsm_pcl) 
      selector <- function(x) {return (x==1)} 
      sampleGOdata <- new("topGOdata",ontology=input$negtype,allGenes=geneList,
                          geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
      incProgress(0.8)
    })
    sampleGOdata
  })
  
  neg_go_table <- reactive({
    sampleGOdata <- neg_go()
    withProgress(message = "Generating negative GO terms\n",detail = paste("Running", toupper(input$negstat),"test"),value=0.1,{
      result <- neg_go_results()
      incProgress(0.7, detail = "Aggregating results")
      allRes <- GenTable(sampleGOdata, pValue = result,
                         orderBy = "pValue",topNodes=1000)
      allRes[,"Fold Enrichment"] <- round(allRes[,"Significant"] / allRes[,"Expected"],1)
      incProgress(0.2)
      allRes <- allRes[c("GO.ID","Term","Annotated","Significant","Expected","Fold Enrichment","pValue")]
      allRes <- allRes[,-c(5)]
      p <- allRes[,"pValue"]
      allRes[,"FDR"] <- round(p.adjust(p, method = "fdr", n = length(p)),5)
      colnames(allRes) <- c("GO.ID", "GO Term","No. of genes","Overlap","Fold Enrichment","P-value","FDR")
      allRes <- allRes[allRes$P-value < 0.05,]
      allRes <- allRes[order(allRes$FDR),]
    })
    allRes
  })

  numposGO <- reactive({
    dim(pos_go_table())[1] 
  })
  numpos <- reactive({
    input$posnodes
  })
  output$posn <- renderUI({
    sliderInput("posnodes",label="No. of significant nodes:",value=min(numposGO(),5),min=1,max=min(numposGO(),15))
  })
  # generate induced subgraphs
  output$pos_go_graph <- renderForceNetwork({
    # table of GO IDs, terms, pvalues
    allRes <- pos_go_table()
    num <- min(numpos(), numposGO())
    # make GO graph
    gn <- GOGraph(allRes$GO.ID[1:num],posmap())
    ig <- igraph.from.graphNEL(gn) # reformat graph
    gd <- get.data.frame(ig, what = "edges") # graph edges
    nod = data.frame(unique(c(unique(gd[,1]),unique(gd[,2])))) # graph nodes
    rownames(nod) <- nod[,1]
    # group nodes based on significant or normal nodes
    nod[,2] = 1 # normal nodes
    nod[which(nod[,1] %in% allRes$GO.ID[1:num]),2] = 2  # sig nodes
    colnames(nod) <- c("Name","Group")
    # transform edges from GO IDs to ranks
    gd[,1] <- match(gd[,1],nod[,1]) - 1
    gd[,2] <- match(gd[,2],nod[,1]) - 1
    # transform nodes from GO IDs to GO terms
    nod[,1] = apply(data.frame(nod[,1]),1,function(x) Term(x))
    pdf(NULL)
    forceNetwork(Links=gd,Nodes=nod,Source="from",Target="to",Value="weight",NodeID="Name",Group="Group",
                 colourScale=JS("d3.scale.category10()"),zoom=TRUE,opacity=0.7,opacityNoHover = TRUE)
  })
  numnegGO <- reactive({
    dim(neg_go_table())[1]
  })
  numneg <- reactive({
    input$negnodes
  })
  output$negn <- renderUI({
    sliderInput("negnodes",label="No. of significant nodes:",value=min(numposGO(),5),min=1,max=min(numnegGO(),15))
  })
  output$neg_go_graph <- renderForceNetwork({
    allRes <- neg_go_table()
    num <- min(numneg(),numnegGO())
    # make GO graph
    gn <- GOGraph(allRes$GO.ID[1:num],negmap())
    ig <- igraph.from.graphNEL(gn) # reformat graph
    gd <- get.data.frame(ig, what = "edges") # graph edges
    nod = data.frame(unique(c(unique(gd[,1]),unique(gd[,2])))) # graph nodes
    rownames(nod) <- nod[,1]
    # group nodes based on significant or normal nodes
    nod[,2] = 1 # normal nodes
    nod[which(nod[,1] %in% allRes$GO.ID[1:num]),2] = 2  # sig nodes
    colnames(nod) <- c("Name","Group")
    # transform edges from GO IDs to ranks
    gd[,1] <- match(gd[,1],nod[,1]) - 1
    gd[,2] <- match(gd[,2],nod[,1]) - 1
    # transform nodes from GO IDs to GO terms
    nod[,1] = apply(data.frame(nod[,1]),1,function(x) Term(x))
    pdf(NULL)
    forceNetwork(Links=gd,Nodes=nod,Source="from",Target="to",Value="weight",NodeID="Name",Group="Group",
                 colourScale=JS("d3.scale.category10()"),zoom=TRUE,opacity=0.7,opacityNoHover = TRUE)
  })

# download options for tables

  output$ptable_dl <- downloadHandler(
    filename = function(){paste0(inlwr(),"-",inupr(),"_pos_genes.csv")},
    content = function(file){write.csv(pos_score(),file)}
  )

  output$ntable_dl <- downloadHandler(
    filename = function(){paste0(inlwr(),"-",inupr(),"_neg_genes.csv")},
    content = function(file){write.csv(neg_score(),file)}
  )

  output$pgo_dl <- downloadHandler(
    filename = function(){paste0(inlwr(),"-",inupr(),"_pos_goterms.csv")},
    content = function(file){write.csv(pos_go_table(),file)}
  )

  output$ngo_dl <- downloadHandler(
    filename = function(){paste0(inlwr(),"-",inupr(),"_neg_goterms.csv")},
    content = function(file){write.csv(neg_go_table(),file)}
  )
})
  
  
