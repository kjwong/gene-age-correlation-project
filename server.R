library("data.table")
library(parallel)
library(shiny)
library(DT)
library(d3heatmap)
library(ggplot2)
library(org.Hs.eg.db)
library(topGO)
options(shiny.maxRequestSize = 1000*1024^2)
#numcores <- detectCores() - 1 # Find no. cores

shinyServer(function(input, output) {
  
  # INPUT

  # reading the entire pcl 
  gsm_pcl <- reactive({
    withProgress(message = "Reading expression data",value=0.1,{
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)   
      gsm_pcl <- data.frame(fread(inFile$datapath, header = input$header1,sep=input$sep1), row.names = 1)
      incProgress(0.9)
    })
    gsm_pcl
  })

  # reading-in sample annotations -INPUT
  # this has the ages and other parameters/filters
  gsm_input <- eventReactive(input$upload2, {
      inFile <- input$file2
      if (is.null(inFile))
        return(NULL)   
      gsm_input <- read.table(inFile$datapath, header=input$header2, sep = input$sep2, row.names=1)
      colnames(gsm_input) <- tolower(colnames(gsm_input))
      
    gsm_input
  })

  # dynamic filters
  output$filters <- renderUI({
    df <- gsm_input()
    df <- df[, !(colnames(df) %in% c("age"))]
    numcol <- length(colnames(df)) 
    lapply(1:numcol,function(i){
      selectInput(inputId=paste0("filter",i),label=paste("Filter samples by",colnames(df)[i]),as.character(unique(df[,i])),multiple=TRUE)
    })
  })

  # list of samples filtered down
  st_samples <- reactive({
    df <- gsm_input()
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
  # list of gene names
  gene_sym <- reactive({
    read.delim("human_gene-info_ncbi.txt", header=T, row.names=1, sep="\t")
  })

  # list of sample ages
  gsm_age <- eventReactive(input$run_samples,{
    gsm_input <- gsm_input()
    st_samples <- st_samples()
    
    gsm_age <- data.frame(gsm_input["age"])
    gsm_age <- data.frame(gsm_age[which(rownames(gsm_age) %in% st_samples),],rownames=FALSE)    

    rownames(gsm_age) <- intersect(rownames(gsm_input),st_samples)
    gsm_age
  })
  
  # read button
  output$upload2 <- renderUI({
    if (is.null(input$file2)) return()
    actionButton("upload2", "Read in file")
  })
  
  # age range
  lwr <- reactive({input$age_range[1]})
  upr <- reactive({input$age_range[2]})
  arng <- reactive({lwr():upr()})
  
  # plot 1 caption
  output$plot_caption <- renderText({
    gsm_age <- gsm_age()
    paste("# of samples selected:\n",dim(gsm_age[which((gsm_age[,1]>=lwr()) & (gsm_age[,1]<=upr())),])[1])
  })
  
  # plot 1 slider
  output$slider_plot <- renderUI({
    gsm_age <- gsm_age()
    quad <- (max(gsm_age[,1]) - min(gsm_age[,1])) / 4
    sliderInput("age_range", label = h6("Choose an age range:"), min = min(gsm_age[,1]), 
                max = max(gsm_age[,1]), 
                value = c(min(gsm_age[,1]) + quad, max(gsm_age[,1] - 2 * quad)))
                
  })
  
  # plot 1
  output$plot <- renderPlot({
    .e <- environment()
    gsm_age <- gsm_age()
    p <- ggplot(data=gsm_age, aes(gsm_age[,1]), environment = .e) + 
      geom_histogram(binwidth=2,fill=I("gray23"),col=I("gray")) + 
      labs(x="age",y="# of GSM samples")
    p + geom_vline(xintercept=as.numeric(lwr()), colour="red",size=0.6, linetype="twodash") +
      geom_vline(xintercept=as.numeric(upr()), colour="red",size=0.6, linetype="twodash") 
  })

  
  
#
#
#
  


  # intersecting col names with sample ids
  st_make <- function(){
    gsm_age <- gsm_age()
    gsm_pcl <- gsm_pcl()
    intersect(rownames(gsm_age)[which((gsm_age[,1]>=lwr()) & (gsm_age[,1]<=upr()))],
              colnames(gsm_pcl))
  }

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
  st_age_pcl <- reactive({
    st_gsm_age <- st_gsm_age()
    st_gsm_pcl <- st_gsm_pcl()
    arng <- arng()
    st_age_pcl <- array(NaN, c(nrow(st_gsm_pcl()), length(lwr():upr())))
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
  })

  # scaling genes (rows) to mean 0 and variance 1
  st_age_pcl_rowz <- reactive({
    st_age_pcl <- st_age_pcl()
    st_gsm_pcl <- st_gsm_pcl()
    idx <- apply(st_age_pcl,1,var)
    st_age_pcl_rowz <- t(apply(st_age_pcl, 1, scale))
    rownames(st_age_pcl_rowz) <- rownames(st_gsm_pcl)
    st_age_pcl_rowz <- st_age_pcl_rowz[-which(idx==0),]
    colnames(st_age_pcl_rowz) <- arng()
    st_age_pcl_rowz
  })

  # boot_rho is the output of the bootstrap runs
  # It contains <nboot> rows and <#genes> columns
  boot_rho <- reactive({
    clust <- makeCluster(10) # Initiate cluster
    arng <- arng()
    nboot <- 10
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
        rho <- parApply(clust, boot_pcl, 1, cor, y=boot_age, method="spearman") 
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
    boot_fisherz <- t(apply(boot_rho(), 1, function(x) { scale(0.5*log((1+x)/(1-x))) })) # Fisher z distribution
    bxs_boot_fisherz <- apply(boot_fisherz, 2, quantile, probs=c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0), na.rm=T)
    bxs_boot_fisherz <- data.frame(bxs_boot_fisherz)
    colnames(bxs_boot_fisherz) <- rownames(st_gsm_pcl)
    bxs_boot_fisherz
  })

  # magnitudes of scores of predictive genes
  abs_scores <- eventReactive(input$runpcl,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    all_predg <- all_predg()
    data.frame(abs(bxs_boot_fisherz[4,all_predg]))
  })

  
  # number of genes
  num_genes <- reactive({
    abs_scores <- abs_scores()
    df<-abs_scores[which(as.numeric(abs_scores[1,])>input$score_mag)]
    ncol(df)
  })


  # histogram plot2 of magnitudes
  output$plot2 <- renderPlot({
    .e <- environment()
    df <- cbind(data.frame(all_predg()), t(abs_scores()))
    lo <- input$score_mag
    hi <- max(abs_scores())
    ggplot(data=df, aes(x=df[,2]), environment = .e) + 
      geom_histogram(binwidth=.05,fill=I("gray23"),col=I("gray")) +
      labs(x="abs(score)",y="# of genes") +
      geom_vline(xintercept=as.numeric(input$score_mag), colour="red",size=0.6, linetype="twodash") +
      geom_vline(xintercept=as.numeric(max(abs_scores())), colour="red",size=0.6, linetype="twodash") 
  })

  # slider for plot2
  output$slider_plot2 <- renderUI({
    lwrbound = signif(min(abs_scores()),3)
    uprbound = signif(max(abs_scores()),3)
    sliderInput("score_mag", label = h6("abs(score):"), ticks=FALSE,min = lwrbound, 
                max = uprbound,value=lwrbound)
  })

  # caption for plot2
  output$plot2_caption <- renderText({
    paste("# of genes selected:\n",num_genes())
  })


#
#
#


  # positive predictors
  pos_predg <- reactive({
    bxs_boot_fisherz <- bxs_boot_fisherz()
    pos_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]>0) & (bxs_boot_fisherz[6,]>=1) &
                                                    bxs_boot_fisherz[4,]>=input$score_mag)]
    pos_predg
  })

  # negative predictors
  neg_predg <- reactive({
    bxs_boot_fisherz <- bxs_boot_fisherz()
    neg_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]<0) & (bxs_boot_fisherz[2,]<=-1) &
                                                    abs(bxs_boot_fisherz[4,])>=input$score_mag)]
  })

  # all predictive genes
  all_predg <- eventReactive(input$runpcl, {
    bxs_boot_fisherz <- bxs_boot_fisherz()
    withProgress(message="Filtering genes with high scores", value = 0.1, {
      all_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]>0 & bxs_boot_fisherz[6,]>=1) | 
                                                      (bxs_boot_fisherz[4,]<0 & bxs_boot_fisherz[2,]<=-1))]
      incProgress(0.9)
    })
    all_predg
  })

  # positive predictive genes greater than slider
  select_pos_predg <- eventReactive(input$rungt,{
    withProgress(message = "Generating positive GO terms\n", detail = "Building most specific GOs", value = 0, {
      bxs_boot_fisherz <- bxs_boot_fisherz()
      all <- colnames(bxs_boot_fisherz)[which(((bxs_boot_fisherz[4,]>0 & bxs_boot_fisherz[6,]>=1)&
                                                 abs(bxs_boot_fisherz[4,])>input$score_mag))]
    })
    all
  })

  # negative predictive genes greater than slider
  select_neg_predg <- eventReactive(input$rungt,{
    withProgress(message = "Generating negative GO terms\n", detail = "Building most specific GOs", value = 0, {
      bxs_boot_fisherz <- bxs_boot_fisherz()
      all <- colnames(bxs_boot_fisherz)[which(((bxs_boot_fisherz[4,]<0 & bxs_boot_fisherz[2,]<=-1)) &
                                                abs(bxs_boot_fisherz[4,])>input$score_mag)]
    })
    all
  })

  # positive pcl table
  pos_pcl <- eventReactive(input$tablepcl,{
    withProgress(message = 'Generating positive data', detail = "0/3", value = 0, {
      bxs_boot_fisherz <- bxs_boot_fisherz()
      st_age_pcl_rowz <- st_age_pcl_rowz()
      incProgress(1/3, detail = "1/3")
      pos_predg <- pos_predg()
      pos_predg <- intersect(pos_predg, rownames(st_age_pcl_rowz))
      incProgress(1/3, detail = "2/3")
      pos_predg_pcl_rowz <- cbind(gene_sym()[pos_predg,],
                                  t(bxs_boot_fisherz[4,pos_predg]), 
                                  apply(st_age_pcl_rowz[pos_predg,], 2, sprintf, fmt="%.6f"))
      colnames(pos_predg_pcl_rowz) <- c("symbol", "name","score", arng())
      incProgress(1/3, detail = "3/3")
    })
    pos_predg_pcl_rowz[order(-pos_predg_pcl_rowz$score,pos_predg_pcl_rowz$name),]
  })

  # negative pcl table
  neg_pcl <- eventReactive(input$tablepcl,{
    withProgress(message = 'Generating negative data', detail = "0/3", value = 0, {
      bxs_boot_fisherz <- bxs_boot_fisherz()
      st_age_pcl_rowz <- st_age_pcl_rowz()
      incProgress(1/3, detail = "1/3")
      neg_predg <- neg_predg()
      neg_predg <- intersect(neg_predg, rownames(st_age_pcl_rowz))
      incProgress(1/3, detail = "2/3")
      neg_predg_pcl_rowz <- cbind(gene_sym()[neg_predg,],
                                  t(bxs_boot_fisherz[4,neg_predg]), 
                                  apply(st_age_pcl_rowz[neg_predg,], 2, sprintf, fmt="%.6f"))
      colnames(neg_predg_pcl_rowz) <- c("symbol", "name","score", arng())
      incProgress(0.25, detail = "3/3")
    })
    neg_predg_pcl_rowz[order(neg_predg_pcl_rowz$score,neg_predg_pcl_rowz$name),]
  })




#
#
#




  #outputs
  output$ptable <- DT::renderDataTable({
    pos <- pos_pcl()
    cap <- paste('Table 1: Genes with the most positive Spearman correlation scores for ages ', lwr(), ' to ', upr())
    DT::datatable(pos[,1:3], rownames=FALSE, caption = cap)
  })

  output$ntable <- DT::renderDataTable({
    neg <- neg_pcl()
    cap <- paste('Table 2: Genes with the most negative Spearman correlation scores for ages ', lwr(), ' to ', upr())
    DT::datatable(neg[,1:3], rownames=FALSE, caption = cap)
  })

  output$ppcl <- DT::renderDataTable({
    pos <- pos_pcl()
    cap <- paste('Table 3: Expression values of genes with the most positive Spearman correlation scores for ages ', lwr(), ' to ', upr())
    DT::datatable(pos, rownames=FALSE, caption = cap)
  })

  output$npcl <- DT::renderDataTable({
    neg <- neg_pcl()
    cap <- paste('Table 4: Expression values of genes with the most negative Spearman correlation scores for ages ', lwr(), ' to ', upr())
    DT::datatable(neg, rownames=FALSE, caption = cap)
  })

  output$posheat <- renderD3heatmap({
    pos_p <- pos_pcl()
    rownames(pos_p) <- pos_p[,1]
    pos_p <- pos_p[,-c(1,2,3)]
    # convert factor to numeric data frame
    indx <- sapply(pos_p, is.factor)
    pos_p[indx] <- lapply(pos_p[indx], function(x) as.numeric(as.character(x)))
    length <- 50
    if (dim(pos_p)[1] < 50) length <- dim(pos_p)[1] 
    d3heatmap(pos_p[1:length,],Colv = FALSE,xaxis_font_size=7,yaxis_font_size=7)
  })

  output$negheat <- renderD3heatmap({
    neg_p <- neg_pcl()
    rownames(neg_p) <- neg_p[,1]
    neg_p <- neg_p[,-c(1,2,3)]
    # convert factor to numeric data frame
    indx <- sapply(neg_p, is.factor)
    neg_p[indx] <- lapply(neg_p[indx], function(x) as.numeric(as.character(x)))
    length <- 50
    if (dim(neg_p)[1] < 50) length <- dim(neg_p)[1] 
    d3heatmap(neg_p[1:length,],Colv = FALSE,xaxis_font_size=7,yaxis_font_size=7)
  })

  # gene ontology of positive predictive genes
  pos_go <- reactive({
    all_predg <- select_pos_predg()
    withProgress(message = "Generating positive GO terms\n", detail = "Compiling tables", value = 0.1, {
      gsm_pcl <- gsm_pcl()
      incProgress(0.1, detail = "Compiling tables")
      ann <- annFUN.org("BP",feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
      incProgress(0.2, detail = "Building GO DAG topology and annotations")
      geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg))
      names(geneList) <- rownames(gsm_pcl)
      selector <- function(x) {return (x==0)}
      sampleGOdata <- new("topGOdata",ontology="BP",allGenes=geneList,
                          geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
      incProgress(0.2, detail = "Running Fisher test")
      resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
      incProgress(0.2, detail = "Running Kolmogorov–Smirnov test")
      resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
      # resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
      incProgress(0.2, detail = "Aggregating results")
      allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                         classicKS = resultKS, 
                         #                    elimKS = resultKS.elim,
                         orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 50)
    })
    allRes
  })

  neg_go <- reactive({
    all_predg <- select_neg_predg()
    withProgress(message = "Generating negative GO terms\n", detail = "Compiling tables", value = 0.1, {
      gsm_pcl <- gsm_pcl()
      incProgress(0.1, detail = "Compiling tables")
      ann <- annFUN.org("BP",feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
      incProgress(0.2, detail = "Building GO DAG topology and annotations")
      geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg))
      names(geneList) <- rownames(gsm_pcl)
      selector <- function(x) {return (x==0)}
      sampleGOdata <- new("topGOdata",ontology="BP",allGenes=geneList,
                          geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
      incProgress(0.2, detail = "Running Fisher test")
      resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
      incProgress(0.2, detail = "Running Kolmogorov–Smirnov test")
      resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
      # resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
      incProgress(0.2, detail = "Aggregating results")
      allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                         classicKS = resultKS, 
                         #                    elimKS = resultKS.elim,
                         orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 50)
    })
    allRes
  })

  output$pos_goterms <- DT::renderDataTable({
    allRes <- pos_go()
    cap <- paste('Table 5: Top 50 GO terms enriched in positively correlated genes')
    DT::datatable(allRes, rownames=TRUE, caption = cap,options = list(columnDefs = list(list(
      targets = 2,
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.length > 20 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
        "}")
    ))))
  })

  # gene ontology of negative predictive genes
  output$neg_goterms <- DT::renderDataTable({
    allRes <- neg_go()
    cap <- paste('Table 6: Top 50 GO terms enriched in negatively correlated genes')
    DT::datatable(allRes, rownames=TRUE, caption = cap,options = list(columnDefs = list(list(
      targets = 2,
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.length > 20 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
        "}")
    ))))
  })

  output$ptable_dl <- downloadHandler(
    filename = function(){"pos_genes.csv"},
    content = function(file){write.csv(pos_pcl()[,1:3],file)}
  )

  output$ntable_dl <- downloadHandler(
    filename = function(){"neg_genes.csv"},
    content = function(file){write.csv(neg_pcl()[,1:3],file)}
  )

  output$ppcl_dl <- downloadHandler(
    filename = function(){"pos_genes_scores.csv"},
    content = function(file){write.csv(pos_pcl(),file)}
  )

  output$npcl_dl <- downloadHandler(
    filename = function(){"neg_genes_scores.csv"},
    content = function(file){write.csv(neg_pcl(),file)}
  )

  output$pgo_dl <- downloadHandler(
    filename = function(){"pos_genes_goterms.csv"},
    content = function(file){write.csv(pos_go(),file)}
  )

  output$ngo_dl <- downloadHandler(
    filename = function(){"neg_genes_goterms.csv"},
    content = function(file){write.csv(neg_go(),file)}
  )
})
  
  
