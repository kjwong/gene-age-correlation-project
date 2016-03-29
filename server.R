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
  # list of samples --INPUT
  st_samples <- reactive({
    st_samples <- read.table("race_healthy1-samples.txt")
    as.vector(st_samples$V1)
  })
  gene_sym <- reactive({
    read.delim("human_gene-info_ncbi.txt", header=T, row.names=1, sep="\t")
  })
  gsm_pcl <- reactive({
    data.frame(fread(paste0('blood.', 'f', '.pcl'), header=T), row.names=1)
  })
  gsm_age <- reactive({
    gsm_age <- read.table(paste0('blood.', 'f', '.txt'), row.names=1)
    colnames(gsm_age) <- c('age','gse')
    gsm_age
  })
#   st_samples <- eventReactive(input$upload4, {
#     inFile <- input$file4
#     if (is.null(inFile))
#       return(NULL)   
#     t <- read.table(inFile$datapath)
#     as.vector(t$V1)
#   })
  # gene ids, symbols, names
#   gene_sym <- eventReactive(input$upload2, {
#     inFile <- input$file2
#     if (is.null(inFile))
#       return(NULL)
#     read.delim(inFile$datapath, sep=input$sep2,header=input$header2,row.names=1) 
#   })
#   sex <- reactive ({
#     if (input$sex == c(1,2)) c('m','f')
#     else if (input$sex == 1) 'm'
#     else 'f'
#   })
  
  # reading-in the entire pcl (both healthy and disease) --INPUT TISSUE,SEX
#   gsm_pcl <- eventReactive(input$upload1, {
#     inFile <- input$file1
#     if (is.null(inFile))
#       return(NULL)   
#     data.frame(fread(inFile$datapath, header = input$header1), row.names = 1)
#   })
  
  # reading-in sample age annotations --INPUT TISSUE,SEX
#   gsm_age <- eventReactive(input$upload3, {
#     inFile <- input$file3
#     if (is.null(inFile))
#       return(NULL)   
#     read.table(inFile$datapath, header=input$header3, row.names=1,col.names=c('age','gse'))
#   })
  lwr <- reactive({input$age_range[1]})
  upr <- reactive({input$age_range[2]})
  arng <- reactive({lwr():upr()})
  
  output$plot_caption <- renderText({
    gsm_age <- gsm_age()
    paste("# of samples selected:\n",dim(gsm_age[which((gsm_age$age>=lwr()) & (gsm_age$age<=upr())),])[1])
  })
  output$slider_plot <- renderUI({
    gsm_age <- gsm_age()
    sliderInput("age_range", label = h6("Age range:"), min = min(gsm_age[,1]), max = max(gsm_age[,1]), value = c(30, 40))
  })
  output$plot <- renderPlot({
    .e <- environment()
    gsm_age <- gsm_age()
    p <- ggplot(data=gsm_age, aes(gsm_age[,1]), environment = .e) + 
      geom_histogram(binwidth=2,fill=I("black"),col=I("gray")) + 
      labs(x="age",y="# of GSM samples") 
    p + geom_vline(xintercept=as.numeric(lwr()), colour="red",size=0.6, linetype="solid") +
      geom_vline(xintercept=as.numeric(upr()), colour="red",size=0.6, linetype="solid") 
#       annotate("text", label = paste("total #\n",
#               dim(gsm_age[which((gsm_age$age>=lwr()) & (gsm_age$age<=upr())),])[1]), 
#             x = lwr() + (upr()-lwr())/2, y = 100, size = 4, colour = "red")
  })
  # intersecting col names with sample ids
  st_make <- function(){
    gsm_age <- gsm_age()
    gsm_pcl <- gsm_pcl()
    intersect(st_samples(),
              intersect(rownames(gsm_age)[which((gsm_age$age>=lwr()) & (gsm_age$age<=upr()))],
                        colnames(gsm_pcl)))
  }
    
  # creating pcl for samples
  st_gsm_pcl <- reactive({
    gsm_pcl <- gsm_pcl()
    st_samples <- st_make()
    gsm_pcl[,st_samples]
  })
  # gsm age for samples
  st_gsm_age <- reactive({
    gsm_age <- gsm_age()
    st_samples <- st_make()
    gsm_age[st_samples,]
  })
  # creating pcl of genes per age with median of gene-expression values per gene
  st_age_pcl <- reactive({
    st_gsm_age <- st_gsm_age()
    st_gsm_pcl <- st_gsm_pcl()
    arng <- arng()
    st_age_pcl <- array(NaN, c(nrow(st_gsm_pcl()), length(lwr():upr())))
    for(j in 1:length(arng)) {
      age_gsm <- rownames(st_gsm_age)[which(st_gsm_age$age==arng[j])]
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
#     print(st_age_pcl["6888",])
#     print(apply(data.frame(st_age_pcl["6888",]),2,scale))
    
    st_age_pcl_rowz <- t(apply(st_age_pcl, 1, scale))
    st_age_pcl_rowz[idx==0] <- st_age_pcl[idx==0]
    colnames(st_age_pcl_rowz) <- arng()
    rownames(st_age_pcl_rowz) <- rownames(st_gsm_pcl)
#     print(st_age_pcl_rowz["6888",])
    st_age_pcl_rowz
  })
  
  # boostrap runs
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
        for(j in 1:length(arng)) {
          age_gsm <- rownames(st_gsm_age)[which(st_gsm_age$age==arng[j])]
          if(length(age_gsm)==1) {
            bag_gsm <- c(bag_gsm, age_gsm)
          } else {
            sel_gsm <- sample(age_gsm, replace=T) # randomness
            bag_gsm <- c(bag_gsm, unique(sel_gsm))
          }
        }
        boot_pcl <- st_gsm_pcl[,bag_gsm]
        boot_age <- st_gsm_age[bag_gsm,]$age
        clusterExport(clust,varlist = c("boot_pcl","boot_age"),envir=environment())
        rho <- parApply(clust, boot_pcl, 1, cor, y=boot_age, method="spearman")
        boot_rho[n,] <- rho
      }
    })
    stopCluster(clust)
    colnames(boot_rho) <- rownames(st_gsm_pcl)
    boot_rho
  })
  # processing rho into comparable correlation scores (fisherz)
  bxs_boot_fisherz <- reactive({
    st_gsm_pcl <- st_gsm_pcl()
    boot_fisherz <- t(apply(boot_rho(), 1, function(x) { scale(0.5*log((1+x)/(1-x))) }))
    bxs_boot_fisherz <- apply(boot_fisherz, 2, quantile, probs=c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0), na.rm=T)
    bxs_boot_fisherz <- data.frame(bxs_boot_fisherz); colnames(bxs_boot_fisherz) <- rownames(st_gsm_pcl)
    bxs_boot_fisherz
  })
  # positive predictors
  pos_predg <- reactive({
    bxs_boot_fisherz <- bxs_boot_fisherz()
    pos_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]>0) & (bxs_boot_fisherz[6,]>=1) &
                                                    bxs_boot_fisherz[4,]>=input$score_mag)]
#     print("6888" %in% pos_predg)
    pos_predg
  })
  # negative predictors
  neg_predg <- reactive({
    bxs_boot_fisherz <- bxs_boot_fisherz()
    neg_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]<0) & (bxs_boot_fisherz[2,]<=-1) &
                                                    abs(bxs_boot_fisherz[4,])>=input$score_mag)]
  })
  # positive pcl table
  pos_pcl <- eventReactive(input$tablepcl,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    st_age_pcl_rowz <- st_age_pcl_rowz()
    pos_predg <- pos_predg()
#     print(gene_sym()["6888",])
#     print(t(bxs_boot_fisherz[4,"6888"]))
#     print(is.matrix(st_age_pcl_rowz))
#     print(st_age_pcl_rowz["6888",])
#     print(apply(st_age_pcl_rowz["6888",], 2, sprintf, fmt="%.6f"))
    pos_predg_pcl_rowz <- cbind(gene_sym()[pos_predg,],
                                t(bxs_boot_fisherz[4,pos_predg]), 
                                apply(st_age_pcl_rowz[pos_predg,], 2, sprintf, fmt="%.6f"))
    colnames(pos_predg_pcl_rowz) <- c("symbol", "name","score", arng())
    pos_predg_pcl_rowz[order(-pos_predg_pcl_rowz$score,pos_predg_pcl_rowz$name),]
  })
  # negative pcl table
  neg_pcl <- eventReactive(input$tablepcl,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    st_age_pcl_rowz <- st_age_pcl_rowz()
    neg_predg <- neg_predg()
    neg_predg_pcl_rowz <- cbind(gene_sym()[neg_predg,],
                                t(bxs_boot_fisherz[4,neg_predg]), 
                                apply(st_age_pcl_rowz[neg_predg,], 2, sprintf, fmt="%.6f"))
    colnames(neg_predg_pcl_rowz) <- c("symbol", "name","score", arng())
    neg_predg_pcl_rowz[order(neg_predg_pcl_rowz$score,neg_predg_pcl_rowz$name),]
  })
  # all predictive genes
  all_predg <- eventReactive(input$runpcl, {
     bxs_boot_fisherz <- bxs_boot_fisherz()
     all_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]>0 & bxs_boot_fisherz[6,]>=1) | 
                                                     (bxs_boot_fisherz[4,]<0 & bxs_boot_fisherz[2,]<=-1))]
  })
  select_pos_predg <- eventReactive(input$rungt,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    colnames(bxs_boot_fisherz)[which(((bxs_boot_fisherz[4,]>0 & bxs_boot_fisherz[6,]>=1)&
                                       abs(bxs_boot_fisherz[4,])>input$score_mag))]
  })
  select_neg_predg <- eventReactive(input$rungt,{
    bxs_boot_fisherz <- bxs_boot_fisherz()
    colnames(bxs_boot_fisherz)[which(((bxs_boot_fisherz[4,]<0 & bxs_boot_fisherz[2,]<=-1)) &
                                       abs(bxs_boot_fisherz[4,])>input$score_mag)]
  })
  output$pos_goterms <- DT::renderDataTable({
    all_predg <- select_pos_predg()
    gsm_pcl <- gsm_pcl()
    ann <- annFUN.org("BP",feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
    geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg))
    names(geneList) <- rownames(gsm_pcl)
    selector <- function(x) {return (x==0)}
    sampleGOdata <- new("topGOdata",ontology="BP",allGenes=geneList,
                        geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
    # resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
    allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                       classicKS = resultKS, 
                       #                    elimKS = resultKS.elim,
                       orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 50)
    cap <- paste('Table 5: Top 50 GO terms enriched in positively correlated genes')
    DT::datatable(allRes, rownames=TRUE, caption = cap)
  })
  output$neg_goterms <- DT::renderDataTable({
    all_predg <- select_neg_predg()
    gsm_pcl <- gsm_pcl()
    ann <- annFUN.org("BP",feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
    geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg))
    names(geneList) <- rownames(gsm_pcl)
    selector <- function(x) {return (x==0)}
    sampleGOdata <- new("topGOdata",ontology="BP",allGenes=geneList,
                        geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
    # resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
    allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                       classicKS = resultKS, 
                       #                    elimKS = resultKS.elim,
                       orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 50)
    cap <- paste('Table 6: Top 50 GO terms enriched in negatively correlated genes')
    DT::datatable(allRes, rownames=TRUE, caption = cap)
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
#   num_pos_genes <- reactive({
#     length(pos_predg())
#   })
#   num_neg_genes <- reactive({
#     length(neg_predg())
#   })
  # histogram plot of magnitudes
  output$plot2 <- renderPlot({
    .e <- environment()
    df <- cbind(data.frame(all_predg()), t(abs_scores()))
    lo <- input$score_mag
    hi <- max(abs_scores())
    ggplot(data=df, aes(x=df[,2]), environment = .e) + 
      geom_histogram(binwidth=.05,fill=I("black"),col=I("gray")) +
      labs(x="abs(score)",y="# of genes") +
      geom_vline(xintercept=as.numeric(input$score_mag), colour="red",size=0.6, linetype="solid") +
      geom_vline(xintercept=as.numeric(max(abs_scores())), colour="red",size=0.6, linetype="solid") 
#       annotate("text", label = paste("total #\n",num_genes()),
#               x = lo + (hi-lo)/2, y = 150, size = 4, colour = "red")
  })
  output$slider_plot2 <- renderUI({
    lwrbound = signif(min(abs_scores()),3)
    uprbound = signif(max(abs_scores()),3)
    sliderInput("score_mag", label = h6("abs(score):"), ticks=FALSE,min = lwrbound, 
                max = uprbound,value=lwrbound)
  })
  output$plot2_caption <- renderText({
    paste("# of genes selected:\n",num_genes())
  })
  output$tablescap <- renderUI({
    h3("Genes with the most significant correlation:")
  })
  output$ptable <- DT::renderDataTable({
    pos <- pos_pcl()
    cap <- paste('Table 1: Genes with the most positive Spearman correlation scores for ages ', lwr(), ' to ', upr(),'.')
    DT::datatable(pos[,1:3], rownames=FALSE, caption = cap)
  })
  output$ntable <- DT::renderDataTable({
    neg <- neg_pcl()
    cap <- paste('Table 2: Genes with the most negative Spearman correlation scores for ages ', lwr(), ' to ', upr(),'.')
    DT::datatable(neg[,1:3], rownames=FALSE, caption = cap)
  })
  output$ppcl <- DT::renderDataTable({
    pos <- pos_pcl()
    cap <- paste('Table 3: Expression values of genes with the most positive Spearman correlation scores for ages ', lwr(), ' to ', upr(),'.')
    DT::datatable(pos, rownames=FALSE, caption = cap)
  })
  output$npcl <- DT::renderDataTable({
    neg <- neg_pcl()
    cap <- paste('Table 4: Expression values of genes with the most negative Spearman correlation scores for ages ', lwr(), ' to ', upr(),'.')
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
  
})
  
  
