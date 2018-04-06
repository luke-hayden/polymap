library(shiny)
library(seqinr)
library(stringr)
library(ggplot2)
library(zoo)

hetdat <- readRDS("data/hetdata.rds")
seqs <- readRDS("data/seqs.rds")

shinyServer(function(input, output) {
  output$text1 <- renderText({
    paste("You have entered: ", input$text)
  })
  output$isfound <- renderText({
    if(input$text %in% hetdat$geneID) "This gene ID has been found" else "Gene ID not found"
  })
  output$msx <- renderText({
    req(input$text)
    gsub("(.{80})", "\\1 ", gsub(" |,", "", toString(getSequence(seqs[names(seqs) %in% input$text[1]][[1]]))))
  })
  output$mdat <- renderTable({
    req(input$text)
    subset(hetdat, hetdat$geneID %in% input$text)
  })
  output$length <- renderText({
    req(input$text)
    paste(subset(hetdat, hetdat$geneID %in% input$text)[1,2], "bp")
  })
  output$l <- renderText({
    as.integer(subset(hetdat, hetdat$geneID %in% input$text)[1,2])
  })
  output$varcode <- renderText({
    req(input$text)
    toString(subset(hetdat, hetdat$geneID %in% input$text)[,4])
  })
  output$pos <- renderText({
    req(input$text)
    sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]])
  })
  output$numpoly <- renderText({
    req(input$text)
    length(sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]]))
  })
  output$polytable <- renderTable({
    req(input$text)
    data.frame(poly=sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]]), num=seq(1,length(sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]])), 1))
  })
  output$polyplot <- renderPlot({
    req(input$text)
    mdata<- data.frame(index = seq(1,subset(hetdat, hetdat$geneID %in% input$text)[1,2],1), 
                      ispoly = seq(1,subset(hetdat, hetdat$geneID %in% input$text)[1,2],1) %in% sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]]),
                      y=1, 
                      polynum = ifelse(seq(1,subset(hetdat, hetdat$geneID %in% input$text)[1,2],1) %in% sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]]), 1, 0.7),
                      seq=getSequence(seqs[names(seqs) %in% input$text[1]][[1]]))
    seqreg <- str_locate(gsub(" |,", "", toString(getSequence(seqs[names(seqs) %in% input$text[1]][[1]]))), gsub(" ", "",input$searchstring))
    mdata$rollm <- c(rep(0.7,24), rollmean(mdata$polynum, 50), rep(0.7,25))
    ggplot(mdata, aes(x=index, y=y, colour=ispoly, shape=ispoly, size=ispoly))+ 
      theme_bw()+
      geom_point()+
      scale_shape_manual(values=c(15,124))+
      scale_colour_manual(values=c("grey", "dark blue"))+
      scale_size_manual(values=c(1,15)) +
      xlab("Position(bp)")+
      ylab("")+
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), panel.grid.major.y =element_blank(), panel.grid.minor.y =element_blank()) +
      theme(legend.position = "none")+
      ylim(0.5,1.5)+
      geom_text(aes(x=index, y=y+0.2, label=seq), size=5)+
      annotate("text", x = max(mdata$index)/4, y = 1.5, label = paste("Polymorphism of", input$text))+
      annotate("text", x=0,y=0.7, label="0%", alpha=(if(input$checkboxpoly == TRUE)0.8 else 0), colour="orange") +
      annotate("text", x=0,y=0.73, label="10%", alpha=(if(input$checkboxpoly == TRUE)0.8 else 0), colour="orange") +
      annotate("text", x=0,y=0.76, label="20%", alpha=(if(input$checkboxpoly == TRUE)0.8 else 0), colour="orange") +
      annotate("text", x=0,y=0.79, label="30%", alpha=(if(input$checkboxpoly == TRUE)0.8 else 0), colour="orange") +
      annotate("rect", 
               xmin = seqreg[1,1], xmax = seqreg[1,2], ymin = 0.85, ymax = 1.15, alpha=0.4, fill="blue")+
      geom_line(aes(x=index, y=rollm), inherit.aes = FALSE, alpha=(if(input$checkboxpoly == TRUE)0.5 else 0), size=2, colour="orange")+
      scale_x_continuous(limits = if(
        input$checkboxseqcut == FALSE)
        c( as.integer(input$startrange),as.integer(input$endrange ))
        else c( seqreg[1,1], seqreg[1,2] ))
  })
  output$msxp <- renderText({
    req(input$text)
    polyseq <- getSequence(seqs[names(seqs) %in% input$text[1]][[1]])
    polyseq[seq(1,subset(hetdat, hetdat$geneID %in% input$text)[1,2],1)%in%sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]])] <- "N"
    gsub("(.{80})", "\\1 ", gsub(" |,", "", toString(polyseq)))
  })
  output$myseq.fasta <- downloadHandler(
    filename = function() { paste(toString(input$text), '.fasta', sep='') },
    content = function(file) {
      write(paste(">", input$text, " ", gsub(" |,", "", toString(getSequence(seqs[names(seqs) %in% input$text[1]][[1]]))), sep=""), file, append=FALSE)
  })
  output$mypolyseq.fasta <- downloadHandler(
    filename = function() { paste(toString(input$text), '.fasta', sep='') },
    content = function(file) {
      polyseq <- getSequence(seqs[names(seqs) %in% input$text[1]][[1]])
      polyseq[seq(1,subset(hetdat, hetdat$geneID %in% input$text)[1,2],1)%in%sub(":...", "", strsplit(toString(subset(hetdat, hetdat$geneID %in% input$text)[,4]), ",")[[1]])] <- "N"
      write(paste(">", input$text, " ", gsub(" |,", "", toString(polyseq)), sep=""), file, append=FALSE)
    })
})
