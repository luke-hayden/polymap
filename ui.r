library(seqinr)
library(ggplot2)
library(shiny)
library(stringr)
library(zoo)

shinyUI(fluidPage(
  titlePanel("Parhyale tools"),
  sidebarLayout(
    sidebarPanel(
      p("Polymorphism mapper"),
      p("version 0.3"),
      p("Created by Luke Hayden"),
      br(),
      br(),
      img(src = "parhy.png", height = 101, width = 199, align="left"),
      br(),
      br(),
      br(),
      br(),
      br(),
      textInput("text", label = h3("Enter a gene model ID below"), value = ""),
      p("Example: phaw_30_tra_m.000001"),
      br()
    ),
    mainPanel(
      h1("Polymorphism mapper"),
      tabsetPanel(
        tabPanel("Sequence",
                 br(),
                 textOutput("text1"),
                 textOutput("isfound"),
                 br(),
                 h4("Sequence:"),
                 textOutput("msx"),
                 br(),
                 p("Sequence Length:"),
                 textOutput("length"),
                 downloadButton('myseq.fasta', 'Download')),
        tabPanel("Polymorphism Plot",
                 plotOutput("polyplot"), 
                 checkboxInput("checkboxpoly", label = "Plot mean polymorphism", value = FALSE), 
                 h4("Set plot limits:"),
                 column(6, 
                        "",
                        textInput("startrange", "Startpoint"), 
                        textInput("endrange", "Endpoint")
                 ),
                 column(6, 
                        "",
                        checkboxInput("checkboxseqcut", label= "Limit plot to sequence of interest", value=FALSE ), 
                        textInput("searchstring", "Sequence of interest" )
                 )                 ),
        tabPanel("Polymorphic sites",
                 h4("Known polymorphic sites:"),
                 textOutput("pos"), 
                 br(),
                 p("Number of polymorphisms:"),
                 textOutput("numpoly")),
        tabPanel("Sequence with Polymorphisms",
                 p("All known polymorphic sites replaced with N"),
                 h4("Sequence:"),
                 textOutput("msxp"),
                 downloadButton('mypolyseq.fasta', 'Download'))
      )
    )
  )
))