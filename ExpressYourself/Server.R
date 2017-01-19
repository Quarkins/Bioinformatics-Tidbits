library(shiny)
library(limma)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(plotly)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("The Hunt for Expression Outliers"),
    
    #A box to search by gene name
    selectizeInput("gene",label="Gene(s)",
                   choices=NULL,multiple=TRUE),
    
    #a box to search by ALL sub-type
    selectizeInput("type",label="Sub-Type",
                   choices=NULL,multiple=FALSE),
    
    
    mainPanel(
        plotlyOutput("expression"),
        plotOutput("boxplot"),
        plotOutput("density"),
        #tableOutput("outliers"),
        tableOutput("classified"),
        plotOutput("MDS")
    )
))