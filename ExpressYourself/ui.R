library(shiny)
library(limma)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("The Hunt for Expression Outliers"),
    
    #A box to search by gene name
    selectizeInput("gene",label="Gene(s)",
                   choices=NULL,multiple=TRUE),
    
    #a box to search by ALL sub-type
    selectizeInput("type",label="Sub-Type",
                   choices=NULL,multiple=FALSE),
    
    selectizeInput("dense_gene",label="Gene for density plot",
                   choices=NULL,multiple=FALSE),
    
    mainPanel(
        plotOutput("expression"),
        plotOutput("boxplot"),
        plotOutput("density"),
        tableOutput("outliers")
    )
))