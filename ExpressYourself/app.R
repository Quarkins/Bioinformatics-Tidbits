
library(shiny)
#library(ggplot2)
library(limma)

#Load in the data
MCRI_data <- readRDS("Data/MCRI.Rds")
MCRI_target <- readRDS("Data/MCRI_target.Rds")


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("The Hunt for Expression Outliers"),
   
   #A box to search by gene name
   selectizeInput("gene",label="Gene(s)",
               choices=NULL,multiple=TRUE),
   
   #a box to search by ALL sub-type
   selectizeInput("type",label="Sub-Type",
              choices=NULL,multiple=FALSE),
   
  
    mainPanel(
           plotOutput("expression"),
           plotOutput("boxplot")
       )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    d2 <- reactive({
        if(!is.null(input$type)){
            selected = MCRI_target[MCRI_target$Type == input$type,]
            m <- match(selected$Sample,colnames(MCRI_data$E))
            m <-m[!is.na(m)]
            MCRI_data$E = MCRI_data$E[,m]
        }
        
        MCRI_data
        
    })
    
    
     output$expression <- renderPlot({
            
            if(length(input$gene) < 1){}   
         
            else{
                if(length(input$gene) >1) {
                stripchart(MCRI_data$E[input$gene,] ~ row.names(MCRI_data$E[input$gene,]) ,ylab="Log CPM",
                       pch=16,method="jitter",col="orange",vertical=TRUE)
                }
                else{
                    stripchart(MCRI_data$E[input$gene,],xlab="Log CPM", ylab=input$gene,
                           pch=16,method="jitter",col="orange")
                }
            }

    })
   
     output$boxplot <- renderPlot({
        
            
            if(length(input$gene) < 1){}
            else{
                sub = t(MCRI_data$E[input$gene,])
                bp = boxplot(sub,ylab="Log CPM",xlab="Gene(s)")
                bp
                #text( bp$group, bp$out,
                #     rownames(sub)[which(sub == bp$out)])
            }
     })
     
     updateSelectizeInput(session,"gene",
                          choices=row.names(MCRI_data$E),server=TRUE)
     
     updateSelectizeInput(session,"type",
                          choices=unique(MCRI_target$Cell_type$Type),server=TRUE)
     
}

# Run the application 
shinyApp(ui = ui, server = server)

