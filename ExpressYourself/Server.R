library(shiny)
library(ggplot2)


#kinase <- readRDS('ShinyExpress/data/kinase.rds')

function(input, output) {
  
  dataset <- reactive({
    
   m < match(tt$Gene,input$gene)
   tt[m,]
    #if(input$x != '')
    #  tt[tt$Gene==input$x,]
  })
  
  output$plot <- renderPlot({
    
    p <- ggplot(dataset(),aes(factor(Gene),LogCounts,fill=factor(Gene))) #So here we colour and factor by gene (so one boxplot per gene)
    p <- p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +xlab("Gene") + ylab("LogCPM") #Make sure to label the plot
    
    
    if (input$jitter)
      p <- p + geom_jitter(size=1,aes(color=factor(Gene),alpha=0.5,text=Sample))
    
    print(p)
    
  }, height=700)
  
}