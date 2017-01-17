library(shiny)
library(limma)


#Load in the data
MCRI_data <- readRDS("Data/MCRI.Rds")
MCRI_target <- readRDS("Data/MCRI_target.Rds")

shinyServer(function(input, output, session) {
    
    dd <- reactive({
        if(input$type != "All"){
            selected = MCRI_target[MCRI_target$Cell_Type == input$type,]
            m <- match(selected$Sample,colnames(MCRI_data))
            m <-m[!is.na(m)]
            MCRI_data[,m,drop=FALSE]
        }
        else{
            MCRI_data
        }
    })
    
    output$expression <- renderPlot({
        
        if(length(input$gene) < 1){}   
        
        else{
            if(length(input$gene) >1) {
                stripchart(dd()[input$gene,] ~ row.names(dd()[input$gene,]) ,ylab="Log CPM",
                           pch=16,method="jitter",col="orange",vertical=TRUE)
            }
            else{
                stripchart(dd()[input$gene,],xlab="Log CPM", ylab=input$gene,
                           pch=16,method="jitter",col="orange")
                
            }
            
        }
     
    })
    
 
    output$boxplot <- renderPlot({
        
        
        if(length(input$gene) < 1){}
        else if(length(input$gene) == 1){
            sub = t(dd()[input$gene,,drop=FALSE])
            bp = boxplot(sub,ylab="Log CPM",xlab=input$gene,pch=16,outcol='red',col=rgb(226/255,9/255,172/255,0.2),border="#e209ac")
            bp 
            if(length(bp$out) > 0){
                text(rep(1,length(bp$out)),bp$out,
                     rownames(sub)[match(bp$out,sub)])
            }
        }
        else{
            sub = t(dd()[input$gene,])
            bp = boxplot(sub,ylab="Log CPM",xlab="Gene(s)",pch=16,outcol='red',col=rgb(226/255,9/255,172/255,0.4),border="#e209ac")
            bp
        }
    })
    
    output$density <- renderPlot({
        sub = dd()[input$dense_gene,,drop=FALSE]
        d<- density(as.numeric(sub[input$dense_gene,]))
        bp<-plot(d,xlab="logFPKM",main = input$dense_gene,col="#4286f4",lwd=2)
        bp 
    })
    
    output$outliers <- renderTable({
        #Get outliers
        sub = t(dd()[input$dense_gene,,drop=FALSE])
        bp = boxplot(sub,plot=FALSE)
        tab = data.frame(outlier_sample = rownames(sub)[match(bp$out,sub)], logFPKM = bp$out)
        tab
    })
    
    
    updateSelectizeInput(session,"gene",
                         choices=row.names(MCRI_data),server=TRUE)
    
    updateSelectizeInput(session,"type",
                         choices=c(unique(MCRI_target$Cell_Type),"All"),
                         server=TRUE,selected='All')
    
    updateSelectizeInput(session,"dense_gene",
                         choices=row.names(MCRI_data),
                         server=TRUE,selected='AHR')
    
})