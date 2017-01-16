library(shiny)
library(limma)


#Load in the data
MCRI_data <- readRDS("Data/MCRI.Rds")
MCRI_target <- readRDS("Data/MCRI_target.Rds")

shinyServer(function(input, output, session) {
    
    dd <- reactive({
        if(input$type != "All"){
            selected = MCRI_target[MCRI_target$Cell_type == input$type,]
            m <- match(selected$Sample,colnames(MCRI_data$E))
            m <-m[!is.na(m)]
            MCRI_data$E[,m,drop=FALSE]
        }
        else{
            MCRI_data$E
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
            bp = boxplot(sub,ylab="Log CPM",xlab=input$gene,pch=16,outcol='red')
            bp 
            if(length(bp$out) > 0){
                text(rep(1,length(bp$out)),bp$out,
                     rownames(sub)[match(bp$out,sub)])
            }
        }
        else{
            sub = t(dd()[input$gene,])
            bp = boxplot(sub,ylab="Log CPM",xlab="Gene(s)",pch=16,outcol='red')
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
        bp = boxplot(sub)
        tab = data.frame(outlier_sample = rownames(sub)[match(bp$out,sub)], logFPKM = bp$out)
        tab
    })
    
    
    updateSelectizeInput(session,"gene",
                         choices=row.names(MCRI_data$E),server=TRUE)
    
    updateSelectizeInput(session,"type",
                         choices=c("Non T-ALL","T-ALL","Excluded (Mostly AML)","All"),
                         server=TRUE,selected='All')
    
    updateSelectizeInput(session,"dense_gene",
                         choices=row.names(MCRI_data$E),
                         server=TRUE,selected='AHR')
    
})