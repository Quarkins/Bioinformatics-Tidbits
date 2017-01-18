library(shiny)
library(limma)
library(RColorBrewer)
library(ggplot2)
library(tidyr)

#Load in the data
MCRI_data <- readRDS("Data/MCRI.Rds")
MCRI_target <- readRDS("Data/MCRI_target.Rds")

shinyServer(function(input, output, session) {
    
    #First reactively subset dataset  and targets based on Cell type
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
    
    dt <- reactive({
        if(input$type != "All"){
            selected = MCRI_target[MCRI_target$Cell_Type == input$type,]
            selected
        }
        else{
            MCRI_target
        }
    })
    
    #Now reactively make the dataset ggplotable!
    dg <- reactive({
        if(length(input$gene) < 1){}
        if(length(input$gene) == 1){
            tidier = data.frame(Sample = colnames(dd()),Genes=input$gene,logFPKM=dd()[input$gene,])
        }
        else{
            ngenes = length(input$gene)
            messy = as.data.frame(t(dd()[input$gene,]))
            messy$Sample = row.names(messy) #Get sample row names as a variable
            tidier = messy %>% gather(Genes,logFPKM,c(1:ngenes))
            tidier
        }
    })
    
    
    output$expression <- renderPlot({
        
        if(length(input$gene) < 1){}   
        
        else{
            q <- ggplot(dg(),aes(x=factor(Genes),y=logFPKM,fill=factor(Genes),col=factor(Genes)))
            q + geom_jitter() + theme_bw() + 
                theme(axis.title.x=element_blank(),legend.position="none")
            
        }
     
    })
    
 
    output$boxplot <- renderPlot({
        if(length(input$gene) < 1){}
        else{
            q <- ggplot(dg(),aes(x=factor(Genes),y=logFPKM,fill=factor(Genes),col=factor(Genes)))
            q + geom_boxplot(alpha=0.8) + theme_bw() + 
                theme(axis.title.x=element_blank(),legend.position="none") 
        }
    })
    
    output$density <- renderPlot({
        if(length(input$gene) < 1){}
        else{
            q <- ggplot(dg())
            q + geom_density(aes(logFPKM,col=factor(Genes),fill=factor(Genes)),alpha=0.1) + theme_bw() + theme(legend.position="none")
        }
    })
    
    output$outliers <- renderTable({
        #Get outliers
        sub = t(dd()[input$dense_gene,,drop=FALSE])
        bp = boxplot(sub,plot=FALSE)
        tab = data.frame(outlier_sample = rownames(sub)[match(bp$out,sub)], logFPKM = bp$out)
        tab
    })
    
    
    output$MDS <- renderPlot({
        pal = brewer.pal(5,"Set1")
        plotMDS(MCRI_data,gene.selection="common",col=pal[as.factor(MCRI_target$Classified)],cex=0.8,pch=16)
        legend('topleft',legend=levels(as.factor(MCRI_target$Classified)),col=pal,cex=1,bty = "n",pch=16)
        
    })
    
    output$classified <- renderTable({
        if(length(input$gene) < 1){ dt()}
        else{
            if(length(input$gene) == 1){ 
                combo_breaker = dt()
                combo_breaker[,input$gene] = dd()[input$gene,]
            }
            else{ combo_breaker = cbind(dt(),t(dd()[input$gene,])) }
            combo_breaker
        }
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