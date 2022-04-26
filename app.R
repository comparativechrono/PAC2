library(shiny)
#load dataset
expression <- read.table("data/PAC2_RNA-seq.tsv", header = TRUE)


# Use a fluid Bootstrap layout
ui <- fluidPage(    
  
  # Give the page a title
  titlePanel(h1("Photoperiodic gene expression in zebrafish PAC2 cells", align = "center"), windowTitle = "whitmorelab.shinyapps.io"),
  
  # Generate a row with a sidebar
  sidebarLayout(      
    
    # Define the sidebar with one input
    sidebarPanel(
      
     # SET DROP DOWN FOR PRIMARY TRANSCRIPT
      
      htmlOutput("gene_selector"),
      htmlOutput("transcript_selector"),
      
      # SELECT CONDTIIONS FOR PRIMARY TRANSCRIPT,
      
      helpText("Pick a primary timecourse to plot (red)"),
      
      selectInput("FILTER", label = "Select a photoperiod", 
                  choices = list("Long days" = "LD", "Short days" = "SD"), 
                  selected = "LD"),
      
      selectInput("GENOTYPE", label = "Select a genotype", 
                  choices = list("Wild-type PAC2 cells" = "WT", "CLOCK dn cells" = "CLOCK"), 
                  selected = 1),
      helpText("Pick an optional secondary timecourse to plot (blue)"),

      # SELECT CONDTIIONS FOR SECONDARY TRANSCRIPT
      
      selectInput("FILTER2", label = "Select a photoperiod", 
                  choices = list("Long days" = "LD", "Short days" = "SD"), 
                  selected = "SD"),
      
      selectInput("GENOTYPE2", label = "Select a genotype", 
                  choices = list("Wild-type PAC2 cells" = "WT", "CLOCK dn cells" = "CLOCK"), 
                  selected = 1)
               ),
   

    
    
    # Create a spot for the plot
    mainPanel(
      plotOutput("genePlot") ,
      
      hr(),
      downloadButton("downloadData", "Download"),
      hr(),
      h3("Funding"),
      p("This data was generated through Leverhulme Trust Research Prject Grant 2017: Can individual cells in culture tell the time of year?"),
      img(src = "leverhulme.png", height = 50),
      h3("Acknowledgement"),
      p("Please cite this webpage if you make use of this data")
      
    )
    
  )
)

# Define a server for the Shiny app
server <- function(input, output) {
  

  output$gene_selector = renderUI({
    # set gene for primary plot
    #make vector with gene names
    selectInput("SEARCH", label = "Start typing your gene here (lowercase)", 
                choices = unique(expression[expression$Gene_name,3]),
                  selected = "cry1aa")
  })
  
  output$transcript_selector = renderUI({ 
    #SEARCH FOR TRANSCRIPTS on primary plot
    available_transcripts <- expression[expression$Gene_name == input$SEARCH & expression$FILTER == "LD" & expression$GENOTYPE == "WT",2]
    
    selectInput("GOI", label = "Select Ensembl Transcript ID:", 
                choices = unique(available_transcripts),
                selected = "ENSDART00000130692")
  })  
  
  
  # Fill in the spot we created for a plot
  output$genePlot <- renderPlot({
  
    # define our exis
    primary <- expression[which(expression$Transcript_ID == input$GOI & expression$FILTER == input$FILTER & expression$GENOTYPE == input$GENOTYPE),c(6:13)]  
    secondary <- expression[which(expression$Transcript_ID == input$GOI & expression$FILTER == input$FILTER2 & expression$GENOTYPE == input$GENOTYPE2),c(6:13)]
    SEMprimary <- expression[(expression$Transcript_ID==input$GOI & expression$FILTER==input$FILTER & expression$GENOTYPE==input$GENOTYPE),c(14:21)]
    SEMsecondary <-  expression[(expression$Transcript_ID==input$GOI & expression$FILTER==input$FILTER2 & expression$GENOTYPE==input$GENOTYPE2),c(14:21)]
    
    
    # Render a plot with secondary on
    plot(x=c(3,9,15,21,27,33,39,45), y=secondary,
            main=input$GOI,
            ylab="Gene expression (Read counts)",
            xlab="Time (h)",
         type="b",
         lwd=2,
         col="blue",
         pch = 21, bg = "blue", cex = 1.5,
         ylim = c(min(primary,secondary)/1.4, max(primary,secondary)*1.4),
         xaxt="n"
    )
    xtick<-seq(0, 48, by=8)
    axis(side=1, at=xtick)
    
    #add error bars to secondary
    arrows(x0=c(3,9,15,21,27,33,39,45), 
           y0=as.numeric(secondary + SEMsecondary),
           x1=c(3,9,15,21,27,33,39,45), 
           y1=as.numeric(secondary - SEMsecondary), code=3, angle=90, length=0.03, col="blue", lwd=1)
   
    
    # render a plot of primary choice
    
    lines(x=c(3,9,15,21,27,33,39,45), y=primary, 
          type="b",
          col="red",
          pch = 21, bg = "red", cex = 1.5,
          lwd=2)
    #add error bars to primary
    arrows(x0=c(3,9,15,21,27,33,39,45), 
           y0=as.numeric(primary + SEMprimary),
           x1=c(3,9,15,21,27,33,39,45), 
           y1=as.numeric(primary - SEMprimary), code=3, angle=90, length=0.03, col="red", lwd=1)
    
    
    }, height = 400, width = 600)
  
  #define my download
  output$downloadData <- downloadHandler(
    filename = "my_graphs.png",
    content = function(file) {
 png(file, height = 400, width = 600)
      
      primary <- expression[which(expression$Transcript_ID == input$GOI & expression$FILTER == input$FILTER & expression$GENOTYPE == input$GENOTYPE),c(6:13)]  
      secondary <- expression[which(expression$Transcript_ID == input$GOI & expression$FILTER == input$FILTER2 & expression$GENOTYPE == input$GENOTYPE2),c(6:13)]
      SEMprimary <- expression[(expression$Transcript_ID==input$GOI & expression$FILTER==input$FILTER & expression$GENOTYPE==input$GENOTYPE),c(14:21)]
      SEMsecondary <-  expression[(expression$Transcript_ID==input$GOI & expression$FILTER==input$FILTER2 & expression$GENOTYPE==input$GENOTYPE2),c(14:21)]
      
      
      plot(x=c(3,9,15,21,27,33,39,45), y=secondary,
           main=input$GOI,
           ylab="Gene expression (Read counts)",
           xlab="Time (h)",
           type="b",
           lwd=2,
           col="blue",
           pch = 21, bg = "blue", cex = 1.5,
           ylim = c(min(primary,secondary)/1.4, max(primary,secondary)*1.2),
           xaxt="n"
      )
      xtick<-seq(0, 48, by=8)
      axis(side=1, at=xtick)
      #add error bars to secondary
      arrows(x0=c(3,9,15,21,27,33,39,45), 
             y0=as.numeric(secondary + SEMsecondary),
             x1=c(3,9,15,21,27,33,39,45), 
             y1=as.numeric(secondary - SEMsecondary), code=3, angle=90, length=0.03, col="blue", lwd=1)
      
      
      # render a plot of primary choice
      
      lines(x=c(3,9,15,21,27,33,39,45), y=primary, 
            type="b",
            col="red",
            pch = 21, bg = "red", cex = 1.5,
            lwd=2)
      #add error bars to primary
      arrows(x0=c(3,9,15,21,27,33,39,45), 
             y0=as.numeric(primary + SEMprimary),
             x1=c(3,9,15,21,27,33,39,45), 
             y1=as.numeric(primary - SEMprimary), code=3, angle=90, length=0.03, col="red", lwd=1)
      
      dev.off()
     
    },
 contentType = 'image/png'
  )

  
}
# Run the app ----
shinyApp(ui = ui, server = server)