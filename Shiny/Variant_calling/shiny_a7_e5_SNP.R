library(shiny)
library(DT)
options(shiny.maxRequestSize=40*1024^2)

ui <- fluidPage( 
  sidebarLayout(
    sidebarPanel(fileInput(inputId = "tablein", label="Upload your file here"),
                  sliderInput(inputId = "1_dp", label = "1 Alelic Depth", min = 0, max = 100, value = 10),
                  sliderInput(inputId = "2_dp", label = "2 Alelic Depth", min = 0, max = 100, value = 10),
                  sliderInput(inputId = "rf", label = "Reference frequency diference", min = 0.3, max = 1, value = 0.2),
                  checkboxGroupInput(
                    "regions", "Display mutations in:", 
                    choices = c("ORF", "3region", "5region"),
                    selected = c("ORF", "3region", "5region")),
                 downloadButton("downloadData", "Download current selection"),
                 width = 3),
    mainPanel(DT::dataTableOutput("table")) 
  )
) 

#snp_table <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/a7_e5_annotated.bed", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

server <- function(input, output) {
  
  download_name <- reactive(gsub(".bed", "_selection.csv", input$tablein$name))
  
  df <- reactive({
    snp_table <- read.csv(input$tablein$datapath, header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
    is.num <- sapply(snp_table, is.numeric)
    snp_table[is.num] <- lapply(snp_table[is.num], round, 2)
    snp_table <- snp_table[,c(-5,-6,-9,-10)]
    snp_table["X1_DP"] <- unlist(lapply(snp_table$X1_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
    snp_table["X2_DP"] <- unlist(lapply(snp_table$X2_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
    freqs_1 <- lapply(snp_table$X1_AD, function(x) as.integer(unlist(strsplit(x,","))))
    freqs_2 <-  lapply(snp_table$X2_AD, function(x) as.integer(unlist(strsplit(x,","))))
    snp_table["X1_RF"] <- unlist(lapply(freqs_1, function(x) round(x[1]/sum(x), 2)))
    snp_table["X2_RF"] <- unlist(lapply(freqs_2, function(x) round(x[1]/sum(x), 2)))
    final_table <- snp_table[snp_table$X1_DP >= input$`1_dp` & 
                snp_table$X2_DP >= input$`2_dp` &
                abs(snp_table$X1_RF - snp_table$X2_RF) >= input$rf &
                snp_table$Region %in% input$regions,]
    return(final_table)
  })
  
  output$table <- DT::renderDataTable(expr = df(), rownames = FALSE, options = list(pageLength = 50))
  
  output$downloadData <- downloadHandler(
    filename = function(){download_name()},
    content = function(file) {
      write.csv(df(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
