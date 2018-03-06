library(shiny)
library(DT)
options(shiny.maxRequestSize=40*1024^2)

ui <- fluidPage( 
  sidebarLayout(
    sidebarPanel(fileInput(inputId = "tablein", label="Upload your file here"),
                  sliderInput(inputId = "a7_dp", label = "A7 Alelic Depth", min = 0, max = 100, value = 10),
                  sliderInput(inputId = "e5_dp", label = "E5 Alelic Depth", min = 0, max = 100, value = 10),
                  sliderInput(inputId = "rf", label = "Reference frequency diference", min = 0.3, max = 1, value = 0.2), width = 3), 
    mainPanel(DT::dataTableOutput("table_a7e5")) 
  )
) 

#snp_table <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/a7_e5_Eli_annotated.bed", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

server <- function(input, output) {
  
  df <- reactive({
    snp_table <- read.csv(input$tablein$datapath, header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
    is.num <- sapply(snp_table, is.numeric)
    snp_table[is.num] <- lapply(snp_table[is.num], round, 2)
    snp_table <- snp_table[,c(-5,-6,-9,-10)]
    snp_table["A7_DP"] <- unlist(lapply(snp_table$A7_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
    snp_table["E5_DP"] <- unlist(lapply(snp_table$E5_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
    freqs_a7 <- lapply(snp_table$A7_AD, function(x) as.integer(unlist(strsplit(x,","))))
    freqs_e5 <-  lapply(snp_table$E5_AD, function(x) as.integer(unlist(strsplit(x,","))))
    snp_table["A7_RF"] <- unlist(lapply(freqs_a7, function(x) round(x[1]/sum(x), 2)))
    snp_table["E5_RF"] <- unlist(lapply(freqs_e5, function(x) round(x[1]/sum(x), 2)))
    return(snp_table)
  })
  
    output$table_a7e5 <- DT::renderDataTable({
        snp_table <- df()
        expr = snp_table[snp_table$A7_DP >= input$a7_dp & 
                          snp_table$E5_DP >= input$e5_dp &
                          abs(snp_table$A7_RF - snp_table$E5_RF) >= input$rf,]
        })
}

shinyApp(ui = ui, server = server)
