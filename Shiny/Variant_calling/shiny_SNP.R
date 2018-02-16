library(shiny)
library(DT)
options(shiny.maxRequestSize=40*1024^2)

ui <- shinyUI(navbarPage("SNPs and Indels",
                   tabPanel("Whole table", 
                            fileInput(inputId = "tablein", label="Upload your file here"),
                            DT::dataTableOutput("table")),
                   navbarMenu("Contrasts",
                              tabPanel("A7 vs E5", sidebarLayout( 
                                                          sidebarPanel(sliderInput(inputId = "a7_dp", label = "A7 Alelic Depth", min = 0, max = 100, value = 10),
                                                                       sliderInput(inputId = "e5_dp", label = "E5 Alelic Depth", min = 0, max = 100, value = 10),
                                                                       sliderInput(inputId = "rf", label = "Reference frequency diference", min = 0, max = 1, value = 0.2), width = 3), 
                                                          mainPanel(DT::dataTableOutput("table_a7e5"))) 
                                
                                      ))
))

server <- function(input, output) {
  
  df <- reactive({
    snp_table <- read.csv(input$tablein$datapath, header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
    snp_table <- snp_table[complete.cases(snp_table),]
    is.num <- sapply(snp_table, is.numeric)
    snp_table[is.num] <- lapply(snp_table[is.num], round, 2)
    colnames(snp_table) <- c("Crom", "Pos", "Ref", "Alt", 
                             "1.2B_GT", "1.2B_AD", "1.2B_GQ", 
                             "10G_GT", "10G_AD", "10G_GQ", 
                             "A7_GT", "A7_AD", "A7_GQ", 
                             "C2_GT", "C2_AD", "C2_GQ",
                             "E5_GT", "E5_AD", "E5_GQ")
    
    return(snp_table)
    })
    

#snp_table <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/all_unifiedgenotyper_Eli_table.txt", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
# is.num <- sapply(snp_table, is.numeric)
# snp_table[is.num] <- lapply(snp_table[is.num], round, 2)
# colnames(snp_table) <- c("Crom", "Pos", "Ref", "Alt",
#                          "1.2B_GT", "1.2B_AD", "1.2B_GQ",
#                          "10G_GT", "10G_AD", "10G_GQ",
#                          "A7_GT", "A7_AD", "A7_GQ",
#                          "C2_GT", "C2_AD", "C2_GQ",
#                          "E5_GT", "E5_AD", "E5_GQ")

  
  output$table <- DT::renderDataTable({expr = df()}, server = TRUE)
  
  output$table_a7e5 <- DT::renderDataTable({
    snp_table <- df()
    snp_table <- snp_table[,c("Crom", "Pos", "Ref", "Alt", "A7_GT", "E5_GT", "A7_AD", "E5_AD", "A7_GQ", "E5_GQ")]
    snp_table["A7_DP"] <- unlist(lapply(snp_table$A7_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
    snp_table["E5_DP"] <- unlist(lapply(snp_table$E5_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
    
    freqs_a7 <- lapply(snp_table$A7_AD, function(x) as.integer(unlist(strsplit(x,","))))
    freqs_e5 <-  lapply(snp_table$E5_AD, function(x) as.integer(unlist(strsplit(x,","))))
    snp_table["A7_RF"] <- unlist(lapply(freqs_a7, function(x) round(x[1]/sum(x), 2)))
    snp_table["E5_RF"] <- unlist(lapply(freqs_e5, function(x) round(x[1]/sum(x), 2)))
    
    expr = snp_table[snp_table$A7_DP >= input$a7_dp & 
                      snp_table$E5_DP >= input$e5_dp &
                      abs(snp_table$A7_RF - snp_table$E5_RF) >= input$rf,]
    })
}

shinyApp(ui = ui, server = server)
