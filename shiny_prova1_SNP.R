library(shiny)

snp_table <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/a7_e5_unifiedgeno_eli.csv", 
            header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

snp_table["A7_DP"] <- unlist(lapply(snp_table$A7_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))
snp_table["E5_DP"] <- unlist(lapply(snp_table$E5_AD, function(x) sum(as.integer(unlist(strsplit(x,","))))))

ui <- fluidPage(
  sliderInput(inputId = "a7_dp", label = "A7 Alelic Depth", min = 0, max = 100, value = 10),
  sliderInput(inputId = "e5_dp", label = "E5 Alelic Depth", min = 0, max = 100, value = 10),
  tableOutput("table")
)

server <- function(input, output) {
  output$table <- renderTable({
    expr = snp_table[snp_table$A7_DP >= input$a7_dp & snp_table$E5_DP >= input$e5_dp,]}, digits = 2, 
    caption = "SNPs and Indels", 
    caption.placement = "top")
}

shinyApp(ui = ui, server = server)
