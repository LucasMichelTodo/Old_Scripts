library(shiny)
library(DT)
options(shiny.maxRequestSize=40*1024^2)


calculate_AD <- function(x){
  return(sum(as.integer(unlist(strsplit(x,",")))))
}

ui <- fluidPage( 
  sidebarLayout(
    sidebarPanel(fileInput(inputId = "tablein", label="Upload your file here"),
                 sliderInput(inputId = "dp", label = "Min Alelic Depth", min = 0, max = 200, value = 10),
                 sliderInput(inputId = "rf", label = "Màx. Reference freq. diference", min = 0.6, max = 1, value = 0.6),
                 checkboxGroupInput(
                   "regions", "Mutated region:", 
                   choices = c("ORF", "3region", "5region"),
                   selected = c("ORF", "3region", "5region")),
                 checkboxGroupInput(
                   "var_types", "Mutation consequences:", 
                   choices = c("inframe_insertion", "inframe_deletion", "synonymous_variant", "missense_variant", "stop_gained", "splice_region_variant", "non-coding"),
                   selected =  c("inframe_insertion", "inframe_deletion", "synonymous_variant", "missense_variant", "stop_gained", "splice_region_variant")),
                 checkboxGroupInput(
                   "exprs", "Differential expression:", 
                   choices = "Only differentially expressed",
                   selected = NULL),
                 downloadButton("downloadData", "Download current selection"), width = 2),
    mainPanel(DT::dataTableOutput("table")) 
  )
) 

#snp_table <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/All_SNP_table_anotated_difs_vepannotated_merged_exprs.txt", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)


server <- function(input, output) {
  
  download_name <- reactive(gsub(".txt", "_selection.csv", input$tablein$name))
  
  df <- reactive({
    snp_table <- read.csv(input$tablein$datapath, header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
    
    #Calculate min DP across cols
    min_dp <- c()
    for (i in 1:dim(snp_table)[1]){
      min_dp[i] <- (min(unlist(lapply(snp_table[i,c(6:10)], FUN = calculate_AD))))
    }
    snp_table["min_dp"] <- min_dp
    
    #Calculate max rf diference across cols
    max_rf <- c()
    for (i in 1:dim(snp_table)[1]){
      max_rf[i] <- max(c(dist(c(snp_table[i,11], snp_table[i,12], snp_table[i,13], snp_table[i,14] ,snp_table[i,15]))))
    }
    snp_table["max_rf"] <- max_rf
    
    final_table <- snp_table[snp_table$min_dp >= input$`dp` &
                              snp_table$max_rf >= input$rf &
                              snp_table$Region %in% input$regions &
                              snp_table$Var_Type %in% input$var_types,-1]
    
    if (!is.null(input$exprs)){
      final_table <- final_table[final_table$Exp_dif == TRUE,]
    }
                              #strsplit(snp_table$Var_Type, ",")[1] %in% input$var_types,]
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


