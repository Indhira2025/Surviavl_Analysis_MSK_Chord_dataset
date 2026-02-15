library(shiny)
library(shinydashboard)
library(data.table)
library(ggplot2)
library(maftools)
library(pheatmap)
library(circlize)
library(survival)
library(survminer)
library(DT)
library(dplyr)
library(plotly)
#install.packages("heatmaply")
library(heatmaply)
# --- Load Data ---
clin_patient <- fread("Data/data_clinical_patient.txt", skip = 4)
clin_sample  <- fread("Data/data_clinical_sample.txt", skip = 4)
muta <- fread("Data/data_mutations.txt")
cna  <- fread("Data/data_cna.txt")
sv   <- fread("Data/data_sv.txt")

# Ensure survival variable
clin_patient$OS_STATUS_BIN <- ifelse(grepl("DECEASED", clin_patient$OS_STATUS), 1, 0)
clin_patient <- as.data.table(clin_patient)
clin_patient <- as.data.table(clin_patient)

# Merge by PATIENT_ID
merged_data <- merge(clin_sample, clin_patient, by = "PATIENT_ID", all.x = TRUE)
 

# --- UI ---
ui <- dashboardPage(
  dashboardHeader(title = "MSK-IMPACT Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      selectInput("cancer_type",
                  "Select Cancer Type:",
                  choices = sort(unique(clin_sample$CANCER_TYPE)),
                  selected = unique(clin_sample$CANCER_TYPE)[1]),
      
      menuItem("Overview", tabName = "overview"),
      menuItem ("Metastatic Sites", tabName = "metastatic_sites"),
      menuItem("Oncoplot", tabName = "oncoplot"),
      menuItem("CNA", tabName = "cna"),
      menuItem("Fusions", tabName = "fusion"),
      menuItem("TMB(Tumor Mutation Burden)", tabName = "tmb"),
      menuItem("Survival", tabName = "survival")
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "overview",
              fluidRow(
                valueBoxOutput("totalPatients"),
                valueBoxOutput("totalSamples"),
                valueBoxOutput("totalMutations")
              ),
              fluidRow(
                box(width = 12, DTOutput("cancerTable"))
              )
      ),
      tabItem(tabName = "metastatic_sites",
              box(width = 12, plotlyOutput("metastaticPlot", height = 700))
      ),
      
      tabItem(tabName = "oncoplot",
              box(width = 12, plotOutput("oncoplotPlot", height = 700),
                  footer = "Shows top mutated genes across selected samples; colors indicate mutation type and height reflects mutation frequency")
      ),
      tabItem(tabName = "cna",
              box(width = 12, plotOutput("cnaPlot", height = 800),
                  footer = "Colors indicate copy number gains (pink) and losses (green) across selected samples")
      ),
      tabItem(tabName = "fusion", 
              box(width = 12, 
                  title = "Gene Fusion Network",
                  plotOutput("fusionPlot", height = 800),
                  footer = "Interpretation: Ribbons represent fusion events between genes. 
                Thickness indicates frequency."
              )
      ),
      tabItem(tabName = "tmb",
              box(width = 12, plotOutput("tmbPlot", height = 500), footer = "TMB measures mutation density (mut/Mb). High TMB (≥10 mut/Mb) often 
         predicts better response to immune checkpoint inhibitors by increasing neoantigen presentation.")
      ),
      tabItem(tabName = "survival",
              box(width = 12, plotOutput("survivalPlot", height = 600))
      )
    )
  )
)

# --- Server ---
server <- function(input, output) {
  
  # --- Reactive filtered samples for selected cancer type ---
  filtered_samples <- reactive({
    merged_data[CANCER_TYPE == input$cancer_type]
  })
  
  # --- Overview Value Boxes ---
  output$totalPatients <- renderValueBox({
    n_pat <- length(unique(filtered_samples()$PATIENT_ID))
    valueBox(n_pat, "Patients", icon = icon("users"), color = "blue")
  })
  
  output$totalSamples <- renderValueBox({
    n_samp <- nrow(filtered_samples())
    valueBox(n_samp, "Samples", icon = icon("vial"), color = "green")
  })
  
  output$totalMutations <- renderValueBox({
    samp_ids <- filtered_samples()$SAMPLE_ID
    n_mut <- muta[Tumor_Sample_Barcode %in% samp_ids, .N]
    valueBox(n_mut, "Mutations", icon = icon("dna"), color = "red")
  })
  
  output$cancerTable <- renderDT({
    datatable(filtered_samples(), options = list(pageLength = 10, scrollX = TRUE))
  })
  # Metastatic Site
  
  output$metastaticPlot <- renderPlotly({
    meta_sub <- merged_data[SAMPLE_ID %in% filtered_samples()$SAMPLE_ID &
                              CANCER_TYPE == input$cancer_type, ]
    
    site_counts <- meta_sub[!is.na(METASTATIC_SITE) & METASTATIC_SITE != "Not Applicable", 
                            .N, by = METASTATIC_SITE]
    
    setorder(site_counts, -N)
    
    if(nrow(site_counts) == 0) return(NULL)
    
    site_counts[, METASTATIC_SITE := factor(METASTATIC_SITE, levels = site_counts$METASTATIC_SITE)]
    
    plot_ly(site_counts, x = ~METASTATIC_SITE, y = ~N, type = 'bar',
            marker = list(color = 'steelblue')) %>%
      layout(title = paste("Metastatic sites for", input$cancer_type),
             xaxis = list(title = "Site"), yaxis = list(title = "Number of Samples"))
  })
  
  
  
  
  
  # --- Oncoplot ---
  output$oncoplotPlot <- renderPlot({
    samp_ids <- filtered_samples()$SAMPLE_ID
    muta_sub <- muta[Tumor_Sample_Barcode %in% samp_ids]
    if(nrow(muta_sub) == 0) return(NULL)
    maf_obj <- read.maf(muta_sub)
    oncoplot(maf_obj, top = 20)
  })
  
  # --- CNA Heatmap ---
  output$cnaPlot <- renderPlot({
    
    # Get selected samples
    samp_ids <- filtered_samples()$SAMPLE_ID
    samp_ids <- samp_ids[samp_ids %in% colnames(cna)]
    if(length(samp_ids) == 0) return(NULL)
    
    # Subset CNA matrix
    cna_sub <- cna[, c("Hugo_Symbol", samp_ids), with = FALSE]
    mat <- as.matrix(cna_sub[, -1, with = FALSE])
    rownames(mat) <- cna_sub$Hugo_Symbol
    
    if(nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
    
    n_samples <- ncol(mat)
    
    # Calculate frequency of amplifications and deletions
    amp_freq <- rowSums(mat > 0, na.rm = TRUE) / n_samples * 100
    del_freq <- rowSums(mat < 0, na.rm = TRUE) / n_samples * 100
    
    # Total alteration frequency
    total_freq <- amp_freq + del_freq
    
    # Remove genes with 0 alterations
    keep <- total_freq > 0
    if(sum(keep) == 0) return(NULL)
    
    amp_freq <- amp_freq[keep]
    del_freq <- del_freq[keep]
    total_freq <- total_freq[keep]
    
    # Select top 20 most frequently altered genes
    top_n <- min(20, length(total_freq))
    top_genes <- names(sort(total_freq, decreasing = TRUE))[1:top_n]
    
    plot_df <- data.frame(
      Gene = factor(top_genes, levels = rev(top_genes)),
      Amplification = amp_freq[top_genes],
      Deletion = -del_freq[top_genes]   # negative for plotting downward
    )
    
    # Plot
    ggplot(plot_df, aes(x = Gene)) +
      geom_bar(aes(y = Amplification), stat = "identity", fill = "lightpink") +
      geom_bar(aes(y = Deletion), stat = "identity", fill = "lightgreen") +
      coord_flip() +
      theme_minimal() +
      ylab("Percentage of Samples (%)") +
      xlab("Gene") +
      ggtitle(paste("Top", top_n, "CNA Frequencies for", input$cancer_type))
  })
  
  
  # --- Fusion Chord Diagram ---
  output$fusionPlot <- renderPlot({
    # 1. Filter and find top pairs
    sv_sub <- sv[Sample_Id %in% filtered_samples()$SAMPLE_ID]
    if(nrow(sv_sub) == 0) return(NULL)
    
    top_pairs <- sv_sub[, .N, by = .(Site1_Hugo_Symbol, Site2_Hugo_Symbol)]
    
    # 2. DATA CLEANING (The likely cause of the 'times' error)
    # - Remove any rows where either gene name is missing or empty string
    # - Remove any rows where N is 0 or NA
    # - Remove "self-fusions" (Site1 == Site2) if they exist, as they can crash chordDiagram
    top_pairs <- top_pairs[
      !is.na(Site1_Hugo_Symbol) & Site1_Hugo_Symbol != "" &
        !is.na(Site2_Hugo_Symbol) & Site2_Hugo_Symbol != "" &
        !is.na(N) & N > 0 &
        Site1_Hugo_Symbol != Site2_Hugo_Symbol
    ]
    
    if(nrow(top_pairs) == 0) return(NULL)
    
    # 3. Take top 20 and convert to STRICT base data.frame
    top_pairs <- top_pairs[order(-N)][1:min(20, .N)]
    
    plot_df <- data.frame(
      from  = as.character(top_pairs$Site1_Hugo_Symbol),
      to    = as.character(top_pairs$Site2_Hugo_Symbol),
      value = as.numeric(top_pairs$N),
      stringsAsFactors = FALSE
    )
    
    # 4. RESET AND PLOT
    # circos.clear() is mandatory in Shiny to prevent the 'times' error 
    # when switching between different cancer types
    circlize::circos.clear() 
    
    # Wrap in tryCatch so the whole app doesn't crash if one specific plot fails
    tryCatch({
      circlize::chordDiagram(
        x = plot_df,
        transparency = 0.5,
        annotationTrack = c("name", "grid")
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Could not render plot with current data selection")
    })
  })
  
  
  # --- Tumor Mutation Burden ---
  output$tmbPlot <- renderPlot({
    data <- filtered_samples()
    if(nrow(data) == 0) return(NULL)
    
    ggplot(data, aes(y = TMB_NONSYNONYMOUS)) +
      geom_boxplot(fill = "steelblue") +
      theme_bw() +
      ylab("TMB (Non-synonymous)") +
      ggtitle(paste("TMB for", input$cancer_type))
  })
  
  # --- Survival Plot ---
  output$survivalPlot <- renderPlot({
    datasurv <- merged_data[merged_data$CANCER_TYPE == input$cancer_type, ]  
    
    if(nrow(datasurv) == 0) return(NULL)
    
    fit <- do.call(survfit, list(
      formula = Surv(as.numeric(OS_MONTHS), OS_STATUS_BIN) ~ 1, 
      data = as.data.frame(datasurv)
    ))
    
    ggsurvplot(fit,
               risk.table = TRUE,
               conf.int = TRUE,
               xlab = "Months",
               ylab = "Survival Probability",
               title = paste("Overall Survival for", input$cancer_type))$plot
  })
}



# --- Run App ---
shinyApp(ui, server)

