## LIBRARIES
library(shiny)
library(shinyjs)
library(bslib)
library(tidyverse)


## FUNCTIONS
# Loading in external functions
source("utils/utils.R")


## APP CONFIGURATION (hard-coded example, to be removed)
# Reading in data
longTPM <- readRDS("data/tpms.RDS")
res.all.sh <- readRDS("data/all_results_shrinkedLFC_curated.RDS")
all.symbols <- res.all.sh$symbol %>% unique()


## UI
ui <- page_sidebar(
  title = "Differential Expression Analysis Visual Representer",
  # theme = bs_theme_update(theme = bs_current_theme(), bg = "#ffffff", fg = "#000", primary = "#00891a", secondary = "#000"),
  sidebar = sidebar(
    
    useShinyjs(),
    h5("Plot settings:"),
    
    selectizeInput(
      inputId = "geneID", 
      label = "Symbol searcher", 
      multiple = FALSE,
      choices = NULL,
      options = list(
        create = FALSE,
        placeholder = "Type in symbol",
        onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
        onType = I("function (str) {if (str === \"\") {this.close();}}"),
        onItemAdd = I("function() {this.close();}")
      )
    ),
    
    # Section for customizing plot download
    checkboxInput("show_comparisons", "LFC comparisons", TRUE),
    checkboxInput("gray_scale", "Gray scale", TRUE),
    checkboxInput("theme_minimal", "Minimal theme", FALSE),
    sliderInput("plot_dims", "Plot dimensions (n x n)", min = 1, max = 5, value = 3),
    radioButtons(
      inputId = "file_format",
      label = "File format",
      choices = list("PDF", "PNG", "JPEG"),
      selected = "PDF"
    ),
    
    # Download button
    downloadButton(outputId = "save_plot", label = "Save plot"),
  ),

  tabsetPanel(
    id = "main_tab_panel",
    tabPanel("Results",
      br(),
      card(
        full_screen = TRUE,
        card_header("TPM plot"),
        plotOutput("tpm"),
      ),
      card(
        card_header("LFC values"),
        card_body("The DEA results for the selected gene will appear here."),
        tableOutput("lfc_results")
      )
    ),
    tabPanel("About",
      br(),
      p("Hi! This is a little web app developed using Shiny in R to aid in the representation of the Differential Expression Analysis results (DEA)."),
      p("For now, it is hardcoded to display the results of my master's thesis at the Barcelona Supercomputing Center (BSC) and Universidad Autónoma de Madrid (UAM). I hope this little web app can help inspire others to start with shiny coding!"),
      tagList("You can check me at", a("LinkedIn!", href = "https://www.linkedin.com/in/xavier-benedicto/")),
      br()
    )
  )
)


## SERVER
server <- function(input, output, session) {
  
  # Increasing max requests size to 10 megabytes
  options(shiny.maxRequestSize = 10*1024^2)
  
  # Generating needed data
  longTPM.filt <- reactive({
    req(input$geneID)
    longTPM %>% filter(symbol == input$geneID) %>% as.data.frame()
  })
  LFCvalues <- reactive({
    req(input$geneID)
    res.all.sh %>% filter(symbol == input$geneID) %>% select(-c(baseMean, lfcSE))
  })
  
  # Adding all possible gene options at this step to avoid freezing
  observe({
    updateSelectizeInput(session, "geneID", selected = "MAPK3", choices = all.symbols, server = TRUE)
  })
  
  output$tpm <- renderPlot({
    TPMplot(longTPM.filt(), input$geneID, LFCvalues(), gray_color = input$gray_scale, minimal_theme = input$theme_minimal, show_comparisons = input$show_comparisons)
  }, res = 200)
  
  output$lfc_results <- renderTable(
    LFCvalues(),
    striped = FALSE,
    hover = TRUE,
    digits = 8,
    width = "100%"
  )
  
  output$save_plot <- downloadHandler(
    filename = function() {
      
      # Reactive to the file format
      selected_format <- reactive({
        switch(
          input$file_format,
          "PDF" = ".pdf",
          "PNG" = ".png",
          "JPEG" = ".jpeg"
        )
      })
      
      # Adding selected file format
      paste0("tpm-", input$geneID, "-", Sys.Date(), selected_format())
    },
    content = function(file) {
      
      # Generating and saving plot
      gg <- TPMplot(longTPM.filt(), input$geneID, LFCvalues(), gray_color = input$gray_scale, minimal_theme = input$theme_minimal, show_comparisons = input$show_comparisons)
      ggsave(file, plot = gg, width = input$plot_dims, height = input$plot_dims)
    }
  )
}


## APP
shinyApp(ui, server)
