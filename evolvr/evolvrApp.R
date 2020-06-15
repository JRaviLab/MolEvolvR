library(shiny)
library(tidyverse)
library(wordcloud2)
library(sunburstR)
library(rmarkdown)
setwd("..")
source("R/cleanup.R")
source("R/summarize.R")
source("R/plotting.R")
source("R/network.R")


#### UI ####
ui <- navbarPage(
  title = "EvolvR",
  id = 'evolvrMenu',
  tabPanel(
    title = "Upload",
    value = "upload",
    fluidRow(

      fileInput(inputId = "fileUpload",
                label = "Choose --- File"
                  ),
      radioButtons(inputId = "fileType", label = "File Type:",
                   choices = c("tsv", "csv")),
      actionButton(inputId = "loadExample", "Example Data"),



      textInput(inputId = "queryInput", label = "Queries:", value = ""),

      verbatimTextOutput("queryVect"),

      tableOutput("uploadedContent")
      #dataTableOutput(outputId = "uploadedContent")
    )
  ),
  source("evolvr/ui/queryDataUI.R")$value,
  source("evolvr/ui/domainArchitectureUI.R")$value,
  source("evolvr/ui/genomicContextUI.R")$value,
  source("evolvr/ui/phylogenyUI.R")$value

)


server <- function(input, output, session)
{

  last_button_val = reactiveVal(0)

  #### Load Data ####
  data <- reactive({
    #req(input$fileUpload)
    if(input$loadExample != last_button_val())
    {
      print(input$loadExample)
      df = read_tsv("data/rawdata_tsv/all_clean.txt")
      last_button_val(last_button_val()+1)
    }
    else
    {
      df <- switch(input$fileType,
                   "tsv" = read_tsv(input$fileUpload$datapath),
                   "csv" = read_csv(input$fileUpload$datapath))
    }

    df

  })

  # Convert query data into vector
  queries <- reactive({
    # Split input by commas or spaces
    unlist(strsplit(input$queryInput, "\\s*,\\s*|\\s+"))
  })

  output$queryVect <- renderPrint({
    queries()
  })

  output$uploadedContent <- renderTable({
    req(input$fileUpload)
    head(data(),n =  5)
  })

  ####  Load Example Data
  observeEvent(input$loadExample,
               {
                 updateTextInput(session, inputId = "queryInput", value =
                                   c("PspA", "Snf7", "PspB","PspC", "PspM", "PspN","DUF3046", "LiaI-LiaF-TM", "Toast-rack",
                                     "Tfu_1009","DUF1700-ahelical", "DUF1707-SHOCT"))

               })


  ### Identify changes to queries and update dropdowns accordingly ###
  observe({
    updateSelectInput(session, "mainSelect", label = "Protein", choices = append("All", queries()))
    updateSelectInput(session, "DASelect", label = "Protein", choices = append("All", queries()))
    updateSelectInput(session, "GCSelect", label = "Protein", choices = append("All", queries()))

  })

  observe({
    if(input$phylo == "Paralog")
    {
      updateSelectInput(session, "PhyloSelect", label = "Protein", choices = queries())
    }
    else
    {
      updateSelectInput(session, "PhyloSelect", label = "Protein", choices = append("All", queries()))
    }
  })

  #### Main Table ####
  output$mainTable = DT::renderDataTable({
    if(input$mainSelect == "All")
    {
      data()
    }
    else
    {
      data() %>% filter(grepl(input$mainSelect, DomArch, ignore.case = T))
    }
  })

  #### Query Heatmap ####
  output$queryHeatmap <- renderPlot({
    req(input$fileUpload)
    req(input$queryInput)
    lineage.Query.plot(data(), queries = queries(), colname = "DomArch", cutoff = 100)
  })


  #### Domain Architecture/ Genomic Context tabs ####
  DA_Prot <- reactive({
    if(input$DASelect == "All")
    {
      data()
    }
    else
    {
      data() %>% filter((grepl(input$DASelect, DomArch, ignore.case = T)))
    }
  })

  GC_Prot <- reactive({
    if(input$GCSelect == "All")
    {
      data()
    }
    else
    {
      data() %>% filter((grepl(input$GCSelect, DomArch, ignore.case = T)))
    }
  })

  ##   #   ###    #    ##
  # Total Counts for DA #
  ##   #   ###    #    ##
  DA_TotalCounts <- reactive({
    prot_tc <- total_counts(DA_Prot(), cutoff = 100, column = "DomArch")
    prot_tc$Lineage = map(prot_tc$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
      unlist()
    prot_tc
  })

  ##   #   ###    #    ##
  # Total Counts for GC #
  ##   #   ###    #    ##
  GC_TotalCounts <- reactive({
    prot_tc <- total_counts(GC_Prot(), cutoff = 100, column = "GenContext")
    prot_tc$Lineage = map(prot_tc$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
      unlist()
    prot_tc
  })



  ###### Heatmap ######
  output$DALinPlot <- renderPlot({
    req(input$fileUpload)
    lineage.DA.plot(DA_Prot(), colname = "DomArch", cutoff = 100, RowsCutoff = F)
  })

  output$GCLinPlot <- renderPlot({
    req(input$fileUpload)
    lineage.DA.plot(GC_Prot(), colname = "GenContext", cutoff = 30, RowsCutoff = F)
  })


  #### Lineage Table ####
  DAlin_count_table <- reactive({
    DA_TotalCounts() %>% group_by(DomArch, totalcount, CumulativePercent) %>%
      summarize(LineageCount = n()) %>%
      select(DomArch, LineageCount, totalcount, CumulativePercent) %>%
      arrange(-totalcount)
  })
  output$DALinTable <- DT::renderDT({
    paged_table(
      {
        DAlin_count_table()
      }
    )
  },
  selection = 'single',
  extensions = c("FixedHeader"),
  options = list(pageLength = 25,
                 scrollX = F,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 0, rightColumns = 0)
  ))


  ## Dialog box that appears when a row is selected
  observeEvent(input$DALinTable_rows_selected, {
    selected_DA = DAlin_count_table()$DomArch[input$DALinTable_rows_selected]
    showModal(modalDialog(
      title = paste("Lineages Of:", selected_DA),
      DT::renderDT({
        DA_TotalCounts() %>% select(DomArch, Lineage, count) %>% arrange(-count) %>%
          filter(DomArch == selected_DA)
        #filter(grepl(regexDAs, DomArch))
      },
      options = list(pageLength = 15,
                     scrollX = T,
                     paging=TRUE,
                     fixedHeader=TRUE,
                     fixedColumns = FALSE)
      ),
      size = "l",
      footer = modalButton(label = "Close"),
      easyClose = T
    ))
  })

  ## GC lin table ##
  GClin_count_table <- reactive({
    GC_TotalCounts() %>% group_by(GenContext, totalcount, CumulativePercent) %>%
      summarize(LineageCount = n()) %>%
      select(GenContext, LineageCount, totalcount, CumulativePercent) %>%
      arrange(-totalcount)
  })
  output$GCLinTable <- DT::renderDT({
    paged_table(
      {
        GClin_count_table()
      }
    )
  },
  selection = 'single',
  extensions = c("FixedHeader"),
  options = list(pageLength = 25,
                 scrollX = F,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 0, rightColumns = 0)
  ))

  ## Dialog box that appears when a row is selected
  observeEvent(input$GCLinTable_rows_selected, {
    selected_GC = GClin_count_table()$GenContext[input$GCLinTable_rows_selected]
    showModal(modalDialog(
      title = paste("Lineages Of:", selected_GC),
      DT::renderDT({
        GC_TotalCounts() %>% select(GenContext, Lineage, count) %>% arrange(-count) %>%
          filter(GenContext == selected_GC)
      },
      options = list(pageLength = 15,
                     scrollX = T,
                     paging=TRUE,
                     fixedHeader=TRUE,
                     fixedColumns = list(leftColumns = 0, rightColumns = 0))
      ),
      size = "l",
      footer = modalButton(label = "Close"),
      easyClose = T
    ))
  })


  #### Upset Plot ####
  output$DAUpsetP <- renderPlot({
    upset.plot(DA_Prot(), colname = "DomArch", cutoff = 100)
  })

  output$GCUpsetP <- renderPlot({
    upset.plot(GC_Prot(), colname = "GenContext", cutoff = 30)
  })

  #### WordCloud ####



  #### Reactive expression determining domain of interest for plotting domain networks
  network_domain_interest <-  reactive({
    if(input$DASelect == "All")
    {
      ".*"
    }
    else
    {
      input$DASelect
    }
  })
  #### Network ####
  output$DANetwork <- renderPlot({
    domain_network(prot = DA_Prot(), column = "DomArch.repeats",
                   domains_of_interest = network_domain_interest(),
                   cutoff = 100)
                   #layout = "layou")
  })




  #### Phyloeny Prot
  phylogeny_prot <- reactive({
    if(input$PhyloSelect == "All")
    {
      data()
    }
    else{
      data() %>% filter(grepl(input$PhyloSelect, DomArch, ignore.case = T))
    }
  })


  #### Sunburst ####
  output$sunburst <- renderSunburst({
    lineage_sunburst(phylogeny_prot(), lineage_column = "Lineage", type = "sunburst", levels = input$levels)
  })


  #### Paralogs ####
  output$ParalogTable <- DT::renderDataTable({
    find_paralogs(phylogeny_prot())
  },extensions = c('FixedColumns'),
  options = list(pageLength = 25,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = FALSE))


}


shinyApp(ui, server,options = options(shiny.maxRequestSize=100*1024^2))


# shiny::runExample("09_upload")