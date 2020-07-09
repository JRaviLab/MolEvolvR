library(shiny)
library(shinyjs)
library(tidyverse)
library(wordcloud2)
library(sunburstR)
library(rmarkdown)
library(shinyBS)
library(shinyauthr)
setwd("..")
source("R/cleanup.R")
source("R/summarize.R")
source("R/plotting.R")
source("R/network.R")
source("R/pre-msa-tree.R")
source("scripts/tree.R")
source("evolvr/components.R")
conflicted::conflict_prefer("strsplit", "base")

###
#### Users ####
###
user_base <- data.frame(
  user = c("pspevolution"),
  password = c("cpathogeno2019"),
  permissions = c("admin"),
  name = c("User One"),
  stringsAsFactors = FALSE,
  row.names = NULL
)



example_data = read_tsv("data/rawdata_tsv/all_clean.txt")



#### UI ####
ui <- tagList(
  shinyjs::useShinyjs(),
  tags$head(includeScript("evolvr/logout-button.js")),
  includeCSS("evolvr/styles.css"),
  includeCSS("evolvr/components.css"),
  # tags$head(includeScript("evolvr/splash-button.js")),
  navbarPage(
    title = actionLink(inputId = "homeButton",
                       style ="text-decoration:none; color:white;"
                       , tags$div(class= "zoom",  "EvolvR")),
    id = 'evolvrMenu',
    inverse = T,
    source("evolvr/ui/splashPageUI.R")$value,
    source("evolvr/ui/uploadUI.R")$value,
    source("evolvr/ui/queryDataUI.R")$value,
    source("evolvr/ui/domainArchitectureUI.R")$value,
    source("evolvr/ui/genomicContextUI.R")$value,
    source("evolvr/ui/phylogenyUI.R")$value,
    bsModal(id = "splashModal", title = "Splash Page", trigger = "splashBtn", size = "large")
  )
)


server <- function(input, output, session)
{

  # call the logout module with reactive trigger to hide/show
  logout_init <- callModule(shinyauthr::logout,
                            id = "logout",
                            active = reactive(credentials()$user_auth))

  # call login module supplying data frame, user and password cols
  # and reactive trigger
  credentials <- callModule(shinyauthr::login,
                            id = "login",
                            data = user_base,
                            user_col = user,
                            pwd_col = password,
                            log_out = reactive(logout_init()))

  # pulls out the user information returned from login module
  user_data <- reactive({credentials()$info})

  output$user_table <- renderTable({
    # use req to only render results when credentials()$user_auth is TRUE
    req(credentials()$user_auth)
    user_data()
  })



  ## Change title tab name (ie: chrome tab name) to EvolvR
  shinyjs::runjs('document.title = "EvolvR"')

  ## Redirect click on the app title to the upload-data page
  observeEvent(input$homeButton, updateNavbarPage(session, "evolvrMenu" ,"home"))

  last_button_val = reactiveVal(0)

  #### Load Data ####
  observeEvent(input$fileUpload,
               {
                 req(credentials()$user_auth)
                 print("Upload")
                 df <- switch(input$fileType,
                              "tsv" = read_tsv(input$fileUpload$datapath),
                              "csv" = read_csv(input$fileUpload$datapath))
                 data(df)
               })


  data <- reactiveVal(data.frame())

  fasta_filepath <- reactiveVal("tmp.fasta")
  aligned_fasta_filepath <- reactiveVal("")


  DACutoff <- reactive({
    input$DA_Cutoff
  })

  GCCutoff <- reactive({
    input$GC_Cutoff
  })

  # Convert query data into vector
  queries <- reactive({
    req(typeof(input$queryInput) == "character")

    # Split input by commas or spaces
    unlist(strsplit(input$queryInput, "\\s*,\\s*|\\s+"))
  })

  output$queryVect <- renderText({
    if(length(queries()) == 0){  }
    else
    {
      paste(queries(), collapse = ",\t")
    }
  })

  output$uploadedContent <- renderTable({
    head(data(),n =  5)
  })

  ####  Load Example Data
  observeEvent(input$loadExample,
               {
                 req(credentials()$user_auth)
                 updateTextInput(session, inputId = "queryInput", value =
                                   c("PspA", "Snf7", "PspB","PspC", "PspM", "PspN","DUF3046", "LiaI-LiaF-TM", "Toast-rack",
                                     "Tfu_1009","DUF1700-ahelical", "DUF1707-SHOCT"))
                 if(last_button_val() < input$loadExample)
                 {
                   data(example_data)

                   last_button_val(last_button_val()+1)
                 }

               })


  #### UI components for Upload tab ####
  output$uploadComponents <- renderUI({
    switch( input$inputType,
            "Full Data" = full_data_ui,
            "Fasta Sequence(s)" = fasta_input_ui,
            "Protein Accession Numbers" = accNum_input_ui,
            "Blast Results" = blast_input_ui,
            "Interproscan Results" = interpro_input_ui
    )
  })



  #### AccNum to Fasta ####
  # Generate a FASTA file from the input AccNum(s)

  # Convert the AccNums to a vector
  accnum_vect <- reactive({
    req(typeof(input$accNumTextInput) == "character")
    # Split input by commas or spaces
    unlist(strsplit(input$accNumTextInput, "\\s*,\\s*|\\s+"))
  })

  # Display the FASTA File
  output$generatedFasta <- renderText({
    fasta()
  })

  output$accnumMSA <- renderText({
    aligned_fasta()
  })


  dataTableFastaAccNums <- reactive({
    req(nrow(data()) > 0)
    reduced_col <- switch(input$fastaRepresentativeType,
                          "One per Species" = "Species",
                          "One per Lineage" = "Lineage"
                          )
    s <- RepresentativeAccNums(data(), reduced = reduced_col , accnum_col = "AccNum")
    print(s)
    s

  })

  output$DF2Fasta <- renderText({
    fasta()
  })
  output$DF2AlignedFasta <- renderText({
    aligned_fasta()
  })

  fasta <- reactiveVal("")
  aligned_fasta <- reactiveVal("")

  #### Observe Input Type ####
  observe({
    aligned_fasta_filepath("")
    fasta("")
    aligned_fasta("")
    # if(input$inputType == "Full Data")
    # {
    #   fasta("")
    #   aligned_fasta("")
    # }
    # else if(input$inputType == "Protein Accession Numbers")
    # {
    #   fasta("")
    #   aligned_fasta("")
    # }
  })


  #### FASTA input Upload ####
  observeEvent(
    input$fastaFileUpload,
    {
      req(credentials()$user_auth)
      f <- read_file(input$fastaFileUpload$datapath)
      fasta(f)

      # write(f, "tmp.fasta")
      # fasta_filepath("tmp.fasta")
      updateTextAreaInput(session, inputId = "fastaTextInput", label = "Enter Fasta Sequence", value = f)
    }
  )



  observeEvent(
    input$fasta2msaBtn,
    {
      write_file(input$fastaTextInput, "tmp.fasta")
      fasta(input$fastaTextInput)
      fasta_filepath("tmp.fasta")

      alignFasta(fasta_filepath(), tool = input$FastaAlignmentTool , outpath = "aligned_fasta.fasta")
      aligned_fasta(read_file("aligned_fasta.fasta"))
      aligned_fasta_filepath("aligned_fasta.fasta")
    }
  )

  output$fasta2msa <- renderText(
    aligned_fasta()
  )




  #### Fasta Buttons Full DF ####
  observeEvent(input$fullDF2Fasta,
               {
                 acc2fa(dataTableFastaAccNums(), out_path = "tmp.fasta")
                 fasta(read_file("tmp.fasta"))
               }
               )

  observeEvent(
    input$DF2msa,
      {
        alignFasta("tmp.fasta", tool = input$DFAlignmentTool, outpath = "aligned_fasta.fasta")
        aligned_fasta(read_file("aligned_fasta.fasta"))

        aligned_fasta_filepath("")
        aligned_fasta_filepath("aligned_fasta.fasta")
      }
  )

  #### FASTA buttons AccNum ####
  observeEvent(
    input$accnum2Fasta,
    {
    acc2fa(accnum_vect(), out_path = "tmp.fasta")
    fasta(read_file("tmp.fasta"))
    }
  )

  observeEvent(
    input$accnum2msa,
    {
      alignFasta("tmp.fasta", tool = input$accnumAlignmentTool, outpath = "aligned_fasta.fasta")
      aligned_fasta(read_file("aligned_fasta.fasta"))
      aligned_fasta_filepath("")
      aligned_fasta_filepath("aligned_fasta.fasta")
    }
  )



  ####  Load Example AccNums
  observeEvent(input$exampleAccNums,
               {
                 req(credentials()$user_auth)
                 lowerbound = sample(1:(nrow(example_data)-5), 1)
                 updateTextInput(session, inputId = "accNumTextInput", value =
                                   paste0(example_data$AccNum[ lowerbound: (lowerbound + 5)])  )
               })




  #### Redirect buttons on splash page ####
  # accnum button
  observeEvent(input$dAccNumBtn,
               {
               updateNavbarPage(session, "evolvrMenu" ,"upload")
               updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                                 choices = c(
                                   "Protein Accession Numbers",
                                   "Fasta Sequence(s)",
                                   "Blast Results",
                                   "Interproscan Results",
                                   "Full Data"
                                 ),
                                 selected = "Protein Accession Numbers"
                                 )
                 })

  # Splash page Fasta button
  observeEvent(input$dFastaBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                                   choices = c(
                                     "Protein Accession Numbers",
                                     "Fasta Sequence(s)",
                                     "Blast Results",
                                     "Interproscan Results",
                                     "Full Data"
                                   ),
                                   selected = "Fasta Sequence(s)"
                 )
               })




  ### Identify changes to queries and update dropdowns accordingly ###
  observe({
    print(input$inputType)
    if(input$inputType == "Full Data")
    {
      updateSelectInput(session, "mainSelect", label = "Protein", choices = append("All", queries()))
      updateSelectInput(session, "DASelect", label = "Protein", choices = append("All", queries()))
      updateSelectInput(session, "GCSelect", label = "Protein", choices = append("All", queries()))
    }
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
    prot_tc <- total_counts(DA_Prot(), cutoff = DACutoff(), column = "DomArch")
    prot_tc$Lineage = map(prot_tc$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
      unlist()
    prot_tc
  })

  ##   #   ###    #    ##
  # Total Counts for GC #
  ##   #   ###    #    ##
  GC_TotalCounts <- reactive({
    prot_tc <- total_counts(GC_Prot(), cutoff = GCCutoff(), column = "GenContext")
    prot_tc$Lineage = map(prot_tc$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
      unlist()
    prot_tc
  })



  ###### Heatmap ######
  output$DALinPlot <- renderPlot({
    lineage.DA.plot(DA_Prot(), colname = "DomArch", cutoff = DACutoff(), RowsCutoff = F,
                    color = input$DA_lin_color)
  })

  output$GCLinPlot <- renderPlot({
    lineage.DA.plot(GC_Prot(), colname = "GenContext", cutoff = GCCutoff(), RowsCutoff = F,
                    color = input$GC_lin_color)
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
    upset.plot(DA_Prot(), colname = "DomArch", cutoff = DACutoff())
  })

  output$GCUpsetP <- renderPlot({
    upset.plot(GC_Prot(), colname = "GenContext", cutoff = GCCutoff())
  })

  #### WordCloud ####
  output$DAwordcloud <- renderPlot({
    wordcloud_element(DA_Prot(), colname = "DomArch", cutoff = DACutoff(), UsingRowsCutoff = F)
  })

  output$GCwordcloud <- renderPlot({
    wordcloud_element(GC_Prot(), colname = "GenContext", cutoff = GCCutoff(), UsingRowsCutoff = F)
  })


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
  output$DANetwork <- renderVisNetwork({
    domain_network(prot = DA_Prot(), column = "DomArch.repeats",
                   domains_of_interest = network_domain_interest(),
                   cutoff = DACutoff(),
                   layout = input$networkLayout)
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

  #### Tree ####
  output$treePlot <- renderPlot({
    seq_tree(fasta_filepath = aligned_fasta_filepath())
  })


  #### Sunburst ####
  output$sunburst <- renderSunburst({
    lineage_sunburst(phylogeny_prot(), lineage_column = "Lineage", type = "sunburst", levels = input$levels)
  })


  #### Paralogs ####
  paralog_table <- reactive({
    find_paralogs(phylogeny_prot())
  })


  output$ParalogTable <- DT::renderDataTable({
    paralog_table()
  },extensions = c('FixedColumns'),
  selection = 'single',
  options = list(pageLength = 25,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = FALSE))

  ## Dialog box for rows selected in the paralog table
  observeEvent(input$ParalogTable_rows_selected,
               {
                 selected_paralogs = paralog_table()$AccNums[input$ParalogTable_rows_selected] %>% unlist()
                 selected_paralogs = paste(selected_paralogs, collapse = "|")
                 showModal(modalDialog(
                   title = "",
                   DT::renderDT({
                     phylogeny_prot() %>% filter(grepl(selected_paralogs,AccNum))
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
               }
  )










  #### Static Flowchart ####
  # output$flowchartImage <- renderImage( list(src = "data/figures/ISMB-Poster-FlowChart-evolvr.png"), deleteFile = F)


}


shinyApp(ui, server,options = options(shiny.maxRequestSize=100*1024^2))


# shiny::runExample("09_upload")