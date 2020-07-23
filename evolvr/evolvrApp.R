library(shiny)
library(shinyjs)
library(tidyverse)
library(wordcloud2)
library(sunburstR)
library(rmarkdown)
library(shinyBS)
library(shinyauthr)
library(pins)

setwd("..")
source("R/cleanup.R")
source("R/summarize.R")
source("R/plotting.R")
source("R/network.R")
source("R/pre-msa-tree.R")
source("scripts/tree.R")
source("evolvr/components.R")
source("evolvr/ui/resultSummaryComponents.R")
conflicted::conflict_prefer("strsplit", "base")
conflicted::conflict_prefer("count", "dplyr")

###########
###########
# Hide the token before making public
###########
###########
board_register_github(repo = "samuelzornchen/evolvR-Pins", token = "6b3b68bcb90353554fb7c84acfae3133abcfb7fa")


## Create Fasta Object so it can easily be pinned
FastaSeq <- setRefClass("FastaSeq", fields = list(sequence = "character"))

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
    source("evolvr/ui/resultSummaryUI.R")$value,
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

  data <- reactiveVal(data.frame())

  fasta_filepath <- reactiveVal("tmp.fasta")
  aligned_fasta_filepath <- reactiveVal("")


  DACutoff <- reactive({
    input$DA_Cutoff
  })

  GCCutoff <- reactive({
    input$GC_Cutoff
  })


  #### Observe what Data is available and hide/show tabs based on that ####
  # observe({
  #   if(nrow(data()) == 0)
  #   {
  #     hideTab(inputId = "QueryData", target = "mainData", session = session)
  #   }
  #   else
  #   {
  #     showTab(inputId = "QueryData", target = "mainData", session = session)
  #   }
  # })


#### Pin Fasta Data ####
  observeEvent(input$pinFasta,
               {

                 pin_name = input$pinName#[1]
                 print(pin_name)
                 if(nchar(pin_name) < 1 || nrow(pin_find(name = pin_name, board = "github")) != 0 )
                 {
                   createAlert(session, "invalidPin", content = "Error: Pin name already taken", append = FALSE)
                 }
                 else
                 {
                   fastaObject <- FastaSeq(sequence = fasta())
                   pin(x = fastaObject, name = pin_name, board = "github")
                 }
               }
               )

  #### UI components for Upload tab ####
  output$uploadComponents <- renderUI({
    switch( input$inputType,
            "Full Data" = full_data_ui,
            # "Fasta Sequence(s)" = fasta_input_ui,
            # "Protein Accession Numbers" = accNum_input_ui,
            "Blast Results" = blast_input_ui,
            "Interproscan Results" = interpro_input_ui,
            "AccNum/FASTA" = acc_fasta_ui
    )
  })

  #### UI Components for the Data tab ####
  output$dataTableComponent <- renderUI({
    if(nrow(data()) == 0)
    {
      noTableComponent
    }
    else
    {
      tableComponent
    }
  })

  output$fastaDataComponent <- renderUI({
    if(fasta() == "")
    {
      noFastaComponent
    }
    else
    {
      fastaComponent
    }
  })

  output$msaDataComponent <- renderUI({
    if(aligned_fasta() == "")
    {
      noMsaComponent
    }
    else
    {
      msaComponent
    }
  })

  observeEvent(
    input$dataTable2Upload,
    {
      updateNavbarPage(session, "evolvrMenu" ,"upload")
      updateSelectInput(session,
                        inputId = "inputType", label = "Input Type:",
                        choices = c(
                          "Blast Results",
                          "Interproscan Results",
                          "Full Data",
                          "AccNum/FASTA"
                        ),
                        selected = "Full Data")
    }

  )


  #### Full Data Upload ####
  observeEvent(input$fileUpload,
               {
                 req(credentials()$user_auth)
                 print("Upload")
                 df <- switch(input$fileType,
                              "tsv" = read_tsv(input$fileUpload$datapath),
                              "csv" = read_csv(input$fileUpload$datapath))
                 data(df)
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


  ## Observe data() to show button to redirect to result summary when available
  observe({
    shinyjs::hide("upload2RS")
    if((nrow(data()) != 0 & input$inputType == "Full Data") | (aligned_fasta() != "" & input$inputType == "AccNum/FASTA") )
    {
      shinyjs::show("upload2RS")
    }
  })

  observeEvent(input$upload2RS,
               {
                 updateNavbarPage(session, "evolvrMenu", selected = "resultSummary")
               }
  )


  ###### FASTA from Data Table ######
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


  #### AccNum Data Upload ####
  # Generate a FASTA file from the input AccNum(s)

  # Upload AccNum
  observeEvent(input$accNumUpload,
               {
                 req(credentials()$user_auth)
                 accNums <- read_file(input$accNumUpload$datapath)
                 updateTextAreaInput(session, "accNumTextInput", label = "Enter Protein Accession Number(s)",
                                     value = accNums)
               })


  # Convert the AccNums to a vector
  accnum_vect <- reactive({
    req(typeof(input$accNumTextInput) == "character")
    # Split input by commas or spaces
    unlist(strsplit(input$accNumTextInput, "\\s*,\\s*|\\s+|\\n"))
  })

  # Display the FASTA File
  # output$generatedFasta <- renderText({
  #   fasta()
  # })

  # output$accnumMSA <- renderText({
  #   aligned_fasta()
  # })



  fasta <- reactiveVal("")
  aligned_fasta <- reactiveVal("")

  #### Observe Input Type ####
  observe({
    aligned_fasta_filepath("")
    fasta("")
    aligned_fasta("")
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

      updateTextAreaInput(session, "msaText", label = "Paste Aligned FASTA Sequence", value = read_file("aligned_fasta.fasta"))
      updateCollapse(session, "accCollapse", open = "msa")
    }
  )

  output$fasta2msa <- renderText(
    aligned_fasta()
  )


  #### Msa Upload ####
  observeEvent(input$msaFileUpload,
               {
                 f = read_file(input$msaFileUpload$datapath)
                 write_file(f, "aligned_fasta.fasta")
                 aligned_fasta(f)
                 aligned_fasta_filepath("aligned_fasta.fasta")
                 updateTextInput(session, "msaText", label = "Paste  Aligned FASTA Sequence", value = f)
               }
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

      updateTextAreaInput(session, inputId = "fastaTextInput",
                          label = "Enter Fasta Sequence",
                          value = fasta())
      updateCollapse(session, "accCollapse", open = "fasta")
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
                 updateCollapse(session, "accCollapse", open = c("accnum"), close = c("msa", "fasta"))
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                                   choices = c(
                                     # "Protein Accession Numbers",
                                     # "Fasta Sequence(s)",
                                     "Blast Results",
                                     "Interproscan Results",
                                     "Full Data",
                                     "AccNum/FASTA"
                                   ),
                                   selected = "AccNum/FASTA"
                 )


               })

  # Splash page Fasta button
  observeEvent(input$dFastaBtn,
               {
                 updateCollapse(session, "accCollapse", open = c("fasta"), close = c("msa", "accnum"))
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                                   choices = c(
                                     # "Protein Accession Numbers",
                                     # "Fasta Sequence(s)",
                                     "Blast Results",
                                     "Interproscan Results",
                                     "Full Data",
                                     "AccNum/FASTA"
                                   ),
                                   selected = "AccNum/FASTA"
                 )

               })





  #### IPRScan Upload ####
  output$IPRScanData <- DT::renderDataTable(
    {
      ipr_cols <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
                    "SignAcc", "SignDesc", "StartLoc", "StopLoc", "Score",
                    "Status", "RunDate",
                    "IPRAcc", "IPRDesc", "GOAnn", "PathAnn")
      iprscanresults <- switch(
        input$fileTypeIPRScan,
        "tsv" = read_tsv(input$interproFileUpload$datapath , col_names = ipr_cols)
      )

      head(iprscanresults)
    }
  )




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


  ####



  #### Main Table ####
  output$mainTable = DT::renderDataTable({
    viewing_cols = c("AccNum", "DomArch", "GenContext","Lineage", "Species", "GeneName", "Length", "GCA_ID")
    if(input$mainSelect == "All")
    {
      data() %>% select(all_of(viewing_cols))
    }
    else
    {
      data() %>% filter(grepl(input$mainSelect, DomArch, ignore.case = T)) %>% select(all_of(viewing_cols))
    }
  })

  #### FASTA Data Output ####
  output$fastaDataText <- renderText({
    fasta()
  })

  output$msaDataText <- renderText({
    aligned_fasta()
  })


  #### Result Summary ####

  output$rs_DomArch_ui <- renderUI(
    {
      cols <- colnames(data())
      if("DomArch.repeats" %in% cols)
      {
        rs_DomArch_component
      }
    }
  )

  output$rs_GenContext_ui <- renderUI(
    {
      cols <- colnames(data())
      if("GenContext" %in% cols)
      {
        rs_GenContext_component
      }
    }
  )

  output$rs_Phylogeny_ui <- renderUI(
    {
      cols <- colnames(data())
      if("Lineage" %in% cols)
      {
        rs_sunburst_component
      }
      else if(fasta() != "")
      {
        rs_tree_component
      }
    }
  )


  # DA network
  output$rs_network <- renderVisNetwork({
    domain_network(data(), column = "DomArch.repeats",
                   domains_of_interest = ".*",
                   cutoff = 95
    )
  })

  # GC heatmap
  ######## Change to a row cutoff
  output$rs_gcHeatmap <- renderPlot({
    lineage.DA.plot(data(), colname = "GenContext", cutoff = 20)
    # RowsCutoff = T)
  })

  # sunburst
  output$rs_sunburst <- renderSunburst({
    lineage_sunburst(data(), "Lineage", levels = 2)
  })

  # tree
  output$rs_tree <- renderPlot({
    seq_tree(fasta_filepath = aligned_fasta_filepath())
  })

  #### action links Result Summary ####
  observeEvent(input$rs2DomArch,
               updateNavbarPage(session, "evolvrMenu" ,"domainArchitecture")
  )

  observeEvent(input$rs2GenContext,
               updateNavbarPage(session, "evolvrMenu", "genomicContext")
  )

  observeEvent(input$rs2Phylogeny,
               updateNavbarPage(session, "evolvrMenu", "phylogeny")
  )

  observeEvent(input$rs2PhylogenyTree,
               {
                 updateNavbarPage(session, "evolvrMenu", "phylogeny")
                 updateTabsetPanel(session, inputId = "phylo", selected = "Tree")
               }
  )


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