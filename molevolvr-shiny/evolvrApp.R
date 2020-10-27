library(shiny)
library(shinyjs)
library(tidyverse)
# library(wordcloud2)
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
source("R/Iprscan.R")
source("scripts/tree.R")
source("evolvr/pins/PinGeneration.R")
source("evolvr/components.R")
source("evolvr/ui/UIOutputComponents.R")
source("evolvr/ui/splashPageComponent.R")
conflicted::conflict_prefer("strsplit", "base")
conflicted::conflict_prefer("count", "dplyr")

conflicted::conflict_prefer("append","base")

###########
###########
# Hide the token before making public
###########
###########
board_register_github(repo = "JRaviLab/MolEvolvR-Pins", token = "YOUR_TOKEN_HERE")


## Create Fasta Object so it can easily be pinned
FastaSeq <- setRefClass("FastaSeq", fields = list(sequence = "character"))
## Create Accession Object so it can easily be pinned
AccNum <- setRefClass("AccNum", fields = list(accessions = "ANY") )

EX_HOMOLOGOUS_ACCNUMS = c("ANY95992.1","APP15780.1","	AGZ55339.1","AGZ57835.1", "AEB64964.1")
EX_FASTA_SEQS = "evolvr/ExampleData/fasta.fa"
EX_MSA = "evolvr/ExampleData/msa.aln"
EX_NON_HOMOLOGOUS_ACCNUMS = c("ANY95992.1")
EX_IPROUTPUT = "evolvr/ExampleData/iprscan5-20201022.txt"

ipr_cols <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
              "SignAcc", "SignDesc", "StartLoc", "StopLoc", "Score",
              "Status", "RunDate",
              "IPRAcc", "IPRDesc", "GOAnn", "PathAnn")


###
#### Users ####
###
user_base <- data.frame(
  user = c("molevolvr"),
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
                       , tags$div(class= "zoom",  "MolEvolvR")),
    id = 'evolvrMenu',
    inverse = T,
    #source("evolvr/ui/splashPageUI.R")$value,
    source("evolvr/ui/splashPageUI2.0.R")$value,
    source("evolvr/ui/uploadUI.R")$value,
    # source("evolvr/ui/analysisUI.R")$value,
    source("evolvr/ui/resultSummaryUI.R")$value,
    source("evolvr/ui/queryDataUI.R")$value,

    source("evolvr/ui/domainArchitectureUI.R")$value,
    # source("evolvr/ui/genomicContextUI.R")$value,
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


  # updateSelectInput(session, inputId = "iprDatabases", label = "DataBases", choices = c())

  ## Change title tab name (ie: chrome tab name) to EvolvR
  shinyjs::runjs('document.title = "EvolvR"')

  ## Redirect click on the app title to the upload-data page
  observeEvent(input$homeButton, updateNavbarPage(session, "evolvrMenu" ,"home"))

  last_button_val = reactiveVal(0)

  # Data from the full upload dataframe
  uFullData <- reactiveVal(data.frame())
  # Data table generated from the blast data
  uBlastData <- reactiveVal(data.frame())
  # Data table generated from the accnum/fasta data
  uAccFaData <- reactiveVal(data.frame())
  # Data table generated from the iprscan results data
  uIprData <- reactiveVal(data.frame())

  data <- reactive({
    switch(input$uploadTabs,
           "FullDataTab" = uFullData(),
           "BlastOutTab" = uBlastData(),
           "IprOutTab" = uIprData(),
           "AccFastaTab" =   uAccFaData()
    )
  })


  uFullFasta <- reactiveVal("")
  uBlastFasta <- reactiveVal("")
  uAccFaFasta <- reactiveVal("")
  uIprFasta <- reactiveVal("")

  uFullAln <- reactiveVal("")
  uBlastAln <- reactiveVal("")
  uAccFaAln <- reactiveVal("")
  uIprAln <- reactiveVal("")

  fasta <- reactive({
    switch(input$uploadTabs,
           "FullDataTab" = uFullFasta(),
           "BlastOutTab" = uBlastFasta(),
           "IprOutTab" = uIprFasta(),
           "AccFastaTab" = uAccFaFasta()
    )
  })

  aligned_fasta <- reactive({
    req(file.exists(aligned_fasta_filepath()))
    read_file(aligned_fasta_filepath())
    # switch(input$inputType,
    #        "Full Data" = uFullAln(),
    #        "Blast Results" = uBlastAln(),
    #        "Interproscan Results" = uAccFaAln(),ovl
    #        "AccNum/FASTA" = uIprAln()
    # )
  })

  aligned_fasta_filepath <- reactive({
    switch(input$uploadTabs,
           "FullDataTab" = "uFullAlnFile.fasta",
           "BlastOutTab" = "uBlastAlnFile.fasta",
           "IprOutTab" = "uIprAlnFile.fasta",
           "AccFastaTab" =   "uAccFaAlnFile.fasta"
    )


    # switch(input$inputType,
    #        "Full Data" = "uFullAlnFile.fasta",
    #        "Blast Results" = "uBlastAlnFile.fasta",
    #        "Interproscan Results" = "uAccFaAlnFile.fasta",
    #        "AccNum/FASTA" = "uIprAlnFile.fasta"
    # )
  })

  fasta_filepath <- reactive({
    switch(input$uploadTabs,
           "FullDataTab" = "uFullFastaFile.fasta",
           "BlastOutTab" = "uFullBlastFastaFile.fasta",
           "IprOutTab" = "uIprFastaFile.fasta",
           "AccFastaTab" =  "uAccFaFastaFile.fasta"
    )
  })




  DACutoff <- reactive({
    input$DA_Cutoff
  })

  GCCutoff <- reactive({
    input$GC_Cutoff
  })


  #### analysis ####
  lastAccAnalysis <- reactiveVal(0)

  analysis_AccNums <- reactive({
    req(typeof(input$accNumBLASTTextInput) == "character")
    # Split input by commas or spaces
    unlist(strsplit(input$accNumBLASTTextInput, "\\s*,\\s*|\\s+|\\n"))
  })

  #### Email ####
  # emailID <- function(from, to, ID)
  # {
  #   from <- paste0("<", from, ">")
  #   to <- paste0("<", from, ">")
  #   subject <- "!Important! MolEvolvr ID"
  #   body <- paste0("Your ID is: ", ID)
  #   sendmail(from, to, subject, body,
  #            control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
  # }

  #### Pin BLAST data ####
  observeEvent(input$pinAccBlast,
               {
                 if(length(analysis_AccNums()) == 0)
                 {
                   #### Makes some type of of alert here
                   print("Please input your Protein Accession Numbers")
                 }
                 # else if()
                 # Was the last submission less than 60 seconds ago?
                 else if( as.numeric(Sys.time()) - lastAccAnalysis() < 60 )
                 {
                   time_till_available = 60- ( as.numeric(Sys.time()) - lastAccAnalysis())
                   showNotification(
                     paste("You can only submit once per minute.\n",
                           time_till_available, "seconds till next available submission")
                     , type = "error")
                 }

                 else
                 {
                   pinName = GeneratePinName(length = 6, postfix = "_AccNum2BLAST", board = "github" )

                   # accNums <- AccNum(accessions = analysis_AccNums())
                   acc2fa(analysis_AccNums(), out_path = "fastaForBlast.fasta")

                   fastaBlast <- FastaSeq(sequence= read_file("fastaForBlast.fasta"))
                   print(fastaBlast$sequence)
                   pin(x = fastaBlast, name = pinName, board = "github")

                   lastAccAnalysis( as.numeric(Sys.time()))

                   createAlert(session, "successAccNumBlastAlert", NULL, title = substr(pinName, 0, 6),
                               content = "Use the above code to receive your results when they are ready"
                               , append = T, dismiss = T)

                 }
               }
  )


  fetched_paths <- reactiveVal(c())
  #### Fetch Pinned Output ####
  observeEvent(input$FetchAnalysisBtn,
               {
                 closeAlert(session, "FetchError")
                 pinName <- input$analysisCode
                 if(nchar(pinName) != 6 | nrow(pin_find( paste0(pinName,"_out"), "github")) == 0 )
                 {
                   # Show alert
                   createAlert(session, "FetchAlert", alertId = "FetchError", title = "Error: incorrect code",
                               content = "Please verify that the code you input is correct"
                               , append = F, dismiss = T, style = "danger")
                 }
                 else
                 {
                   # Set these outputs to some reactive vals
                   fetched_paths(pin_get( paste0(pinName,"_out") ,board = "github"))


                 }
               }
  )

  observe(
    {
      if(any(grepl("nr\\.1e-5\\.txt$", fetched_paths())))
      {
        showTab("FetchedDataTabs", target = "FetchedBlastTab")
      }
      else
      {
        hideTab("FetchedDataTabs", "FetchedBlastTab")
      }

      if(any(grepl("iprscan\\.tsv$", fetched_paths())))
      {
        showTab("FetchedDataTabs", target = "FetchedIprTab")
      }
      else
      {
        hideTab("FetchedDataTabs", "FetchedIprTab")
      }
    }
  )

  output$fetchedBlastOut <- DT::renderDataTable(
    {
      req(any(grepl("nr\\.1e-5\\.txt$", fetched_paths())))
      blastfile = fetched_paths()[grep("nr\\.1e-5\\.txt$", fetched_paths())]
      read_tsv(blastfile, col_names = F)
    }
  )

  output$fetchedIprOut <- DT::renderDataTable(
    {
      req(any(grepl("iprscan\\.tsv$", fetched_paths())))
      iprscanfile = fetched_paths()[grep("iprscan\\.tsv$", fetched_paths())]
      read_tsv(iprscanfile , col_names = F)
    }
  )


  output$splashUIComponent <- renderUI({
    req(credentials()$user_auth)
    splashUIComponent
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
      updateTabsetPanel(session, "uploadTabs", "FullDataTab")
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
                 uFullData(df)
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
    head(uFullData(),n =  5)
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
                   uFullData(example_data)

                   last_button_val(last_button_val()+1)
                 }

               })


  ## Observe data() to show button to redirect to result summary when available
  # observe({
  #
  #   if((nrow(data()) != 0 && input$inputType == "Full Data") || (aligned_fasta() != "" && input$inputType == "AccNum/FASTA") ||
  #      (nrow(ipr_data()) != 0 && input$inputType == "Interproscan Results"))
  #   {
  #     shinyjs::show("upload2RS")
  #   }
  #   else
  #   {
  #     shinyjs::hide("upload2RS")
  #   }
  # })

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






  #### FASTA input Upload ####
  observeEvent(
    input$fastaFileUpload,
    {
      req(credentials()$user_auth)
      f <- read_file(input$fastaFileUpload$datapath)
      uAccFaFasta(f)

      updateTextAreaInput(session, inputId = "fastaTextInput", label = "Enter Fasta Sequence", value = f)
    }
  )

  #### FASTA Load Example ####
  observeEvent(
    input$exampleFASTA,
    {
      f <- read_file(EX_FASTA_SEQS)
      uAccFaFasta(f)
      updateTextAreaInput(session, inputId = "fastaTextInput", label = "Enter Fasta Sequence", value = f)
    }
  )


  #### MSA Load Example ####
  ############## NEED A WAY OF CONFIRMING PASTES IN HERE
  observeEvent(
    input$exampleMSA,
    {
      f <- read_file(EX_MSA)
      uAccFaAln(f)
      write_file(f, aligned_fasta_filepath())

      updateTextAreaInput(session, "msaText", label = "Paste Aligned FASTA Sequence", value = f)
    }
  )


  observeEvent(
    input$fasta2msaBtn,
    {
      write_file(input$fastaTextInput, fasta_filepath())
      uAccFaFasta(input$fastaTextInput)

      alignFasta(fasta_filepath(), tool = input$FastaAlignmentTool , outpath = aligned_fasta_filepath())
      uAccFaAln(read_file(aligned_fasta_filepath()))
      # uAccFaAlnFile("uploadAccFaAln.fasta")


      print(uAccFaAln)

      updateTextAreaInput(session, "msaText", label = "Paste Aligned FASTA Sequence", value = read_file(aligned_fasta_filepath()))
      updateCollapse(session, "accCollapse", open = "msa")
    }
  )

  output$fasta2msa <- renderText(
    aligned_fasta()
  )





  # Full Fasta
  #### Fasta Buttons Full DF ####
  observeEvent(input$fullDF2Fasta,
               {
                 acc2fa(dataTableFastaAccNums(), out_path = fasta_filepath())
                 uFullFasta(read_file(fasta_filepath()))
               }
  )

  observeEvent(
    input$DF2msa,
    {
      alignFasta(fasta_filepath(), tool = input$DFAlignmentTool, outpath = aligned_fasta_filepath())
      uFullDataAln(read_file(aligned_fasta_filepath()))
    }
  )



  #### FASTA Accnum Upload ####
  #### FASTA buttons AccNum ####
  # Accnum to FASTA (AccNumFasta)
  observeEvent(
    input$accnum2Fasta,
    {
      acc2fa(accnum_vect(), out_path = fasta_filepath())
      uAccFaFasta(read_file(fasta_filepath()))

      print(fasta())

      updateTextAreaInput(session, inputId = "fastaTextInput",
                          label = "Enter Fasta Sequence",
                          value = fasta())
      updateCollapse(session, "accCollapse", open = "fasta")
    }
  )

  # Fasta to multiple sequence alignment (AccNumFasta)
  observeEvent(
    input$accnum2msa,
    {
      alignFasta(fasta_filepath() , tool = input$accnumAlignmentTool, outpath = aligned_fasta_filepath())
      uAccFaAln(read_file(aligned_fasta_filepath()))
    }
  )



  ####  Load Example AccNums
  observeEvent(input$exampleAccNums,
               {
                 req(credentials()$user_auth)
                 lowerbound = sample(1:(nrow(example_data)-5), 1)
                 updateTextInput(session, inputId = "accNumTextInput", value =  paste0(EX_HOMOLOGOUS_ACCNUMS, collapse =", "))
               })

  #### Msa Upload
  observeEvent(input$msaFileUpload,
               {
                 f = read_file(input$msaFileUpload$datapath)
                 write_file(f, aligned_fasta_filepath())
                 uAccFaAln(f)
                 updateTextInput(session, "msaText", label = "Paste  Aligned FASTA Sequence", value = f)
               }
  )



  #### Redirect buttons on splash page ####
  # accnum button
  observeEvent(input$dAccNumBtn,
               {
                 updateCollapse(session, "accCollapse", open = c("accnum"), close = c("msa", "fasta"))
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateTabsetPanel(session, "uploadTabs", selecte = "AccFastaTab")
                 # updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                 #                   choices = c(
                 #                     # "Protein Accession Numbers",
                 #                     # "Fasta Sequence(s)",
                 #                     "Blast Results",
                 #                     "Interproscan Results",
                 #                     "Full Data",
                 #                     "AccNum/FASTA"
                 #                   ),
                 #                   selected = "AccNum/FASTA"
                 # )


               })

  # Splash page Fasta button
  observeEvent(input$dFastaBtn,
               {
                 updateCollapse(session, "accCollapse", open = c("fasta"), close = c("msa", "accnum"))
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateTabsetPanel(session, "uploadTabs", "AccFastaTab")
                 # updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                 #                   choices = c(
                 #                     # "Protein Accession Numbers",
                 #                     # "Fasta Sequence(s)",
                 #                     "Blast Results",
                 #                     "Interproscan Results",
                 #                     "Full Data",
                 #                     "AccNum/FASTA"
                 #                   ),
                 #                   selected = "AccNum/FASTA"
                 # )

               })
  # Splashpage Ipr button
  observeEvent(input$dIprScanBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateTabsetPanel(session, "uploadTabs", "IprOutTab")
                 # updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                 #                   choices = c(
                 #                     # "Protein Accession Numbers",
                 #                     # "Fasta Sequence(s)",
                 #                     "Blast Results",
                 #                     "Interproscan Results",
                 #                     "Full Data",
                 #                     "AccNum/FASTA"
                 #                   ),
                 #                   selected = "Interproscan Results"
                 # )
               })
  # Splashpage BLAST button
  observeEvent(input$dBlastBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateTabsetPanel(session, "uploadTabs", "BlastOutTab")

                 # updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                 #                   choices = c(
                 #                     # "Protein Accession Numbers",
                 #                     # "Fasta Sequence(s)",
                 #                     "Blast Results",
                 #                     "Interproscan Results",
                 #                     "Full Data",
                 #                     "AccNum/FASTA"
                 #                   ),
                 #                   selected = "Blast Results"
                 # )
               })
  # Splashpage Full Data button
  observeEvent(input$dFullBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"upload")
                 updateTabsetPanel(session, "uploadTab", "FullDataTab")
                 # updateSelectInput(session, inputId = "inputType", label = "Input Type:",
                 #                   choices = c(
                 #                     # "Protein Accession Numbers",
                 #                     # "Fasta Sequence(s)",
                 #                     "Blast Results",
                 #                     "Interproscan Results",
                 #                     "Full Data",
                 #                     "AccNum/FASTA"
                 #                   ),
                 #                   selected = "Full Data"
                 # )
               })

  # Splashpage Domain Architecture button
  observeEvent(input$dDomArchBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"domainArchitecture")
               })
  # Splashpage Genomic Context button
  observeEvent(input$dGenContextBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"genomicContext")
               })
  # Splashpage Phylogeny button
  observeEvent(input$dPhyloBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"phylogeny")
               })
  # Splashpage Homology button
  observeEvent(input$dHomologBtn,
               {
                 updateNavbarPage(session, "evolvrMenu" ,"datatable")
               })

  #### IPRScan Upload ####


  ipr_data = reactiveVal(data.frame())
  ipr_filepath = reactiveVal("")

  observe(
    {
      if(input$evolvrMenu == "resultSummary")
      {
        # req("Analysis" %in% colnames(ipr_data()))
        options <- ipr_data()$Analysis %>% unique()
        # updateSelectizeInput(session, inputId = "iprDatabases", label = "DataBases", choices = options, selected = options)
        updateSelectInput(session, inputId = "iprDatabases", label = "DataBases", choices = options, selected = options)
      }
      else if(input$evolvrMenu == "domainArchitecture")
      {
        options <- ipr_data()$Analysis %>% unique()
        # updateSelectizeInput(session, inputId = "iprDatabases", label = "DataBases", choices = options, selected = options)
        updateSelectInput(session, inputId = "da_iprDatabases", label = "DataBases", choices = options, selected = options)
      }
    }
  )

  observeEvent(input$interproFileUpload, {
    ipr_data(switch(
      input$fileTypeIPRScan,
      "tsv" = read_tsv(input$interproFileUpload$datapath , col_names = ipr_cols)
    ))
    write_tsv(ipr_data(), "ipr_out_file.txt")
    ipr_filepath("ipr_out_file.txt")
  })

  #### Load IPR Example Data ####
  observeEvent(input$loadIPRExample,
               {
                 ipr_data(read_tsv(EX_IPROUTPUT, col_names = ipr_cols))
                 ipr_filepath(EX_IPROUTPUT)
               }
               )


  output$IPRScanData <- DT::renderDataTable(
    {
      req(nrow(ipr_data()) != 0)
      ipr_data()
    }
  )




  ### Identify changes to queries and update dropdowns accordingly ###
  observe({
    if(input$uploadTabs == "FullDataTab")
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
      else if(aligned_fasta() != "")
      {
        rs_tree_component
      }
    }
  )

  output$rs_IprScan_ui <- renderUI(
    {
      if(input$uploadTabs ==  "IprOutTab")
      {
        rs_iprscan_component
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


  # IprVis
  output$rs_IprGenes <- renderPlot({
    req(nrow(ipr_data()) != 0)
    ipr2domarch(infile_ipr = ipr_filepath(),
                PfamClans_path = "evolvr/TestData/Pfam-A.clans.txt",
                analysis = input$iprDatabases, group_by = input$iprVisType,
                topn = 20 ##### What does this parameter do?
    )
  })

  output$da_IprGenes <- renderPlot({
    req(nrow(ipr_data()) != 0)
    ipr2domarch(ipr_filepath(),
                PfamClans_path = "evolvr/TestData/Pfam-A.clans.txt",
                analysis = input$da_iprDatabases, group_by = input$da_iprVisType,
                topn = 20 ##### What does this parameter do?
                )
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


  #### Handle button dependencies ####


  # BLAST CLUST
  observeEvent(input$Acc_blastclustCB,
               {
                 # Selecting blast clust requires DeltaBlast
                 if(input$Acc_blastclustCB)
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_deltaCB", label = "DeltaBLAST", value = T)
                   updatePrettyCheckbox(session, inputId = "Acc_acc2faCB", label = "Acc2FASTA", value = T)
                 }
                 # Unselecting it means that rps blast and ipr scan can't be run
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_rpsblastCB", label = "RPSBLAST", value = F)
                   updatePrettyCheckbox(session, inputId = "Acc_iprscanCB", label = "InterproScan", value = F)
                 }
               })

  # DeltaBLAST
  observeEvent(input$Acc_deltaCB,
               {
                 if(!input$Acc_deltaCB)
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = F)
                   updatePrettyCheckbox(session, inputId = "Acc_acc2faCB", label = "Acc2FASTA", value = F)
                 }

               })


  observeEvent(input$Acc_acc2faCB,
               {
                 if(input$Acc_acc2faCB)
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_deltaCB", label = "DeltaBLAST", value = T)
                 }
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = F)
                 }
               }
  )
  # Iprscan
  observeEvent(input$Acc_iprscanCB,
               {
                 if(input$Acc_iprscanCB){
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = T)
                 }
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_ipr2daCB", label = "InterproScan2DA", value = F)
                 }
               })

  # rpsblast
  observeEvent(input$Acc_rpsblastCB,
               {
                 if(input$Acc_rpsblastCB){
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = T)
                 }
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_rps2daCB", label = "RPSBLAST2DA", value = F)
                 }
               })
  # rps2da
  observeEvent(input$Acc_rps2daCB,
               {
                 if(input$Acc_rps2daCB){
                   updatePrettyCheckbox(session, inputId = "Acc_rpsblastCB", label = "RPSBLAST", value = T)
                 }
               })

  # ipr2da
  observeEvent(input$Acc_ipr2daCB,
               {
                 if(input$Acc_ipr2daCB){
                   updatePrettyCheckbox(session, inputId = "Acc_iprscanCB", label = "InterproScan", value = T)
                 }
               })




}


shinyApp(ui, server,options = options(shiny.maxRequestSize=100*1024^2))