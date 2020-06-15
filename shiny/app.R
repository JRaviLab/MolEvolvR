library(shiny)
library(tidyverse)
library(DT)
library(rmarkdown)
library(shinydashboard)
library(shinyjqui)
library(wordcloud)
library(shinyjs)
library(shinyauthr)
library(V8)
conflicted::conflict_prefer("intersect", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("strsplit", "base")
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("box", "shinydashboard")
conflicted::conflict_prefer("upset", "UpSetR")
setwd("..")
source("shiny/PSP_Web_Data.R")
source("R/plotting.R")
source("R/network.R")
# source("R/GC_network_directed.R")
#source("R/cleanup.R")
#source("R/reverse_operons.R")
source("shiny/shinyfunctions.R")
source("shiny/initialCutoff.R")
source("shiny/legend_text.R")

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



####
#### UI ####
####
ui <- tagList(
  shinyjs::useShinyjs(),
  tags$head(includeScript("shiny/logout-button.js")),
  tags$head(includeScript("shiny/github-button.js")),
  #tags$head(includeScript("shiny/navbar-title-hyperlink.js")),
  tags$head(
    tags$style(HTML("
                      .innerbox {
                        /*border: 2px solid black;*/
                        box-shadow: 2px 2px 3px 3px #ccc;
                        margin: auto;
                        padding: 20px;
                      }

                      .bord {
                        margin: auto;
                        padding: 20px;
                      }

                      .lightblue-link{
                        color:#11aad9;
                      }

                      .iMargin{
                        margin: 5px;
                      }

                     a{
                        color: #141414;
                        text-decoration:none;
                      }


                      a:visited{
                        color:none;
                      }

                      .noDec{
                        color:white;
                        text-decoration:none;
                      }

                      .zoom:hover {
                        /*color:white;*/
                        transform: scale(1.3);
                      }

                      "))
  ),
  navbarPage(
    # must turn shinyjs on
    title =  actionLink(inputId = "homeButton",
                        style ="text-decoration:none; color:white;"
                        , tags$div(class= "zoom",  "Psp")),
    id = "pspMenu",
    inverse = TRUE,
    selected = "appInfo",
    # collapsible = T,
    #App title
    # header,
    #body
    source("shiny/ui/appInfoUI.R")$value,
    source("shiny/ui/queryDataUI.R")$value,
    source("shiny/ui/domainArchitectureUI.R")$value,
    source("shiny/ui/genomicContextUI.R")$value,
    source("shiny/ui/phylogenyUI.R")$value,
    source("shiny/ui/aboutUI.R")$value
  )
)

####
#### Server ####
server <- function(input, output,session){

  ## Change title tab name (ie: chrome tab name) to Psp-Evolution
  shinyjs::runjs('document.title = "Psp-Evolution;"')

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


  observeEvent(input$homeButton, updateNavbarPage(session, "pspMenu" ,"appInfo"))

  ### Reorder all to be ordered by query
  reorder_all <- function(prot, queries)
  {
    all_reordered <- c()
    for(q in queries)
    {
      all_reordered <- append(all_reordered, grep(q, all$DomArch))
    }

    all_reordered <- unique(all_reordered)

    all_re <- all[all_reordered,]
    return(all_re)
  }

  #Reactive expression used to determine which data table should be displayed
  #based on selected input
  pspTable<- reactive({
    req(credentials()$user_auth)
    switch(input$proSelec,
           "All" = all %>% select(viewing_cols) %>% distinct() %>% reorder_all(queries),
           "DUF1700-ahelical" = DUF1700%>% select(viewing_cols) %>% distinct(),
           "DUF1707-SHOCT" = DUF1707%>% select(viewing_cols) %>% distinct(),
           "PspA-Snf7" = pspa_snf7%>% select(viewing_cols) %>% distinct(),
           "PspA" = pspa%>% select(viewing_cols) %>% distinct(),
           "Snf7" = snf7%>% select(viewing_cols) %>% distinct(),
           #"Psp-AA" = psp_aa%>% select(viewing_cols) %>% distinct(),
           "PspB" = pspb%>% select(viewing_cols) %>% distinct(),
           "PspC" = pspc%>% select(viewing_cols) %>% distinct(),
           "PspM" = pspm%>% select(viewing_cols) %>% distinct(),
           "PspN" = pspn%>% select(viewing_cols) %>% distinct(),
           "DUF3046" = DUF3046 %>% select(viewing_cols) %>% distinct(),
           "LiaI-LiaF-TM" = liai_liaf%>% select(viewing_cols) %>% distinct(),
           "Toast-rack" = toast_rack%>% select(viewing_cols) %>% distinct(),
           "Tfu-1009" = tfu_1009%>% select(viewing_cols) %>% distinct()
    )
  })
  #Render the Data table for selected protein
  output$proTable <- DT::renderDT({
    req(credentials()$user_auth)
    d = DT::datatable(
      paged_table(pspTable() ), extensions = c('FixedColumns',"FixedHeader"),
      options = list(pageLength = 100,
                     #The below line seems to disable other pages and the search bar
                     #dom = 't',
                     scrollX = TRUE,
                     paging=TRUE,
                     fixedHeader= TRUE,
                     fixedColumns = list(leftColumns = 2, rightColumns = 0)))
  })

  #### Query Heatmap ####
  output$queryHeatmap <- renderPlot({
    req(credentials()$user_auth)
    lineage.Query.plot(query_data = all, queries = queries, colname = "DomArch", cutoff = 100)
  })



  #### Reactive expression determining data for which plotting is based. Controlled by dropdown
  plotting_prot <-  reactive({
    if(input$pspMenu == "domainArchitecture"){
      switch(input$DAlinSelec,
             "All" = all,
             "DUF1700-ahelical" = DUF1700,
             "DUF1707-SHOCT" = DUF1707,
             "PspA" = pspa,
             "Snf7" = snf7,
             "PspB" = pspb,
             "PspC" = pspc,
             "PspM" = pspm,
             "PspN" = pspn,
             "LiaI-LiaF-TM" = liai_liaf,
             "Toast-rack" = toast_rack,
             "Tfu-1009" = tfu_1009)
    }
    # else if(input$mainTabs == "genomicContext")
    # {
    else{
      switch(input$GClinSelec,
             "All" = all,
             "DUF1700-ahelical" = DUF1700,
             "DUF1707-SHOCT" = DUF1707,
             # "PspA-Snf7" = pspa,
             "PspA" = pspa,
             "Snf7" = snf7,
             "PspB" = pspb,
             "PspC" = pspc,
             "PspM" = pspm,
             "PspN" = pspn,
             "LiaI-LiaF-TM" = liai_liaf,
             "Toast-rack" = toast_rack,
             "Tfu-1009" = tfu_1009)
    }
  })

  prot_word_percents <- reactive({
    if(input$pspMenu == "domainArchitecture"){
      #### Should these be on DomArch.repeats for Network?
      max_word_percents(plotting_prot(), "DomArch")
    }
    # else if(input$mainTabs == "genomicContext"){
    #   max_word_percents(plotting_prot(), "GenContext")
    # }
  })


  DA_cutoff_val <- reactive({
    if(cutoff_status() == "Percent")
    {
      ## Percent cutoff, no need for conversion
      input$DA_Cutoff
    }
    else
    {
      # Cutoff type is Row, convert it to the percentage
      if(input$DALin_data == "Heatmap")
      {

        rownumber_to_cutoff(plotting_prot(),input$DA_Cutoff, col = "DomArch")
      }

      else
      {
        #### Do some other row cutoff based on the words
        #input$DA_Cutoff
        100-max_word_percents(plotting_prot(), "DomArch")[input$DA_Cutoff, 'MaxPercent']
      }

    }
  })

  GC_cutoff_val <- reactive({
    input$GC_Cutoff
  })

  ### Adaptive cutoff value to be used for plotting functions
  ### Converts the cutoff slider's value into respective percentage value
  ### if the cutoff status is "Row"
  cutoff_val <- reactive({
    if(cutoff_status() == "Percent")
    {
      # Cutoff type is Percentage already, so use original slider value
      if(input$mainTabs == "domainArchitecture"){
        input$DA_cutoff
      }
      else if(input$mainTabs == "genomicContext")
      {
        input$GC_Cutoff
      }
    }
    else
    {
      # Cutoff type is Row, convert it to the percentage
      if(input$lin_data == "Heatmap")
      {
        if(input$linSelec == "All")
        {
          if(input$DA_GC== "Domain Architecture"){
            100 - query_DA_row_CutoffPercs$maxPercent[input$cutoff]
          }
          else
          {
            100 - query_GC_row_CutoffPercs$maxPercent[input$cutoff]
          }
        }
        else{
          if(input$DA_GC== "Domain Architecture"){
            column = "DomArch"
            rmAstrk = F
          }
          else
          {
            column = "GenContext"
            rmAstrk = T
          }
          rownumber_to_cutoff(plotting_prot(),input$cutoff, col = column)
        }
      }

      else
      {
        # Cutoff to percent for upset and wordclouds and network
        100 - prot_word_percents()$MaxPercent[input$cutoff]
      }
    }
  })

  # #### Total Number of Rows the currently selected protein has
  plotting_prot_maxRows <- reactive(
    {
      if(input$DALin_data == "Heatmap")
      {
        # Current tab is Heatmap

        RowNums(plotting_prot(), column = "DomArch")

      }
      else
      {
        # Current tab is not Heatmap
        nrow(max_word_percents(plotting_prot(), "DomArch"))
        # 100
      }
    }
  )

  # #### Initial number of rows that should be selected For the selected protein
  init_rowCutoff <- reactive({
    min(10,plotting_prot_maxRows())
  })

  # ### Current type of cutoff that is being used either
  cutoff_status <- reactiveVal("Percent")

  rows_cutoff <- reactiveVal(FALSE)

  #Observer used to determine initial heatmap slider whenever protein changes for DA
  observe({
    if(cutoff_status() == "Percent"){
      label_val = "Percent Cutoff"
      max_val = 100

      #### Below is broken?
      cutoff_init <- top_n_rows_cutoff(plotting_prot(), column = "DomArch", n_rows = 10)
    }
    else
    {
      # Use Row Cutoffs: update
      label_val = "Row Cutoff"
      max_val = plotting_prot_maxRows()
      cutoff_init = init_rowCutoff()
    }
    updateSliderInput(session, label = label_val, inputId = "DA_Cutoff",min=1, max= max_val, value=cutoff_init)
  })

  #Observer used to determine initial linplot slider whenever protein changes for GC
  observe({
    init_val = top_n_rows_cutoff(plotting_prot(), "GenContext", 10)
    updateSliderInput(session, "Percent Cutoff", inputId = "GC_Cutoff", min = 1, max = 100, value = init_val)
  })



  # ### Observe when the CutoffSwitch is pressed, and toggle the text and set to the correct status
  observeEvent(input$DACutoffSwitch,
               {
                 if(cutoff_status() == "Percent" )
                 {
                   # The current cutoff status is Percent, switch to Rows
                   cutoff_status("Row")
                   rows_cutoff(TRUE)

                   # Update Action button text to "Row Cutoff"
                   updateActionButton(session, "DACutoffSwitch",
                                      label = ("Percent Cutoff"))
                 }
                 else if(cutoff_status() == "Row")
                 {
                   cutoff_status("Percent")
                   rows_cutoff(FALSE)

                   # Update Action button text to "Percent Cutoff"
                   updateActionButton(session, "DACutoffSwitch",
                                      label = ("Row Cutoff"))
                 }
               })

  rows_cutoff <- reactiveVal(FALSE)



  ####
  ##### Render the heatmap #####
  ####

  output$DALinPlot <- renderPlot({
    req(credentials()$user_auth)
    lineage.DA.plot(plotting_prot(), colname = "DomArch", cutoff = DA_cutoff_val(), RowsCutoff = rows_cutoff())

  })

  output$GCLinPlot <- renderPlot({
    req(credentials()$user_auth)
    lineage.DA.plot(plotting_prot(), colname = "GenContext", cutoff = GC_cutoff_val(), RowsCutoff = rows_cutoff())

  })



  ##   #   ###    #    ##
  # Total Counts for DA #
  ##   #   ###    #    ##
  DA_TotalCounts <- reactive({
    prot_tc <- total_counts(plotting_prot(), cutoff = DA_cutoff_val(), column = "DomArch")
    prot_tc$Lineage = map(prot_tc$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
      unlist()
    prot_tc
  })



  ####
  ##### Render the datatable for the lineage counts #####
  ####
  DAlin_count_table <- reactive({
    DA_TotalCounts() %>% group_by(DomArch, totalcount, CumulativePercent) %>%
      summarize(LineageCount = n()) %>%
      select(DomArch, LineageCount, totalcount, CumulativePercent) %>%
      arrange(-totalcount)
  })
  output$DALinTable <- DT::renderDT({
    req(credentials()$user_auth)
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
        # regexDAs = paste0( DAlin_count_table()$DomArch[input$DALinTable_rows_selected],collapse = "$|^") %>%
        #   str_replace_all(pattern = "\\+", replacement = "\\\\+") %>%
        #   str_replace_all(pattern = "\\(", replacement = "\\\\(") %>%
        #   str_replace_all(pattern = "\\)", replacement = "\\\\)") %>%
        #   unlist()
        # regexDAs = paste0("^", regexDAs, "$")

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




  ##   #   ###    #    ##
  # Total Counts for GC #
  ##   #   ###    #    ##
  GC_TotalCounts <- reactive({
    prot_tc <- total_counts(plotting_prot(), cutoff = GC_cutoff_val(), column = "GenContext")
    prot_tc$Lineage = map(prot_tc$Lineage, function(x) str_replace_all(string = x,pattern = ">", replacement = "_")) %>%
      unlist()
    prot_tc
  })



  ####
  ##### Render the datatable for the lineage counts #####
  ####
  GClin_count_table <- reactive({
    GC_TotalCounts() %>% group_by(GenContext, totalcount, CumulativePercent) %>%
      summarize(LineageCount = n()) %>%
      select(GenContext, LineageCount, totalcount, CumulativePercent) %>%
      arrange(-totalcount)
  })
  output$GCLinTable <- DT::renderDT({
    req(credentials()$user_auth)
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

  #### Upset Plots ####
  output$DAUpsetP <-renderPlot({
    req(credentials()$user_auth)
    upset.plot(plotting_prot(), cutoff = DA_cutoff_val(), colname = "DomArch", RowsCutoff = rows_cutoff())
  })

  output$GCUpsetP <- renderPlot({
    req(credentials()$user_auth)
    upset.plot(plotting_prot(), cutoff = GC_cutoff_val(), colname = "GenContext", RowsCutoff = rows_cutoff())
  })


  #### Reactive expression determining domain of interest for plotting domain networks
  network_domain_interest <-  reactive({
    switch(input$DAlinSelec,
           "All" = ".*",  #c("DUF1700-ahelical","DUF1707-SHOCT","PspA", "Snf7","PspB","PspC", "PspM","PspN","LiaI-LiaF-TM","Toast-rack"),
           "DUF1700-ahelical" = "DUF1700-ahelical",
           "DUF1707-ahelical" = "DUF1707-SHOCT",
           "PspA-Snf7" = c("PspA","Snf7"),
           "PspA" = "PspA",
           "Snf7" = "Snf7",
           "PspB" = "PspB",
           "PspC" = "PspC",
           "PspM" = "PspM",
           "PspN" = "PspN",
           "LiaI-LiaF-TM" = "LiaI-LiaF-TM",
           "Toast-rack" = "Toast-rack",
           "Tfu-1009" = "Tfu_1009")
  })


  #### Network Output ####
  output$DANetwork <- renderPlot({
    req(credentials()$user_auth)
    domain_network(prot = plotting_prot(), column = "DomArch.repeats",
                   domains_of_interest = network_domain_interest(),
                   cutoff = DA_cutoff_val(),
                   layout = "auto")
  })


  #### Wordcloud ####
  output$DAwordcloud <- renderWordcloud2({
    req(credentials()$user_auth)
    wordcloud_element(query_data = plotting_prot(), colname = "DomArch",
                      cutoff = DA_cutoff_val(), UsingRowsCutoff = rows_cutoff())
  })

  output$GCwordcloud <-renderWordcloud2({
    req(credentials()$user_auth)
    wordcloud_element(query_data = plotting_prot(), colname = "GenContext",
                      cutoff = GC_cutoff_val(), UsingRowsCutoff = rows_cutoff())
  })

  # output$DAwordcloud <- renderPlot({
  #   wordcloud_element(query_data = plotting_prot(), colname = "DomArch", cutoff = DA_cutoff_val())
  # })
  #
  # output$GCwordcloud <-renderPlot({
  #   wordcloud_element(query_data = plotting_prot(), colname = "GenContext", cutoff = GC_cutoff_val())
  # })

  # output$wordcloud <- renderWordcloud2({
  #   req(credentials()$user_auth)
  #   if(input$DA_GC == "Genomic Context"){
  #     wordcloud_element(query_data = plotting_prot(), colname = "GenContext", cutoff = cutoff_val())
  #   }
  #   else{
  #     wordcloud_element(query_data = plotting_prot(), colname = "DomArch", cutoff = cutoff_val())
  #   }
  # })

  #### Legends ####

  ### Observe changes in tabset panel to determine legend text ###
  DA_legend_txt <- reactive({
    switch( input$DALin_data,
            "Heatmap" = ifelse(input$DAlinSelec == "All", AllLineageHeatmap_txt, LineageHeatmap_txt),
            "Upset" = Upset_txt,
            "Lintable" = LinTable_txt,
            "Network_WC" = paste(Network_txt, Wordcloud_txt, sep = "\n")
    )
  })

  GC_legend_txt <- reactive({
    switch( input$GCLin_data,
            "Heatmap" = ifelse(input$GClinSelec == "All", AllLineageHeatmap_txt, LineageHeatmap_txt),
            "Upset" = Upset_txt,
            "Lintable" = LinTable_txt,
            "Network_WC" = paste(Network_txt, Wordcloud_txt, sep = "\n")
    )
  })

  output$DALegend <- renderText({
    req(credentials()$user_auth)
    DA_legend_txt()
  })

  output$GCLegend <- renderText({
    req(credentials()$user_auth)
    GC_legend_txt()
  })


  phylogeny_prot <- reactive({
    switch(input$alignSelec,
           "All" = all,
           "DUF1700" = DUF1700,
           "DUF1707" = DUF1707,
           "PspA-Snf7" = pspa,
           "Psp-AA" = psp_aa,
           "PspB" = pspb,
           "PspC" = pspc,
           "PspM" = pspm,
           "PspN" = pspn,
           "LiaI-LiaF-TM" = liai_liaf,
           "Toast-rack" = toast_rack,
           "Tfu-1009" = tfu_1009
    )
  })

  observe({
    ### PspA, PspA-Snf7, PspB, PspB, PspC, Snf7, Toast-rack, LiaI-LiaF
    if(input$phylo == "Tree"){
      choices = c(
        "LiaI-LiaF-TM.allfa50",
        "LiaI-LiaF-TM_LiaFN.2",
        "LiaI-LiaF-TM_LiaI.1",
        "LiaI-LiaF-TM_PspC.3",
        "PspA Only",
        "PspA Snf7 Gismo",
        "PspA Snf7",
        "PspB Gismo",
        "PspC Gismo",
        "Snf7 Only",
        "Toast-rack DUF2154-LiaF",
        "Toast-rack DUF2807",
        "Toast-rack DUF4097",
        "Toast-rack PspC-Cterm"
      )
      selected = "PspA Snf7"
    }
    else{
      choices = c("PspA-Snf7",
                  "PspB",
                  "PspC",
                  "PspN",
                  "LiaI-LiaF-TM",
                  "Toast-rack"
      )
      selected = "PspA-Snf7"
    }
    updateSelectInput(session, inputId = "alignSelec",
                      choices = choices, selected = selected)
  })

  #### Tree ####
  output$msaTree <- renderUI({
    print(input$phylo)
    req(credentials()$user_auth)
    switch(input$alignSelec,
           "LiaI-LiaF-TM.allfa50" = tags$iframe(style="height:600px; width:100%", src="FigTrees/LiaI-LiaF-TM.allfa50_tree.pdf", seamless=T),
           "LiaI-LiaF-TM_LiaFN.2"= tags$iframe(style="height:600px; width:100%", src="FigTrees/LiaI-LiaF-TM_LiaFN.2_tree.pdf", seamless=T),
           "LiaI-LiaF-TM_LiaI.1" = tags$iframe(style="height:600px; width:100%", src="FigTrees/LiaI-LiaF-TM_LiaI.1_tree.pdf", seamless=T),
           "LiaI-LiaF-TM_PspC.3" = tags$iframe(style="height:600px; width:100%", src="FigTrees/LiaI-LiaF-TM_PspC.3_tree.pdf", seamless=T),
           "PspA Only" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspa_only.1_tree.pdf", seamless=T),
           "PspA Snf7" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspa_snf7_tree.pdf", seamless=T),
           "PspA Snf7 Gismo" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspa_snf7.gismo_tree.pdf", seamless=T),
           "PspB Gismo" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspb.gismo_tree.pdf", seamless=T),
           "PspC Gismo" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspc.gismo_tree.pdf", seamless=T),
           "Snf7 Only" = tags$iframe(style="height:600px; width:100%", src="FigTrees/snf7_only.1_tree.pdf", seamless=T),
           "Toast-rack DUF2154-LiaF" = tags$iframe(style="height:600px; width:100%", src="FigTrees/Toast-rack_DUF2154-LiaF.1_tree.pdf", seamless=T),
           "Toast-rack DUF2807" = tags$iframe(style="height:600px; width:100%", src="FigTrees/Toast-rack_DUF2807.2_tree.pdf", seamless=T),
           "Toast-rack DUF4097" = tags$iframe(style="height:600px; width:100%", src="Toast-rack_DUF4097-LiaG.3.allfa.2nd_tree.pdf", seamless=T),
           "Toast-rack PspC-Cterm" = tags$iframe(style="height:600px; width:100%", src="FigTrees/Toast-Rack_PspC-Cterm.4_tree.pdf", seamless=T)
    )
  })


  #### Sunburst ####

  # output$sund2b <- renderSund2b({
  #   req(credentials()$user_auth)
  #   lineage_sunburst(pspa, lineage_column = "Lineage", type = "sund2b")
  # })

  output$sunburst <- renderSunburst({
    req(credentials()$user_auth)
    req(input$levels)
    lineage_sunburst(phylogeny_prot(), lineage_column = "Lineage", type = "sunburst", levels = input$levels)
  })


  output$msaPlot <- renderUI({
    req(credentials()$user_auth)
    tags$iframe(style="height:600px; width:100%", src="pspa_reduced.fasta.pdf", seamless=T)
  })





  #### Paralogs ####
  paralog_table <- reactive({
    find_paralogs(phylogeny_prot())
  })

  output$ParalogTable <- DT::renderDataTable({
    req(credentials()$user_auth)
    paralog_table()
  },extensions = c('FixedColumns'),
  selection = 'single',
  options = list(pageLength = 25,
                 #The below line seems to disable other pages and the search bar
                 #dom = 't',
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




  #Reactive expresion to change file name depending on which protein is selected
  fileNam <- reactive({
    if(input$downloadType == "tsv"){
      paste(input$proSelec, ".txt", sep = "")
    }
    else  paste(input$proSelec, ".csv", sep = "")
  })

  ##### Download main data table #####
  # Should it download cleaned data? or all rows?
  output$downloadData <- downloadHandler(
    filename <-  function() {
      fileNam()
    },
    content = function(file) {

      if(input$downloadType == "tsv"){
        write_tsv(pspTable(), file)
      }
      else  write.csv(pspTable(), file, row.names = FALSE)
    }
  )


  #### About Text ####
  output$aboutApp <- renderUI({
    req(credentials()$user_auth)
    tagList(
      h4("Background"),
      "This web app was built to provide a visual and interactive supplement for the Psp Evolution paper.",



      h4("Code"),
      tags$span(
        "Code and data used to generate this Shiny app are available on ",
        a(href = "https://github.com/JRaviLab/the-approach",
          class ="lightblue-link", "GitHub.")
      )
    )
  })

  output$aboutAbstract <- renderUI({
    req(credentials()$user_auth)

    tagList(
    h1("Phage-shock-protein (Psp) Envelope Stress Response:
                           Evolutionary History & Discovery of Novel Players"),

    h4("Janani Ravi", tags$sup("1,2*"),", Vivek Anantharaman", tags$sup("3")
       , ", Samuel Zorn Chen", tags$sup("1"), ", Pratik Datta", tags$sup("2"),
       ", L Aravind", tags$sup("3*"),", Maria Laura Gennaro", tags$sup("2*"),"."),
    tags$sup("1"),"Pathobiology and Diagnostic Investigation, Michigan State University, East Lansing, MI;",
    tags$sup("2"), "Public Health Research Institute, Rutgers University, Newark, NJ; ",
    tags$sup("3"),"National Center for Biotechnology Information, National Institutes of Health, Bethesda, MD.",
    tags$br(),
    "*Corresponding authors. janani@msu.edu; aravind@nih.gov; marila.gennaro@rutgers.edu ",

    h2('Abstract:'),

    p(abstract, style = "font-size:120%"))

  })



  # # #                             # # #
  ##### Download Data for Lin Table #####
  # # #                             # # #

  ## DA lin table download
  output$DAdownloadCounts <- downloadHandler(
    req(credentials()$user_auth),
    filename = function(){
      extension <- "-da_lin_counts"
      if(input$DAcountDownloadType == "tsv"){
        paste(input$DAlinSelec,extension,".txt", sep = "")
      }
      else{  paste(input$DAlinSelec,extension , ".csv", sep = "")}
    },
    content = function(file){
      selected <- DA_TotalCounts()
      if(input$DAcountDownloadType == "tsv"){
        write_tsv(selected, file )
      }
      else if(input$DAcountDownloadType == "csv"){
        write.csv(selected,file)
      }
    }
  )


  ## GC lin table download
  output$GCdownloadCounts <- downloadHandler(
    req(credentials()$user_auth),
    filename = function(){
      extension <- "-gc_lin_counts"
      if(input$GCcountDownloadType == "tsv"){
        paste(input$GClinSelec,extension,".txt", sep = "")
      }
      else{  paste(input$GClinSelec,extension , ".csv", sep = "")}
    },
    content = function(file){

      selected <- GC_TotalCounts()
      if(input$GCcountDownloadType == "tsv"){
        write_tsv(selected, file )
      }
      else if(input$GCcountDownloadType == "csv"){
        write.csv(selected,file)
      }
    }
  )


}

#Call to shiny app
shinyApp(ui = ui, server = server)

