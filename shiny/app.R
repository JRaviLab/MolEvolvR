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
source("R/GC_network_directed.R")
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
  stringsAsFactors = FALSE
)

###
#### Header ####
###
header <- dashboardHeader(title = "PSP Data",
                          ### Github Icon in header
                          tags$li(
                            a(href = "https://github.com/jananiravi/psp-evolution",
                              icon("github-square",
                                   "fa-2x",
                                   lib = "font-awesome")
                            ),
                            class = "dropdown"
                          )
)

###
#### Sidebar ####
###
sidebar<- dashboardSidebar(
  width = 180,

  sidebarMenu(id = "mainTabs",
              menuItem("Data Table", tabName = "datatable"),
              #menuItem("Lineage Plots", tabName = "lineagePlots"),
              menuItem("Domain Architecture", tabName = "domainArchitecture"),
              menuItem("Genomic Contect", tabName = "genomicContext"),
              menuItem("Phylogeny", tabName = "phylogeny"),
              #menuItem("Usage", tabName="usage"),
              menuItem("About",tabName="about")
  )
)

####
#### Body ####
####
body <- dashboardBody(
  # must turn shinyjs on
  shinyjs::useShinyjs(),
  # add logout button UI
  div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
  # add login panel UI function
  shinyauthr::loginUI(id = "login"),

  uiOutput("user_table"),
  uiOutput("testUI"),
  tabItems(
    # Source main data tab and query heatmap
    source("shiny/UI/queryDataUI.R")$value,



    #lineagePlots contains a heatmap, datatable, and upset plot
    ### Load Lineage Tab ###
    # source("shiny/ui/lineagePlots.R")$value,

    ### Load Domain Architecture Tab ###
    tabItem(tabName = "domainArchitecture",
            fluidRow(
              column(width = 3, offset = 0,
                     fluidRow(
                       sidebarPanel(
                         width = 12,
                         #dropdown to select protein for plots
                         selectInput(inputId =  "DAlinSelec", label = "Protein",
                                     choices = c("All", "PspA-Snf7", "PspB", "PspC","PspN", "LiaI-LiaF-TM","Toast-rack")
                                     , selected = "PspA"),


                         column(width = 9, offset = 3,
                                fluidRow(
                                  actionButton(inputId = "DACutoffSwitch", label = tags$b("Row Cutoff")),
                                )
                         ),

                         #Slider input to determine cutoff value for totalcounts
                         sliderInput(inputId = "DA_Cutoff", label = "Total Count Cutoff:", min = 0, max = 100, value = 95),


                         textOutput("DALegend")
                       )
                     )
              )
              ,


              column(width = 9, offset = 0,
                     #mainpanel dictates what is displayed on screen depending on which tabset panel is selected
                     mainPanel(
                       width = 12,
                       tabsetPanel(
                         id= 'DALin_data',
                         tabPanel("Heatmap", value = "Heatmap",
                                  plotOutput(outputId = "DALinPlot", height = '600px' )),
                         tabPanel("Table", value = "LinTable",
                                  DT::dataTableOutput(outputId = "DALinTable"),
                                  column(downloadButton(outputId = "DAdownloadCounts", label = "Download"),radioButtons(inputId = "DAcountDownloadType", label = "Download Type:",
                                                                                                                        choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                         tabPanel("Upset Plot", value = "Upset",
                                  plotOutput(outputId = "DAUpsetP", height = '600px')),

                         tabPanel("Network",
                                  value = "Network_WC",
                                  fluidRow(
                                    box(width = 12,
                                        plotOutput(outputId = "DANetwork")
                                    ),
                                    box(width = 12,
                                        #plotOutput(outputId = "DAwordcloud")
                                        wordcloud2Output(outputId = "DAwordcloud")
                                    )
                                  )


                         )

                       )
                     )
              )
            )
    ),
    ### Genomic Context tab
    tabItem(tabName = "genomicContext",
            fluidRow(
              column(width = 3, offset = 0,
                     fluidRow(
                       sidebarPanel(
                         width = 12,
                         #dropdown to select protein for plots
                         selectInput(inputId =  "GClinSelec", label = "Protein",
                                     choices = c("All", "PspA-Snf7", "PspB", "PspC","PspN", "LiaI-LiaF-TM","Toast-rack")
                                     , selected = "PspA"),
                         #Slider input to determine cutoff value for totalcounts
                         sliderInput(inputId = "GC_Cutoff", label = "Total Count Cutoff:", min = 0, max = 100, value = 95),
                         textOutput("GCLegend")
                       )
                     )
              )
              ,
              column(width = 9, offset = 0,
                     #mainpanel dictates what is displayed on screen depending on which tabset panel is selected
                     mainPanel(
                       width = 12,
                       tabsetPanel(
                         id= 'GCLin_data',
                         tabPanel("Heatmap", value = "Heatmap",
                                  plotOutput(outputId = "GCLinPlot", height = '600px' )),
                         tabPanel("Table", value = "LinTable",
                                  DT::dataTableOutput(outputId = "GCLinTable"),
                                  column(downloadButton(outputId = "GCDownloadCounts", label = "Download"),radioButtons(inputId = "GCcountDownloadType", label = "Download Type:",
                                                                                                                        choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                         tabPanel("Upset Plot", value = "Upset",
                                  plotOutput(outputId = "GCUpsetP", height = '600px')),
                         tabPanel("WordCloud", value = "Network_WC",
                                  fluidRow(
                                    # box(width = 12,
                                    #     plotOutput(outputId = "GCNetwork")
                                    # ),
                                    box(width = 12,
                                        #plotOutput(outputId = "GCwordcloud")
                                        wordcloud2Output(outputId = "GCwordcloud")
                                    )
                                  )

                         )
                       )
                     )
              )
            )
    ),

    source("shiny/ui/phylogenyUI.R")$value,

    ### Load About Tab ###
    source("shiny/ui/aboutUI.R")$value
  ),
  HTML('<div data-iframe-height></div>')
)
####
#### UI ####
####
ui <- dashboardPage(
  skin = "green",
  #App title
  header,
  #Create Sidebar with inputs and download button output
  sidebar,
  body
)

####
#### Server ####
server <- function(input, output,session){
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

  observe({
    if(credentials()$user_auth) {
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    } else {
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
  })


  #Reactive expression used to determine which data table should be displayed
  #based on selected input
  pspTable<- reactive({
    req(credentials()$user_auth)
    switch(input$proSelec,
           "All" = all %>% select(viewing_cols) %>% distinct(),
           "DUF1700" = DUF1700%>% select(viewing_cols) %>% distinct(),
           "DUF1707" = DUF1707%>% select(viewing_cols) %>% distinct(),
           "PspA-Snf7" = pspa%>% select(viewing_cols) %>% distinct(),
           "Psp-AA" = psp_aa%>% select(viewing_cols) %>% distinct(),
           "PspB" = pspb%>% select(viewing_cols) %>% distinct(),
           "PspC" = pspc%>% select(viewing_cols) %>% distinct(),
           "PspM" = pspm%>% select(viewing_cols) %>% distinct(),
           "PspN" = pspn%>% select(viewing_cols) %>% distinct(),
           "LiaI-LiaF-TM" = liai_liaf%>% select(viewing_cols) %>% distinct(),
           "Toast-rack" = toast_rack%>% select(viewing_cols) %>% distinct(),
           "Tfu-1009" = tfu_1009%>% select(viewing_cols) %>% distinct()
    )
  })
  #Render the Data table for selected protein
  output$proTable <- DT::renderDT({
    req(credentials()$user_auth)
    ##### Can use select here to determine what columns shown
    paged_table(pspTable() )}, extensions = c('FixedColumns',"FixedHeader"),
    options = list(pageLength = 100,
                   #The below line seems to disable other pages and the search bar
                   #dom = 't',
                   scrollX = TRUE,
                   paging=TRUE,
                   fixedHeader=TRUE,
                   fixedColumns = list(leftColumns = 2, rightColumns = 0)))

  #### Query Heatmap ####
  output$queryHeatmap <- renderPlot({
    req(credentials()$user_auth)
    lineage.Query.plot(query_data = all, queries = queries, colname = "DomArch", cutoff = 100)
  })



  #### Reactive expression determining data for which plotting is based. Controlled by dropdown
  plotting_prot <-  reactive({
    if(input$mainTabs == "domainArchitecture"){
      switch(input$DAlinSelec,
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
             "Tfu-1009" = tfu_1009)
    }
    # else if(input$mainTabs == "genomicContext")
    # {
    else{
      switch(input$GClinSelec,
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
             "Tfu-1009" = tfu_1009)
    }
  })

  # prot_word_percents <- reactive({
  #   if(input$mainTabs == "domainArchitecture"){
  #     #### Should these be on DomArch.repeats for Network?
  #     max_word_percents(plotting_prot(), "DomArch")
  #   }
  #   # else if(input$mainTabs == "genomicContext"){
  #   #   max_word_percents(plotting_prot(), "GenContext")
  #   # }
  # })


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

  # ### Adaptive cutoff value to be used for plotting functions
  # ### Converts the cutoff slider's value into respective percentage value
  # ### if the cutoff status is "Row"
  # cutoff_val <- reactive({
  #   if(cutoff_status() == "Percent")
  #   {
  #     # Cutoff type is Percentage already, so use original slider value
  #     if(input$mainTabs == "domainArchitecture"){
  #       input$DA_cutoff
  #     }
  #     else if(input$mainTabs == "genomicContext")
  #     {
  #       input$GC_Cutoff
  #     }
  #   }
  #   else
  #   {
  #     # Cutoff type is Row, convert it to the percentage
  #     if(input$lin_data == "Heatmap")
  #     {
  #       if(input$linSelec == "All")
  #       {
  #         if(input$DA_GC== "Domain Architecture"){
  #           100 - query_DA_row_CutoffPercs$maxPercent[input$cutoff]
  #         }
  #         else
  #         {
  #           100 - query_GC_row_CutoffPercs$maxPercent[input$cutoff]
  #         }
  #       }
  #       else{
  #         if(input$DA_GC== "Domain Architecture"){
  #           column = "DomArch"
  #           rmAstrk = F
  #         }
  #         else
  #         {
  #           column = "GenContext"
  #           rmAstrk = T
  #         }
  #         rownumber_to_cutoff(plotting_prot(),input$cutoff, col = column)
  #       }
  #     }
  #
  #     else
  #     {
  #       # Cutoff to percent for upset and wordclouds and network
  #       100 - prot_word_percents()$MaxPercent[input$cutoff]
  #     }
  #   }
  # })

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


  # output$LinPlot <- renderPlot({
  #   req(credentials()$user_auth)
  #   if(input$DA_GC == "Domain Architecture"){
  #     if(input$linSelec != "All"){
  #       lineage.DA.plot(plotting_prot(), colname = "DomArch", cutoff = cutoff_val(), RowsCutoff = rows_cutoff())
  #     }
  #     else{
  #       lineage.Query.plot(plotting_prot(), queries = queries, colname = "DomArch", cutoff = cutoff_val())
  #     }
  #   }
  #   else{
  #     if(input$linSelec != "All"){
  #       lineage.DA.plot(plotting_prot(), colname = "GenContext", cutoff = cutoff_val(), RowsCutoff = rows_cutoff())
  #     }
  #     else{
  #       lineage.Query.plot(plotting_prot(), queries = queries, colname = "GenContext", cutoff = cutoff_val())
  #     }
  #   }
  # }, height = "auto")

  ####
  ##### Render the datatable for the lineage counts #####
  ####
  output$DALinTable <- DT::renderDT({
    req(credentials()$user_auth)
    paged_table(
      total_counts(plotting_prot(), cutoff = DA_cutoff_val(),column = "DomArch")
    )
  },
  extensions = c('FixedColumns',"FixedHeader"),
  options = list(pageLength = 15,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))

  output$GCLinTable <- DT::renderDT({
    req(credentials()$user_auth)
    paged_table(
      total_counts(plotting_prot(), cutoff = GC_cutoff_val(), column = "GenContext")
    )

  },
  extensions = c('FixedColumns',"FixedHeader"),
  options = list(pageLength = 15,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))


  # output$LinTable <- DT::renderDT({
  #   req(credentials()$user_auth)
  #   paged_table(
  #     if(input$DA_GC == "Domain Architecture"){
  #       total_counts(plotting_prot(), cutoff = cutoff_val(), column = "DomArch", RowsCutoff = rows_cutoff()) # %>% arrange(CumulativePercent)
  #     }
  #     else{
  #       total_counts(plotting_prot(), cutoff = cutoff_val(), column = "GenContext", RowsCutoff = rows_cutoff()) # %>% arrange(CumulativePercent)
  #     }
  #   )
  # }, extensions = c('FixedColumns',"FixedHeader"),
  # options = list(pageLength = 15,
  #                scrollX = TRUE,
  #                paging=TRUE,
  #                fixedHeader=TRUE,
  #                fixedColumns = list(leftColumns = 2, rightColumns = 0)))



  #### Upset Plots ####
  output$DAUpsetP <-renderPlot({
    req(credentials()$user_auth)
    upset.plot(plotting_prot(), cutoff = DA_cutoff_val(), colname = "DomArch", RowsCutoff = rows_cutoff())
  })

  output$GCUpsetP <- renderPlot({
    req(credentials()$user_auth)
    upset.plot(plotting_prot(), cutoff = GC_cutoff_val(), colname = "GenContext", RowsCutoff = rows_cutoff())
  })

  #   output$upsetP <- renderPlot({
  #     req(credentials()$user_auth)
  #     selected <- input$linSelec
  #     DA_or_GC <- input$DA_GC
  #     if(DA_or_GC == "Domain Architecture"){
  #       upset.plot(plotting_prot(), cutoff = cutoff_val(), colname = "DomArch", RowsCutoff = rows_cutoff())
  #     }
  #     else{
  #       upset.plot(plotting_prot(), cutoff = cutoff_val(), colname = "GenContext", RowsCutoff = rows_cutoff())
  #     }
  #   }, height = 550)


  #### Reactive expression determining domain of interest for plotting domain networks
  network_domain_interest <-  reactive({
    switch(input$DAlinSelec,
           "All" = c("DUF1700-ahelical","DUF1707-SHOCT","PspA", "Snf7","PspB","PspC","PspM","PspN","LiaI-LiaF-TM","Toast-rack"),
           "DUF1700" = "DUF1700-ahelical",
           "DUF1707" = "DUF1707-SHOCT",
           "PspA-Snf7" = c("PspA","Snf7"),
           "PspB" = "PspB",
           "PspC" = "PspC",
           "PspM" = "PspM",
           "PspN" = "PspN",
           "LiaI-LiaF-TM" = "LiaI-LiaF-TM",
           "Toast-rack" = "Toast-rack")
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

  # output$wordcloud <- renderPlot({
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
    lineage_sunburst(phylogeny_prot(), lineage_column = "Lineage", type = "sunburst")
  })


  output$msaPlot <- renderUI({
    req(credentials()$user_auth)
    tags$iframe(style="height:600px; width:100%", src="pspa_reduced.fasta.pdf", seamless=T)
  })





  #### Paralogs ####

  output$ParalogTable <- DT::renderDataTable({
    req(credentials()$user_auth)
    find_paralogs(phylogeny_prot())
  },extensions = c('FixedColumns'),
  options = list(pageLength = 10,
                 #The below line seems to disable other pages and the search bar
                 #dom = 't',
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))






  #Reactive expresion to change file name depending on which protein is selected
  fileNam <- reactive({
    if(input$downloadType == "tsv"){
      paste(input$proSelec, ".txt", sep = "")
    }
    else  paste(input$proSelec, ".csv", sep = "")
  })

  ###### Downloads the data from datatable #####
  ### Should it download cleaned data? or all rows?
  output$downloadData <- downloadHandler(
    filename = function() {
      fileNam()
    },
    content = function(file) {

      if(input$downloadType == "tsv"){
        write_tsv(pspTable(), file)
      }
      else  write.csv(pspTable(), file, row.names = FALSE)
    }
  )

  #####
  ##### Download Data for Lin Table #####
  #####
  output$downloadCounts <- downloadHandler(
    req(credentials()$user_auth),
    filename = function(){
      if(input$DA_GC == "Domain Architecture"){extension <- "-da_lin_counts"}
      else if(input$DA_GC == "Genomic Context"){ extension <- "-gc_lin_counts"}
      if(input$downloadType == "tsv"){
        paste(input$linSelec,extension,".txt", sep = "")
      }
      else{  paste(input$linSelec,extension , ".csv", sep = "")}
    },
    content = function(file){
      if(input$DA_GC == "Domain Architecture"){
        selected <- switch(input$linSelec,
                           "PspA" = pspa_totalC,
                           "PspB" = pspb_totalC,
                           "PspC" = pspc_totalc,
                           "LiaF" = liaf_totalC,
                           "LiaG" = liag_totalC,
                           "LiaI" = liai_totalC)
      }
      else if(input$DA_GC == "Genomic Context"){
        selected <- switch(input$linSelec,
                           "PspA" = pspa_cum,
                           "PspB" = pspb_cum,
                           "PspC" = pspc_cum,
                           "LiaF" = liaf_cum,
                           "LiaG" = liag_cum,
                           "LiaI" = liai_cum)
      }
      if(input$countDownloadType == "tsv"){
        write_tsv(selected, file )
      }
      else if(input$countDownloadType == "csv"){
        write.csv(selected,file)
      }
    }
  )
}

#Call to shiny app
shinyApp(ui = ui, server = server)

