library(shiny)
library(tidyverse)
library(DT)
library(rmarkdown)
library(shinydashboard)
library(shinyjqui)
library(wordcloud)
library(shinyauthr)
library(svgPanZoom)
conflicted::conflict_prefer("intersect", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("strsplit", "base")
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("box", "shinydashboard")
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

#########
##Users##
#########
user_base <- data.frame(
  user = c("pspevolution"),
  password = c("cpathogeno2019"),
  permissions = c("admin"),
  name = c("User One"),
  stringsAsFactors = FALSE
)

########
#Header#
########
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

#########
#Sidebar#
#########
sidebar<- dashboardSidebar(
  width = 180,

  sidebarMenu(id = "mainTabs",
              menuItem("Data Table", tabName = "datatable"),
              menuItem("Lineage Plots", tabName = "lineagePlots"),
              menuItem("Phylogeny", tabName = "phylogeny")
              #,menuItem("Usage", tabName="usage")
              ,menuItem("About",tabName="about")
  )
)
######
#Body#
######
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
    #Datatable tab contains all protein datatables
    tabItem("datatable",
            fluidPage(
              column(width = 2,
                     #Dropdown to select protein for viewing
                     selectInput(inputId =  "proSelec", label = "Protein",
                                 choices = c( "All","DUF1700", "DUF1707",
                                              "PspA-Snf7","Psp-AA", "PspB", "PspC", "PspM", "PspN",
                                              "LiaI-LiaF-TM","Toast-rack", "Tfu-1009" )
                                 , selected = "All")
              ),
              #Buttons to select which file type to download
              column( width = 3, offset= 1,
                      #Radiobuttons to select what to download data table as: tab separated or comma seperated
                      radioButtons(inputId = "downloadType", label = "Download Type:",
                                   choices= c("tsv", "csv"), selected = "tsv" ),
                      #Output download button
                      downloadButton(outputId = "downloadData", label = "Download")),
              #Create mainpanel where dataTable is displayed
              column(
                DT::dataTableOutput(outputId = "proTable"), width = 12)))
    ,

    #lineagePlots contains a heatmap, datatable, and upset plot
    ### Load Lineage Tab ###
    source("shiny/ui/lineagePlots.R")$value,

    source("shiny/ui/phylogenyUI.R")$value,

    ### Load About Tab ###
    source("shiny/ui/aboutUI.R")$value
  ),
  HTML('<div data-iframe-height></div>')
)
####
#UI#
####
ui <- dashboardPage(
  skin = "green",
  #App title
  header,
  #Create Sidebar with inputs and download button output
  sidebar,
  body
)

########
#Server#
########
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
  #Render the Data table for selected protein
  output$proTable <- DT::renderDT({
    req(credentials()$user_auth)
    ##### Can use select here to determine what columns shown
    paged_table(pspTable() )}, extensions = c('FixedColumns',"FixedHeader"),
    options = list(pageLength = 10,
                   #The below line seems to disable other pages and the search bar
                   #dom = 't',
                   scrollX = TRUE,
                   paging=TRUE,
                   fixedHeader=TRUE,
                   fixedColumns = list(leftColumns = 2, rightColumns = 0)))


  #### Reactive expression determining data for which plotting is based. Controlled by dropdown ####
  plotting_prot <-  reactive({
    switch(input$linSelec,
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
  })

  #Observer used to determine initial heatmap slider
  observe({
    if(input$DA_GC== "Domain Architecture"){
      cutoff_init <- top_n_rows_cutoff(plotting_prot(), column = "DomArch", n_rows = 10)
    }
    else{
      cutoff_init <- top_n_rows_cutoff(plotting_prot(), column = "GenContext", n_rows = 10)
    }
    updateSliderInput(session,inputId = "cutoff",min=0, max=100, value=cutoff_init)
  })



  ####
  ##### Render the heatmap #####
  ####
  output$LinPlot <- renderPlot({
    req(credentials()$user_auth)
    if(input$DA_GC == "Domain Architecture"){
      if(input$linSelec != "All"){
        lineage.DA.plot(plotting_prot(), colname = "DomArch", type = "da2doms", cutoff = input$cutoff)
      }
      else{
        lineage.Query.plot(plotting_prot(), queries = queries, colname = "DomArch", cutoff = input$cutoff)
      }
    }
    else{
      if(input$linSelec != "All"){
        lineage.DA.plot(plotting_prot(), colname = "GenContext", type = "gc2da", cutoff = input$cutoff)
      }
      else{
        lineage.Query.plot(plotting_prot(), queries = queries, colname = "GenContext", cutoff = input$cutoff)
      }
    }
  }, height = 500)

  ####
  ##### Render the datatable for the lineage counts #####
  ####
  output$LinTable <- DT::renderDT({
    req(credentials()$user_auth)
    paged_table(
      if(input$DA_GC == "Domain Architecture"){
        total_counts(plotting_prot(), cutoff = input$cutoff, column = "DomArch")
      }
      else{
        total_counts(plotting_prot(), cutoff = input$cutoff, column = "GenContext")
      }
    )
  }, extensions = c('FixedColumns',"FixedHeader"),
  options = list(pageLength = 15,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))


  # upset_vals <- reactiveValues(GC.DA.wc=pspa$GC.DA %>% words2wc(), DA.doms.wc=pspa$DA.doms %>% words2wc())
  # observeEvent(input$linSelec,{
  #     upset_vals$GC.DA.wc <- switch(input$linSelec,
  #                             "PspA" = pspa$GC.DA %>% words2wc(),
  #                             "PspB" = pspb$GC.DA %>% words2wc(),
  #                             "PspC" = pspc$GC.DA %>% words2wc(),
  #                             "PspN" = pspn$GC.DA %>% words2wc(),
  #                             "LiaF" = liai_liaf$GC.DA %>% words2wc(),
  #                             "Toast-rack" = toast_rack$GC.DA %>% words2wc(),
  #                             "Liag" = liag$GC.DA %>% words2wc())
  #     upset_vals$DA.doms.wc <- switch(input$linSelec,
  #                               "PspA" = pspa$DA.doms %>% words2wc(),
  #                               "PspB" = pspb$DA.doms %>% words2wc(),
  #                               "PspC" = pspc$DA.doms %>% words2wc(),
  #                               "PspN" = pspn$DA.doms %>% words2wc(),
  #                               "LiaF" = liai_liaf$DA.doms %>% words2wc(),
  #                               "Toast-rack" = toast_rack$DA.doms %>% words2wc(),
  #                               "Liag" = liag$DA.doms %>% words2wc())
  # }
  # )

  #Change from: "PspA"= lineage.upset(pspa,input$cutoff, "da2doms",upset_vals$DA.doms.wc)
  # to filtering the data with cutoff of Total count first then just running it through?
  #Renders the upsetPlot
  output$upsetP <- renderPlot({
    req(credentials()$user_auth)
    selected <- input$linSelec
    DA_or_GC <- input$DA_GC
    if(DA_or_GC == "Domain Architecture"){
      upset.plot(plotting_prot(), cutoff = input$cutoff, colname = "DomArch" )
    }
    else{
      upset.plot(plotting_prot(), cutoff = input$cutoff, colname = "GenContext")
    }
  }, height = 550)


  #### Reactive expression determining domain of interest for plotting domain networks
  network_domain_interest <-  reactive({
    switch(input$linSelec,
           "All" = "DUF1700-ahelical|DUF1707-SHOCT|pspa|snf7|pspb|pspc|pspm|pspn|LiaI-LiaF-TM|Toast-rack",
           "DUF1700" = "DUF1700-ahelical",
           "DUF1707" = "DUF1707-SHOCT",
           "PspA-Snf7" = "pspa|snf7",
           "PspB" = "pspb",
           "PspC" = "pspc",
           "PspM" = "pspm",
           "PspN" = "pspn",
           "LiaI-LiaF-TM" = "LiaI-LiaF-TM",
           "Toast-rack" = "Toast-rack")
  })


  ### Network Output ###
  output$network <- renderPlot({
    if(input$DA_GC == "Domain Architecture"){
      domain_network(plotting_prot(), column = "DomArch.repeats", cutoff = input$cutoff, layout = "auto",
                     domains_of_interest = network_domain_interest())
    }
    else{
      gc_directed_network(plotting_prot(), column = "GenContext.repeats",
                          cutoff = input$cutoff)
    }
  })


  #Render Wordcloud
  #reverse_operons probably not necessary anymore
  output$wordcloud <- renderPlot({
    req(credentials()$user_auth)
    if(input$DA_GC == "Genomic Context"){
      wordcloud_element(type = "gc2da", query_data = plotting_prot(), cutoff = input$cutoff)
    }
    else{
      wordcloud_element(type = "da2doms", query_data = plotting_prot(), cutoff = input$cutoff)
    }

  }, height = 550)

  ### Observe changes in tabset panel to determine legend text ###
  legend_txt <- reactive({
    # switch(input$lin_data,
    #   "Upset" = Upset_txt,
    #   "Heatmap" = LineageHeatmap_txt
    # )
    if(input$lin_data == "Heatmap" & input$linSelec != "All"){
      LineageHeatmap_txt
    }
    else if(input$lin_data == "Heatmap" & input$linSelec == "All"){
      AllLineageHeatmap_txt
    }
    else if(input$lin_data == "Upset"){
      Upset_txt
    }
    else if(input$lin_data == "Lintables"){
      LinTable_txt
    }
    else if(input$lin_data == "Network_WC"){
      paste(Network_txt, Wordcloud_txt, sep = "\n")
    }


  })

  output$Legend <- renderText({
    legend_txt()
  })




  observe({
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
                  "PspC")
      selected = "PspA-Snf7"
    }
    updateSelectInput(session, inputId = "alignSelec",
                      choices = choices, selected = selected)
  })

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
             #tags$embed(src="pspa_snf7_tree.pdf"), #, "data-zoom-image" ="pspa_snf7_tree.pdf"),#, height=1024,width=800),
           "PspB Gismo" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspb.gismo_tree.pdf", seamless=T),
           "PspC Gismo" = tags$iframe(style="height:600px; width:100%", src="FigTrees/pspc.gismo_tree.pdf", seamless=T),
           "Snf7 Only" = tags$iframe(style="height:600px; width:100%", src="FigTrees/snf7_only.1_tree.pdf", seamless=T),
           "Toast-rack DUF2154-LiaF" = tags$iframe(style="height:600px; width:100%", src="FigTrees/Toast-rack_DUF2154-LiaF.1_tree.pdf", seamless=T),
           "Toast-rack DUF2807" = tags$iframe(style="height:600px; width:100%", src="FigTrees/Toast-rack_DUF2807.2_tree.pdf", seamless=T),
           "Toast-rack DUF4097" = tags$iframe(style="height:600px; width:100%", src="Toast-rack_DUF4097-LiaG.3.allfa.2nd_tree.pdf", seamless=T),
           "Toast-rack PspC-Cterm" = tags$iframe(style="height:600px; width:100%", src="FigTrees/Toast-Rack_PspC-Cterm.4_tree.pdf", seamless=T)
    )
  })

  output$msaPlot <- renderUI({
    req(credentials()$user_auth)
    tags$iframe(style="height:600px; width:100%", src="pspa_reduced.fasta.pdf", seamless=T)
  })

  vals <- reactiveValues(btn = 0, tab = "home")


  observeEvent(input$myBtn,{
    if( input$phylo =="Tree" ){
      vals$btn <- 1
      vals$tab <- input$phylo
    }
  }
  )
  observeEvent(input$phylo,{
    if( input$phylo !="Tree"){
      vals$tab <- "notphylo"
      vals$btn <- 0
    }
  }
  )


  observe({
    if(vals$btn == 1 && input$mainTabs == "phylogeny"){
      session$sendCustomMessage(type = 'testmessage'
                                ,message = list(value="ZoomOn")
      )
    }
    else{
      session$sendCustomMessage(type = 'testmessage'
                                ,message = list(value="ZoomOff")
      )
    }
  })


  output$ParalogTable <- DT::renderDataTable({
    req(credentials()$user_auth)
    switch(input$alignSelec,
           "PspA-Snf7"= find_paralogs(pspa)
           ,"PspB"= find_paralogs(pspb),
           "PspC"= find_paralogs(pspc)
           # "DUF1700"= find_paralogs(all%>% filter(Query=="DUF1700-alpha-helical")),
           # "DUF1707"= find_paralogs(all%>% filter(Query=="DUF1707-SHOCT")),
           # "Toast-rack"= find_paralogs(all%>% filter(Query=="Toast-rack")),
           # "LiaI-LiaF-TM"= find_paralogs(all%>% filter(Query=="LiaI-LiaF-TM"))
    )
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

