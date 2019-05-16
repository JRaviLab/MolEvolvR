library(shiny)
library(tidyverse)
library(DT)
library(rmarkdown)
library(shinydashboard)
library(shinyjqui)

source("PSP_Web_Data.R")

source("shinyfunctions.R")

#########
#Sidebar#
#########
sidebar<- dashboardSidebar(
  sidebarMenu(id = "tabs",
              menuItem("Data Table", tabName = "datatable"),
              menuItem("Lineage Plots", tabName = "lineagePlots")
              #              ,
              #              menuItem("Phylogeny", tabName = "phylogeny")
  )
)
######
#Body#
######
body <- dashboardBody(
  tabItems(
    #Datatable tab contains all protein datatables
    tabItem("datatable",
            fluidPage(
              column(width = 2,
                     #Dropdown to select protein for viewing
                     selectInput(inputId =  "proSelec", label = "Protein",
                                 choices = c( #"DUF1700", "DUF1707",
                                   "PspA", "PspB", "PspC", "PspM", "PspN","Liai","Liaf", "Liag")
                                 , selected = "PspA")
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
    tabItem("lineagePlots",
            fluidPage(
              sidebarLayout(
                sidebarPanel(
                  #dropdown to select protein
                  selectInput(inputId =  "linSelec", label = "Protein",
                              choices = c( "PspA", "PspB", "PspC","Liaf","Liag","Liai")
                              , selected = "PspA"),
                  #Radiobuttons to selext domain architecture and genomic context
                  radioButtons(inputId = "DA_GC", label = "Lineage by:"
                               , choices= c("Domain Architecture", "Genomic Context")
                               , selected = "Domain Architecture"),
                  #Slider input to determine cutoff value for totalcounts
                  sliderInput(inputId = "cutoff", label = "Total Count Cutoff:", min = 0, max = 500, value = 30)
                ),
                #mainpanel dictates what is displayed on screen depending on which tabset panel is selected
                mainPanel(
                  tabsetPanel(
                    id= 'lin_data',
                    tabPanel("Heatmap",jqui_resizable(plotOutput(outputId = "LinPlot", height = '500px' ))),
                    tabPanel("Table", DT::dataTableOutput(outputId = "LinTable")),
                    tabPanel("Upset Plot", jqui_resizable(plotOutput(outputId = "upsetP")))
                  )
                )


              )
            )
    )#,
    #    tabItem("phylogeny",
    #            fluidPage(sidebarLayout(
    #              sidebarPanel(
    #                selectInput(inputId =  "alignSelec", label = "Protein",
    #                            choices = c( "PspA", "PspB", "PspC","PspM","PspN")
    #                            , selected = "PspN"),
    #                radioButtons(inputId = "plottype", label = "Plot Type",
    #                            choices = c("Tree1", "Tree2", "Tree3"),selected = "Tree1")),
    #              mainPanel(
    #                tabsetPanel(
    #                  id= 'lin_data',
    #                  tabPanel("Align1",plotOutput(outputId = "msaPlot" ))
    #    tabPanel("Table", DT::dataTableOutput(outputId = "LinTable")),
    #   tabPanel("Upset Plot", plotOutput(outputId = "upsetP"))
    #                )
    #              )
    #            )
    #           )
    #    )
  )
)
####
#UI#
####
ui <- dashboardPage(skin = "green",
                    #App title
                    dashboardHeader(title = "PSP Data"),
                    #Create Sidebar with inputs and download button output
                    sidebar,
                    body
)

########
#Server#
########
server <- function(input, output){


  #Reactive expression used to determine which data table should be displayed
  #based on selected input
  pspTable<- reactive({
    switch(input$proSelec,
           "DUF1700" = DUF1700_table,
           "DUF1707" = DUF1707_table,
           "PspA" = pspa_table,
           "PspB" = pspb_table,
           "PspC" = pspc_table,
           "PspM" = pspm_table,
           "PspN" = pspn_table,
           "Liai" = liai_table,
           "Liaf"= liaf_table,
           "Liag" = liag_table)
  })
  #Render the Data table for selected protein
  output$proTable <- DT::renderDT({paged_table(pspTable())}, extensions = c('FixedColumns',"FixedHeader"),
                                  options = list(pageLength = 10,
                                                 #The below line seems to disable other pages and the search bar
                                                 #dom = 't',
                                                 scrollX = TRUE,
                                                 paging=TRUE,
                                                 fixedHeader=TRUE,
                                                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))

  #Renders the heatmap
  output$LinPlot <- renderPlot({
    if(input$DA_GC == "Domain Architecture"){
      switch(input$linSelec,
             "PspA" = lineage.DA.plot(pspa_DA_Lin, pspa_totalC,type = "da2doms",cutoff = input$cutoff),
             "PspB" = lineage.DA.plot(pspb_DA_Lin, pspb_totalC,type ="da2doms", cutoff = input$cutoff),
             "PspC" = lineage.DA.plot(pspc_DA_Lin, pspc_totalC,type = "da2doms",cutoff = input$cutoff),
             "Liaf" = lineage.DA.plot(liaf_DA_lin, liaf_totalC, type = "da2doms", cutoff = input$cutoff),
             "Liag" = lineage.DA.plot(liag_DA_lin, liag_totalC, type = "da2doms", cutoff = input$cutoff),
             "Liai" = lineage.DA.plot(liai_DA_lin, liai_totalC, type = "da2doms", cutoff = input$cutoff))
    }
    else{
      switch(input$linSelec,
             "PspA" = lineage.DA.plot(pspa_cum, pspa_cum,type = "gc2da", cutoff =input$cutoff),
             "PspB" = lineage.DA.plot(pspb_cum, pspb_cum,type = "gc2da",cutoff =input$cutoff),
             "PspC" = lineage.DA.plot(pspc_cum, pspc_cum,type = "gc2da",cutoff = input$cutoff),
             "Liaf" = lineage.DA.plot(liaf_cum, liaf_cum, type = "gc2da", cutoff = input$cutoff),
             "Liag" = lineage.DA.plot(liag_cum, liag_cum, type = "gc2da", cutoff = input$cutoff),
             "Liai" = lineage.DA.plot(liai_cum, liai_cum, type = "gc2da", cutoff = input$cutoff))
    }
  }, height = 500)

  #Renders the datatable for the lineage counts
  output$LinTable <- DT::renderDT({paged_table(
    if(input$DA_GC == "Domain Architecture"){
      switch(input$linSelec,
             "PspA" = filter(pspa_totalC,totalcount >= input$cutoff),
             "PspB" = filter(pspb_totalC,totalcount >= input$cutoff),
             "PspC" = filter(pspc_totalC,totalcount >= input$cutoff),
             "Liaf" = filter(liaf_totalC, totalcount >= input$cutoff),
             "Liai" = filter(liai_totalC, totalcount >= input$cutoff),
             "Liag" = filter(liag_totalC, totalcount >= input$cutoff)
      )}
    else{
      switch(input$linSelec,
             "PspA" = filter(pspa_cum,totalcount >= input$cutoff),
             "PspB" = filter(pspb_cum,totalcount >= input$cutoff),
             "PspC" = filter(pspc_cum,totalcount >= input$cutoff),
             "Liaf" = filter(liaf_cum,totalcount >= input$cutoff),
             "Liag" = filter(liag_cum,totalcount >= input$cutoff),
             "Liai" = filter(liai_cum,totalcount >= input$cutoff)
      )
    }
  )
  }, extensions = c('FixedColumns',"FixedHeader"),
  options = list(pageLength = 15,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))


  #Renders the upsetPlot
  output$upsetP <- renderPlot({
    selected <- input$linSelec
    DA_or_GC <- input$DA_GC
    if(DA_or_GC == "Domain Architecture"){
      switch(selected,
             "PspA"= upset.plot(pspa_table,pspa.DA.doms.wc,input$cutoff, "da2doms"),
             "PspB"= upset.plot(pspb_table,pspb.DA.doms.wc,input$cutoff, "da2doms"),
             "PspC"= upset.plot(pspc_table,pspc.DA.doms.wc, input$cutoff, "da2doms"),
             "Liaf"= upset.plot(liaf_table,liaf.DA.doms.wc,input$cutoff, "da2doms"),
             "Liag"= upset.plot(liag_table,liag.DA.doms.wc,input$cutoff, "da2doms"),
             "Liai"= upset.plot(liai_table,liai.DA.doms.wc,input$cutoff, "da2doms")
      )
    }
    else{
      switch(selected,
             "PspA"= upset.plot(pspa_table, pspa.GC.doms.wc, input$cutoff, "gc2da"),
             "PspB"= upset.plot(pspb_table, pspb.GC.doms.wc, input$cutoff, "gc2da"),
             "PspC"= upset.plot(pspc_table, pspc.GC.doms.wc, input$cutoff, "gc2da"),
             "Liaf"=upset.plot(liaf_table, liaf.GC.doms.wc, input$cutoff, "gc2da"),
             "Liag"=upset.plot(liag_table, liag.GC.doms.wc, input$cutoff, "gc2da"),
             "Liai"=upset.plot(liai_table, liai.GC.doms.wc, input$cutoff, "gc2da")
      )
    }
  }, height = 550)

  msaPlotType <-reactive({
    switch(input$plottype,
           "Tree1"= "apeTree",
           "Tree2" = "ggTree",
           "Tree3"= "msaTree")
  })

  #  output$msaPlot <- renderPlot({
  #    switch(input$alignSelec,
  #           "PspN"= phylo.plots("data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt",msaPlotType()))
  #  }, height = 600
  #  )

  #Reactive expresion to change file name depending on which protein is selected
  fileNam <- reactive({
    if(input$downloadType == "tsv"){
      paste(input$proSelec, ".txt", sep = "")
    }
    else  paste(input$proSelec, ".csv", sep = "")
  })

  #Downloads the data
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

}

#Call to shiny app
shinyApp(ui = ui, server = server)


