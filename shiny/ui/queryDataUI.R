# Datatable tab contains all protein datatables
tabPanel( title = "Data/Query",
  value = "datatable",
  loginUI("login"),
        fluidPage(
          tabsetPanel(id = "QueryData",
                      tabPanel("Data Table", value = "mainData",
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
                                 DT::dataTableOutput(outputId = "proTable"), width = 12)
                      ),


                      ## Heatmap by query tab
                      tabPanel("Query Heatmap",value = "queryHeatmapTab",
                               column(width = 12,
                                      plotOutput(outputId = "queryHeatmap", height = "600px")
                               )

                      )
          )


        ))