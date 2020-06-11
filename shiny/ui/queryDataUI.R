# Datatable tab contains all protein datatables
tabPanel( title = "Data/Query",
  value = "datatable",
        fluidPage(
          tabsetPanel(id = "QueryData",
                      tabPanel("Data Table", value = "mainData",
                               column(width = 2,
                                      #Dropdown to select protein for viewing
                                      selectInput(inputId =  "proSelec", label = "Protein",
                                                  choices = c("All", "PspA-Snf7", "PspB","PspC", "LiaI-LiaF-TM",
                                                               "Toast-rack", "PspM", "PspN","Tfu-1009",
                                                               "DUF1700", "DUF1707"), selected = "PspA-Snf7")

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