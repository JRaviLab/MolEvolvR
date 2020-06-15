# Datatable tab contains all protein datatables
tabPanel( title = "Data/Query",
  value = "datatable",
        fluidPage(
          tabsetPanel(id = "QueryData",
                      tabPanel("Data Table", value = "mainData",
                               column(width = 2,
                                      #Dropdown to select protein for viewing
                                      selectInput(inputId =  "mainSelect", label = "Protein",
                                                  choices = c("All"), selected = "All")

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
                                 DT::dataTableOutput(outputId = "mainTable"), width = 12)
                      ),


                      ## Heatmap by query tab
                      tabPanel("Query Heatmap",value = "queryHeatmapTab",
                               column(width = 12,
                                      plotOutput(outputId = "queryHeatmap", height = "600px")
                               )

                      )
          )


        ))