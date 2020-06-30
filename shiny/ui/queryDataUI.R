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
                                                               "DUF1700", "DUF1707"), selected = "All")

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
                               fluidRow( tags$span(
                                 div(class = "queryLegend", style = "margin: 20px;",
                                 strong("The phyletic spread of query proteins."),
                                 p(" The color gradient in the heatmap represents the number of homologs within each lineage. Rows: Query proteins/domains. Columns: Key lineages within the three kingdoms of life.")
                                 )
                                 )),
                               column(width = 12,
                                      plotOutput(outputId = "queryHeatmap", height = "600px")
                               )

                      )
          )


        ))