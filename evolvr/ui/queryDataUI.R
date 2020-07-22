# Datatable tab contains all protein datatables
tabPanel( title = "Data",
          icon = icon("table", lib = "font-awesome"),
          value = "datatable",
          fluidPage(
            tabsetPanel(id = "QueryData",
                        tabPanel("Data Table", value = "mainData",

                                 uiOutput(outputId = "dataTableComponent")

                        ),
                        tabPanel(
                          "FASTA",
                          value = "fastaData",
                          fluidRow(
                            column(width = 12,
                              uiOutput(outputId = "fastaDataComponent")
                            )
                          ),
                          fluidRow(
                            column(width = 12,
                                   uiOutput(outputId = "msaDataComponent")

                                   )
                          )

                        )
            )

          ))
