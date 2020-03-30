### Genomic Context tab
tabItem(tabName = "GCPlots",
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
                     sliderInput(inputId = "GCcutoff", label = "Total Count Cutoff:", min = 0, max = 100, value = 95),
                     textOutput("Legend")
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
                     tabPanel("Heatmap",plotOutput(outputId = "GCLinPlot", height = '500px' )),
                     tabPanel("Table", DT::dataTableOutput(outputId = "GCLinTable"),
                              column(downloadButton(outputId = "GCDownloadCounts", label = "Download"),radioButtons(inputId = "GCcountDownloadType", label = "Download Type:",
                                                                                                                  choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                     tabPanel("Upset Plot",plotOutput(outputId = "GCUpsetP")),
                     tabPanel("Network",
                              fluidRow(
                                box(width = 12,
                                    plotOutput(outputId = "GCNetwork")
                                ),
                                box(width = 12,
                                    plotOutput(outputId = "GCwordcloud")
                                )
                              )

                     )
                   )
                 )
          )
        )
)