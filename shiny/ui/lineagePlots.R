# UI component containing Lineage Plots


tabItem("lineagePlots",
        fluidRow(
          column(width = 3, offset = 0,
                 fluidRow(
                   sidebarPanel(
                     width = 12,
                     #dropdown to select protein for plots
                     selectInput(inputId =  "linSelec", label = "Protein",
                                 choices = c("All", "PspA-Snf7", "PspB", "PspC","PspN", "LiaI-LiaF-TM","Toast-rack")
                                 , selected = "PspA"),

                     radioButtons(inputId = "DA_GC", label = "DA or GC:",
                                  choices= c("Domain Architecture", "Genomic Context"), selected = "Domain Architecture" ),
                     #Slider input to determine cutoff value for totalcounts
                     sliderInput(inputId = "cutoff", label = "Total Count Cutoff:", min = 0, max = 100, value = 95),
                     #sliderInput(inputId = "nr", label="Rows",min=0,max=10,value=0)


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
                     id= 'lin_data',
                     tabPanel("Heatmap",
                              value = "Heatmap",
                              plotOutput(outputId = "LinPlot", height = '500px' )),
                     tabPanel("Table",
                              value = "Lintables",
                              DT::dataTableOutput(outputId = "LinTable"),
                              column(downloadButton(outputId = "downloadCounts", label = "Download"),radioButtons(inputId = "countDownloadType", label = "Download Type:",
                                                                                                                  choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                     tabPanel("Upset Plot",
                              value = "Upset",
                              plotOutput(outputId = "upsetP")),
                     tabPanel("Network",
                              value = "Network_WC",
                              fluidRow(
                                box(width = 12,
                                    plotOutput(outputId = "network")
                                ),
                                box(width = 12,
                                    plotOutput(outputId = "wordcloud")
                                )
                              )
                     )
                   )
                 )
          )
        )
)