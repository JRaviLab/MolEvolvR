#### UI for the Domain Architecture tab ####
#### Contains plots and tables pertainent to Domain Architecture

tabItem("domainArchitecture",
        fluidRow(
          column(width = 3, offset = 0,
                 fluidRow(
                   sidebarPanel(
                     width = 12,
                     useShinyjs(),
                     extendShinyjs(text = jscode),
                     #dropdown to select protein for plots
                     selectInput(inputId =  "linSelec", label = "Protein",
                                 choices = c("All", "PspA-Snf7", "PspB", "PspC","PspN", "LiaI-LiaF-TM","Toast-rack")
                                 , selected = "PspA"),

                     column(width = 9, offset = 3,
                            fluidRow(
                              actionButton(inputId = "CutoffSwitch", label = tags$b("Row Cutoff")),
                            )
                     ),
                     sliderInput(inputId = "cutoff", label = "Percent Cutoff:", min = 1, max = 100, value = 95),
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
                              plotOutput(outputId = "LinPlot", height = "650px", width = "100%")),
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
                                box(width = 12,height = '100%',
                                    plotOutput(outputId = "network", width = '100%')
                                ),
                                box(width = 12,
                                    wordcloud2Output(outputId = "wordcloud")
                                    # plotOutput(outputId = "wordcloud")
                                )
                              )
                     )
                   )
                 )
          )
        )
)