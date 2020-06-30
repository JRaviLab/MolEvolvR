# Domain Architecture Panel
tabPanel( title = "Domain Architecture",
          value = "domainArchitecture",
          fluidRow(

            #column(width = 10, offset = 1,
            box( width = 12,

                 column(width = 4,
                        #dropdown to select protein for plots
                        selectInput(inputId =  "DAlinSelec", label = "Protein",
                                    choices = c("All", "PspA","Snf7", "PspB", "PspC", "LiaI-LiaF-TM",
                                                "Toast-rack","PspM","PspN", "DUF1700-ahelical",
                                                "DUF1707-SHOCT", "Tfu-1009")
                                    , selected = "PspA", width = "100%")
                 ),

                 column(width = 4,
                        actionButton(inputId = "DACutoffSwitch", label = tags$b("Row Cutoff")),


                        #Slider input to determine cutoff value for totalcounts
                        sliderInput(inputId = "DA_Cutoff", label = "Total Count Cutoff:", min = 0, max = 100, value = 95)
                 ),

                 column(width = 4,
                        textOutput("DALegend", inline = T)
                 )

            )



            ,


            column(width = 12, offset = 0,
                   #mainpanel dictates what is displayed on screen depending on which tabset panel is selected
                   mainPanel(
                     width = 12,
                     tabsetPanel(
                       id= 'DALin_data',
                       tabPanel("Heatmap", value = "Heatmap",
                                selectInput(inputId = "DA_lin_color", label = "Theme", choices =
                                              c("default", "viridis", "inferno","magma", "plasma", "cividis"), selected = "default"
                                            ),
                                plotOutput(outputId = "DALinPlot", height = '600px' )),
                       tabPanel("Table", value = "LinTable",
                                p("Select a row to see which lineages the domain architecture is present in", style = "color:#242320;", class = "note-box"),
                                DT::dataTableOutput(outputId = "DALinTable"),
                                column(downloadButton(outputId = "DAdownloadCounts", label = "Download"),radioButtons(inputId = "DAcountDownloadType", label = "Download Type:",
                                                                                                                      choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                       tabPanel("Upset Plot", value = "Upset",
                                plotOutput(outputId = "DAUpsetP", height = '600px')),

                       tabPanel("Network",
                                value = "Network_WC",
                                fluidRow(
                                  tags$div(class = "note-box",
                                           "Messy graph? Try re-arranging the vertices by
                                  clicking and dragging the vertices! Also try zooming in and out using your scroll wheel!"
                                  ),
                                  tags$div(class = "bord",
                                           tags$div(class = "innerbox",
                                                    selectInput(inputId = "networkLayout", label = "Layout:",
                                                                choices = c("nice", "grid", "circle", "random"),
                                                                selected = "nice"),
                                                    visNetworkOutput(outputId = "DANetwork")
                                                    #plotOutput(outputId = "DANetwork")
                                           )
                                  ),
                                  tags$div( class = "bord",
                                            tags$div(class = "innerbox",
                                                     # plotOutput(outputId = "DAwordcloud")
                                                    wordcloud2Output(outputId = "DAwordcloud")
                                            )
                                  )
                                )


                       )

                     )
                   )
            )
          )
)