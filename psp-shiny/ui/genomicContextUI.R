### Genomic Context tab
tabPanel(title = "Genomic Context",
         value = "genomicContext",
         fluidRow(

           column(width = 4,
                      #dropdown to select protein for plots
                      selectInput(inputId =  "GClinSelec", label = "Protein",
                                  choices = c("All", queries)
                                  , selected = "PspA")
           ),
           column(width = 4,
                      #Slider input to determine cutoff value for totalcounts
                      sliderInput(inputId = "GC_Cutoff", label = "Total Count Cutoff:", min = 0, max = 100, value = 95)
           ),
           column(width = 4,
                  htmlOutput("GCLegend")
           )



           ,
           column(width = 12, offset = 0,
                  #mainpanel dictates what is displayed on screen depending on which tabset panel is selected
                  mainPanel(
                    width = 12,
                    tabsetPanel(
                      id= 'GCLin_data',
                      tabPanel("Heatmap", value = "Heatmap",
                               selectInput(inputId = "GC_lin_color", label = "Color Scheme", choices =
                                             c("default", "viridis", "inferno","magma", "plasma", "cividis"), selected = "default"
                               ),
                               plotOutput(outputId = "GCLinPlot", height = '600px', width = "1440px" )),
                      tabPanel("Summary Table", value = "LinTable",
                               p("Select a row to view lineage-level stats", style = "color:#242320;", class = "note-box"),
                               DT::dataTableOutput(outputId = "GCLinTable"),
                               column(downloadButton(outputId = "GCdownloadCounts", label = "Download"),radioButtons(inputId = "GCcountDownloadType", label = "Download Type:",
                                                                                                                     choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                      tabPanel("Upset Plot", value = "Upset",
                               plotOutput(outputId = "GCUpsetP", height = '600px')),
                      tabPanel("WordCloud", value = "Network_WC",
                               fluidRow(
                                 # box(width = 12,
                                 #     plotOutput(outputId = "GCNetwork")
                                 # ),
                                 box(width = 12,
                                     # plotOutput(outputId = "GCwordcloud")
                                     wordcloud2Output(outputId = "GCwordcloud")
                                 )
                               )

                      )
                    )
                  )
           )
         )
)