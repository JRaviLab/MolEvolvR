tabPanel( title = "Phylogeny",
  value = "phylogeny",
        fluidPage(
          column(width = 2,
                 selectInput(inputId =  "alignSelec", label = "Protein",
                             choices = c(
                               "LiaI-LiaF-TM.allfa50",
                               "LiaI-LiaF-TM_LiaFN.2",
                               "LiaI-LiaF-TM_LiaI.1",
                               "LiaI-LiaF-TM_PspC.3",
                               "PspA Only",
                               "PspA Snf7 Gismo",
                               "PspA Snf7",
                               "PspB Gismo",
                               "PspC Gismo",
                               "Snf7 Only",
                               "Toast-rack DUF2154-LiaF",
                               "Toast-rack DUF2807",
                               "Toast-rack DUF4097",
                               "Toast-rack PspC-Cterm"
                               )
                             , selected = "PspA Snf7")),
          column(width = 4, offset = 6,
                 # textOutput("DALegend", inline = T)
                 htmlOutput("PhyloLegend", inline = T)
          ),
          column(width = 12,
                 tabsetPanel(
                   id= "phylo",
                   tabPanel("Sunburst", value ="sunburst",
                            p("Hover over the colored segments for expanded lineage info.", style = "text-align:center", class = "note-box"),
                            numericInput(inputId = "levels",label = "Number of Levels:" ,value = 2, min = 1, max = 5),
                            sunburstOutput(outputId = "sunburst")
                   ),
                   tabPanel("Tree", value="Tree",
                            htmlOutput(outputId = "msaTree" )
                   ),

                   tabPanel("MSA", value="MSA",
                            htmlOutput(outputId="msaPlot")),
                   tabPanel("Paralog Table", value="Paralog",
                            p("Select a row for more information about the paralogs", style = "text-align:center", class = "note-box"),
                            DT::dataTableOutput(outputId = "ParalogTable")
                   )
                 )
          )
        )
)