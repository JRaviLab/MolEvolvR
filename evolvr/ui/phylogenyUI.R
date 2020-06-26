tabPanel( title = "Phylogeny",
  value = "phylogeny",
        fluidPage(
          column(width = 2,
                 selectInput(inputId =  "PhyloSelect", label = "Protein",
                             choices = c("All"), selected = "All")),
          column(width = 12,
                 tabsetPanel(
                   id= "phylo",
                   tabPanel("Sunburst", value ="sunburst",
                            numericInput(inputId = "levels",label = "Number of Levels:" ,value = 2, min = 1),
                            sunburstOutput(outputId = "sunburst")
                   ),
                   tabPanel("Tree", value="Tree",
                            htmlOutput(outputId = "msaTree" )
                   ),

                   tabPanel("MSA", value="MSA",
                            htmlOutput(outputId="msaPlot")),
                   tabPanel("Paralog Table", value="Paralog",
                            p("Select a row for more information about those paralogs", style = "text-align:center", class = "note-box"),
                            DT::dataTableOutput(outputId = "ParalogTable")
                   )
                 )
          )
        )
)