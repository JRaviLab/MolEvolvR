tabItem("phylogeny",
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
          column(width = 12,
                 tabsetPanel(
                   id= "phylo",
                   tabPanel("Tree", value="Tree",
                            htmlOutput(outputId = "msaTree" ),
                            #sunburstOutput(outputId = "sunburst")
                            # fluidRow(
                            #   column(width = 2,
                            #          p("Note: the sunburst plot shown is for PspA and currently does not change with the dropdown", style = "font-size:160%")) ,
                            #   column(width = 5, sund2bOutput(outputId = "sund2b")),
                            #   column(width = 5, sunburstOutput(outputId = "sunburst"))
                            #   )
                   ),
                   tabPanel("Sunburst", value ="sunburst",
                           sunburstOutput(outputId = "sunburst")
                            ),
                   tabPanel("MSA", value="MSA",
                            htmlOutput(outputId="msaPlot")),
                   tabPanel("Paralog Table", value="Paralog",
                            DT::dataTableOutput(outputId = "ParalogTable",width = 1000)
                   )
                 )
          )
        )
)