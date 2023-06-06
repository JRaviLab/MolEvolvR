library(shiny)
library(shinyBS)
library(shinyWidgets) # Used to set the background image
source("components.R")

ui <-
  tagList(
    includeCSS("components.css"),
    shinyUI(fluidPage(

      # setBackgroundImage(src = "the_approach_flowchart.png"),
      fluidRow(
        column(width = 3, offset = 1,
               # diagram_box("Query", image = "None", button_name = "dQueryBtn", description = "Description"),
               basic_box("Query Protein", image = "none", description = "Description", height = '70px'),
               sub_box("Protein Accession Number", image_name = "none", description = "", button_name = "dAccNumBtn"),
               sub_box("Fasta", image_name = "none", description = "", button_name = "dFastaBtn")
        ),
        column(width = 1, offset = 0,
               arrow_r()
        ),
        column(width = 3, offset = 0,
               diagram_box("Molecular Characterization of Target Protein", image = "mol_characterization", button_name = "dMolTargetBtn", description = "Description")
        ),
        column(width = 1, offset = 0,
               arrow_r()
        ),
        column(width = 3, offset = 0,
               diagram_box("Identify Homologs", image = "none", button_name = "dHomoBtn", description = "Description")
        )
      ),
      fluidRow(
        column(width = 3, offset = 9, arrow_d())
      ),
      fluidRow(
        column(width = 3, offset = 5,
               diagram_elipse("Analysis", image = "none") ),
        column(width = 1, offset = 0, arrow_l()),
        column(width = 3, offset = 0,
               diagram_box("Molecular Characterization of Homologs", image = "mol_characterization", button_name = "dMolHomoBtn", description = "Description") )
      ),
      fluidRow(
        column(width = 2, offset = 4, arrow_dl()),
        column(width = 1, offset = 0, arrow_d()),
        column(width = 2, offset = 0, arrow_dr())
      ),
      fluidRow(

        column(width = 3, offset = 2,
               diagram_box("Phyletic Spread", image = "heatmap", button_name = "dPhyleticBtn", description = "Description")
        ),
        column(width = 3, offset = 0,
               diagram_box("Proximity Network", image = "domain_network", button_name = "dNetworkBtn", description = "Description")
        ),
        column(width = 3, offset = 0,
               diagram_box("Structure Based Sequence Alignment", image = "phylo_tree", button_name = "dAlignmentBtn", description = "Description")
        )
      ),

      ##############################################################
      fluidRow(
        column(width = 4, offset = 4,
               diagram_box("Query", image = "none", button_name = "dQueryBtn", description = "Description")
        )),
      fluidRow(
        column(width = 4, offset = 4, arrow_d())
      ),

      #### Molecular Characterization ####
      fluidRow(
        column(width = 4, offset = 4,
               diagram_box("Molecular Characterization of Target", image = "none", button_name = "dMolTargetBtn", description = "Description") )
      ),
      fluidRow(
        column(width = 4, offset = 4, arrow_d())
      ),
      fluidRow(
        column(width = 4, offset = 4,
               diagram_box("Identify Homologs", image = "none", button_name = "dHomologsBtn", description = "Description") )
      ),
      fluidRow(
        column(width = 4, offset = 4, arrow_d())
      ),
      fluidRow(
        column(width = 4, offset = 4,
               diagram_box("Molecular Characterization of Homologs", image = "none", button_name = "dMolHomoBtn", description = "Description") )
      ),
      fluidRow(
        column(width = 2, offset = 3, arrow_dl()),
        column(width = 2, offset = 0, arrow_d() ),
        column(width = 2, offset = 0, arrow_dr())
      ),

      #### Results ####
      fluidRow(
        column(width = 4, offset = 0,
               diagram_box("Phyletic Spread", image = "none", button_name = "dPhyleticBtn", description = "Description")
        ),
        column(width = 4, offset = 0,
               diagram_box("Proximity Network", image = "none", button_name = "dNetworkBtn", description = "Description")
        ),
        column(width = 4, offset = 0,
               diagram_box("Structure Based Sequence Alignment", image = "none", button_name = "dAlignmentBtn", description = "Description")
        )

      )



      ,
      tags$head(

        #     tags$style(HTML('#query{background-color:orange;}
        #/*
        #                     .transparentButton{
        #                       background-color: Transparent;
        #
        #           }
        # */
        #
        # .flowchart-image {
        #   background-image: url("the_approach_flowchart.png");
        #   background-color: #cccccc;
        #   height: 100%;
        #   width: 100%;
        #   background-position: center;
        #   background-repeat: no-repeat;
        #   /*background-size: cover;*/
        #   background-size: 100% 1000px;
        #   position: absolute; /*relative;*/
        #   overflow-y: scroll;
        # }
        #
        #                     '))
      )
    )))
server <- shinyServer(function(input, output) {

})
shinyApp(ui, server)


