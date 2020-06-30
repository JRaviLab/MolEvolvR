
tabPanel("Home", value = "home", icon = icon("home", class = "fa-1x", lib = "font-awesome"),
fluidPage(
  #### Input Query ####
  fluidRow(
    column(width = 3, offset = 1,
           diagram_box("Query", image = "none", button_name = "dQueryBtn", description = "Description")
    ),
    column(width = 1, offset = 0,
           arrow_r()
    ),
    column(width = 3, offset = 0,
           diagram_box("Molecular Characterization of Target Protein", image = "none", button_name = "dMolTargetBtn", description = "Description")
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
           diagram_box("Molecular Characterization of Homologs", image = "none", button_name = "dMolHomoBtn", description = "Description") )
  ),
  fluidRow(
    column(width = 2, offset = 4, arrow_dl()),
    column(width = 1, offset = 0, arrow_d()),
    column(width = 2, offset = 0, arrow_dr())
  ),
  fluidRow(
    column(width = 3, offset = 2,
           diagram_box("Phyletic Spread", image = "none", button_name = "dPhyleticBtn", description = "Description")
    ),
    column(width = 3, offset = 0,
           diagram_box("Proximity Network", image = "none", button_name = "dNetworkBtn", description = "Description")
    ),
    column(width = 3, offset = 0,
           diagram_box("Structure Based Sequence Alignment", image = "none", button_name = "dAlignmentBtn", description = "Description")
    )
  )
))