
tabPanel("Home", value = "home", icon = icon("home", class = "fa-1x", lib = "font-awesome"),
fluidPage(
  fluidRow(
    column(width = 4, offset = 4,
           diagram_box("Query", image = "none", button_name = "dQueryBtn", description = "Description")
    )),
  fluidRow(
    column(width = 4, offset = 4, arrow_d())
  ),
  fluidRow(
    column(width = 4, offset = 4,
           diagram_box("Molecular Characterization", image = "none", button_name = "dMolecularBtn", description = "Description") )
  ),
  fluidRow(
    column(width = 4, offset = 2, arrow_dl()),
    column(width = 4, offset = 0, arrow_dr())
  ),
  fluidRow(
    column(width = 4, offset = 0,
           diagram_box("Phyletic Spread", image = "none", button_name = "dPhyleticBtn", description = "Description")
    ),
    column(width = 4, offset = 4,
           diagram_box("Structure Based Sequence Alignment", image = "none", button_name = "dAlignmentBtn", description = "Description")
    )
  )

))