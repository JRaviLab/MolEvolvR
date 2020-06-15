

tabPanel(title = "", icon = icon("home", class = "fa-2x", lib = "font-awesome"),
         value = "appInfo",
         loginUI("login"),
         fluidRow(
           column(width = 10, offset = 1,
                  htmlOutput("aboutApp")
           )
         )

)