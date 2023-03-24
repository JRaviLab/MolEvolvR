

tabPanel(title = "", icon = icon("home", class = "fa-2x", lib = "font-awesome"),
         value = "appInfo",
         loginUI("login"),
         fluidRow(
           column(width = 10, offset = 1,
                  h4("Background"),
                  "This web app was built to provide a visual and interactive supplement for the Psp Evolution paper.",



                  h4("Code"),
                  tags$span(
                    "Code and data used to generate this Shiny app are available on ",
                    a(href = "https://github.com/JRaviLab/the-approach",
                      class ="lightblue-link", "GitHub.")
                  )
           )
         )

)