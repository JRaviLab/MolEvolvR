

tabPanel(title = "", icon = icon("home", class = "fa-2x", lib = "font-awesome"),
         value = "appInfo",
         loginUI("login"),


         fluidRow(
         ),
         fluidRow(column(10,
                         htmlOutput("aboutAbstract")

         )),
         fluidRow(column(10,
                         tags$a(href = "https://github.com/jananiravi/psp-evolution",icon("github-square", class = "zoom iMargin fa-3x", lib = "font-awesome")),
                         # Link to JRaviLab github page using cpathogeno icon
                         tags$a(href = "https://github.com/JRaviLab",
                                icon("users",class = "zoom iMargin fa-3x",lib = "font-awesome")
                         ),

                         tags$a(href = "mailto:janani@msu.edu",
                           icon(name = "envelope",class = "zoom iMargin fa-3x",
                                lib = "font-awesome")),
                         tags$a(href = "https://twitter.com/JRaviLab",
                           icon(name = "twitter", class = "zoom iMargin fa-3x", lib = "font-awesome"))


         ))
         # fluidRow(
         #   column(width = 10, offset = 1,
         #          htmlOutput("aboutApp")
         #   )
         # )

)