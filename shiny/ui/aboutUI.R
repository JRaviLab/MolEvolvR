source("shiny/ui/aboutText.R")
library("here")

tabPanel(title = "About",
         value = "about",
         # loginUI("login"),
        fluidRow(
        ),
        fluidRow(column(10,
                        htmlOutput("aboutAbstract")

        )),
        fluidRow(column(10,
                        # h2("Links:"),
                        # # Link to the github repo using github icon
                        # p("GitHub", style = "font-size:120%"),
                        # #tags$br(),
                        tags$a(href = "https://github.com/jananiravi/psp-evolution",icon("github-square", class = "zoom iMargin fa-3x", lib = "font-awesome")),
                        # tags$br(),
                        # p("Paper", style = "font-size:120%"),
                        # tags$br(),
                        # p("link to lab", style = "font-size:120%"),


                        # Link to JRaviLab github page using cpathogeno icon
                        tags$a(href = "https://github.com/JRaviLab",
                               icon("users",class = "zoom iMargin fa-3x",lib = "font-awesome")
                        ),

                        a(href = "mailto:janani@msu.edu",
                          icon(name = "envelope",class = "zoom iMargin fa-3x",
                               lib = "font-awesome"))


        ))
        # fluidRow(column(10,
        #                 h2("Contacts:"),
        #                 a(href = "janani@msu.edu",
        #                   icon(name = "envelope",class = "zoom fa-3x",
        #                        lib = "font-awesome"))
        # ))

        #put link to paper, lab, github, other readings

)