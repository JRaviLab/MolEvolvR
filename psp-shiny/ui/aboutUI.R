source("shiny/ui/aboutText.R")
library("here")

tabPanel(title = "About",
         value = "about",

         fluidRow(
                 column(width = 10, offset = 1,
                        htmlOutput("aboutApp")
                 )
         )


        # fluidRow(
        # ),
        # fluidRow(column(10,
        #                 htmlOutput("aboutAbstract")
        #
        # )),
        # fluidRow(column(10,
        #                 tags$a(href = "https://github.com/jananiravi/psp-evolution",icon("github-square", class = "zoom iMargin fa-3x", lib = "font-awesome")),
        #                 # Link to JRaviLab github page using cpathogeno icon
        #                 tags$a(href = "https://github.com/JRaviLab",
        #                        icon("users",class = "zoom iMargin fa-3x",lib = "font-awesome")
        #                 ),
        #
        #                 a(href = "mailto:janani@msu.edu",
        #                   icon(name = "envelope",class = "zoom iMargin fa-3x",
        #                        lib = "font-awesome")),
        #                 a(href = "https://twitter.com/JRaviLab",
        #                   icon(name = "twitter", class = "zoom iMargin fa-3x", lib = "font-awesome"))
        #
        #
        # ))
)