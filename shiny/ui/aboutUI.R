source("shiny/ui/aboutText.R")
library("here")

tabPanel(title = "About",
         value = "about",
        fluidRow(column(10,
                        h2('Citation:')
                        #put citations here


        )),
        fluidRow(column(10,
                        h2('Abstract:'),

                        p(abstract, style = "font-size:120%")

        )),
        fluidRow(column(10,
                        h2("Links:"),
                        # Link to the github repo using github icon
                        p("GitHub", style = "font-size:120%"),
                        #tags$br(),
                        tags$a(href = "https://github.com/jananiravi/psp-evolution", icon("github-square", "fa-3x", lib = "font-awesome")),
                        tags$br(),
                        p("Paper", style = "font-size:120%"),
                        tags$br(),
                        p("link to lab", style = "font-size:120%"),
                        icon("users","fa-3x",lib = "font-awesome"),

                        # Link to JRaviLab github page using cpathogeno icon
                        tags$a(href = "https://github.com/JRaviLab",
                                  tags$img(src = "/icons/jananiravi-worklogo-grey.png",
                                               alt = "group",
                                               width = "100",
                                               height = "100")
                        )




                        # tags$dl(
                        #   tags$dt("Links:"),
                        #   tags$dd("GitHub"),
                        #   tags$dd("Paper"),
                        #   tags$dd("link to lab"))


        )),
        fluidRow(column(10,
                        h2("Contacts:")


        ))

        #put link to paper, lab, github, other readings

)