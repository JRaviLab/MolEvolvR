source("shiny/ui/aboutText.R")
library("here")

tabPanel(title = "About",
         value = "about",
         # loginUI("login"),
        fluidRow(
        ),
        fluidRow(column(10,
                        h1("Phage-shock-protein (Psp) Envelope Stress Response:
                           Evolutionary History & Discovery of Novel Players"),

                        h4("Janani Ravi", tags$sup("1,2*"),", Vivek Anantharaman", tags$sup("3")
                        , ", Samuel Zorn Chen", tags$sup("1"), ", Pratik Datta", tags$sup("2"),
                        ", L Aravind", tags$sup("3*"),", Maria Laura Gennaro", tags$sup("2*"),"."),
tags$sup("1"),"Pathobiology and Diagnostic Investigation, Michigan State University, East Lansing, MI;",
tags$sup("2"), "Public Health Research Institute, Rutgers University, Newark, NJ; ",
tags$sup("3"),"National Center for Biotechnology Information, National Institutes of Health, Bethesda, MD.",
tags$br(),
"*Corresponding authors. janani@msu.edu; aravind@nih.gov; marila.gennaro@rutgers.edu ",

                        h2('Abstract:'),

                        p(abstract, style = "font-size:120%")

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