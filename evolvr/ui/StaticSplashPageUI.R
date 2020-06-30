tabPanel("Home", value = "home", icon = icon("home", class = "fa-1x", lib = "font-awesome"),
         fluidPage(
           fluidRow(
           div(class = "SplashText", style = "text-align:center;",
           "The EvolvR webapp provides a streamlined approach for phylogenetic
                    and molecular evolution analysis."
           )
           ),
           fluidRow(
             div(class = "diagram-image",  style="background-image: url(ISMB-Poster-FlowChart-evolvr.png);
          background-size: auto 100%; background-position: center; background-repeat: no-repeat;")
             )
           )
         )