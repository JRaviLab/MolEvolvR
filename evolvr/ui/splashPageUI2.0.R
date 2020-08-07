
tabPanel("Home", value = "home", icon = icon("home", class = "fa-1x", lib = "font-awesome"),
         fluidPage(
           loginUI("login"),
           #### Input Query ####
           uiOutput("splashUIComponent")


         ))