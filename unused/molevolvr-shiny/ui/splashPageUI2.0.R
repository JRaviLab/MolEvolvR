loginUI2 <- function(id, title = "Please log in", user_title = "User Name", pass_title = "Password",
                   login_title = "Log in", error_message = "Invalid username or password!") {
  ns <- shiny::NS(id)

  shiny::div(id = ns("panel"), style = "width: 500px; max-width: 100%; margin: 0 auto; padding: 20px;",
             shiny::wellPanel(
               shiny::tags$h2(title, class = "text-center", style = "padding-top: 0;"),

               shiny::textInput(ns("user_name"), shiny::tagList(shiny::icon("user"), user_title)),

               shiny::passwordInput(ns("password"), shiny::tagList(shiny::icon("unlock-alt"), pass_title)),

               shiny::div(
                 style = "text-align: center;",
                 shiny::actionButton(ns("button"), login_title, class = "btn-primary", style = "color: white;")
               ),

               shiny::div(
                # style = "display: inline;
   # vertical-align: top;", #float: left;",
                  shiny::tags$p("Contact ", style = "display:inline-block;"),
                  shiny::tags$a(href = "mailto:janani@msu.edu", "janani@msu.edu"),
                  shiny::tags$p( "for the login credentials.", style = "display:inline-block;")
               ),

               shinyjs::hidden(
                 shiny::div(id = ns("error"),
                            shiny::tags$p(error_message,
                                          style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
               )
             )
  )
}


tabPanel("Home", value = "home", icon = icon("home", class = "fa-1x", lib = "font-awesome"),
         fluidPage(
           loginUI2("login"),
           #### Input Query ####
           uiOutput("splashUIComponent")


         ))