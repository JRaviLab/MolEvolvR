library(shiny)
library(shinyBS) #modals
library(shinyWidgets)


diagram_box <- function(title, image, button_name, description)
{
  div(class = "diagram-box",

      div(title, class = "box-title"),
      div(description, class = "diagram-box-description"),

      actionButton(button_name, NULL, class="diagram-button")



      )



}

arrow_d <- function() {
  div(
  icon(name = "arrow-right", class = "arrow arrow-d fa-3x", lib = "font-awesome"),
  )
}

arrow_dl <- function(){
  icon(name = "arrow-right", class = "arrow arrow-dl fa-3x", lib = "font-awesome")
}

arrow_dr <- function(){
  icon(name = "arrow-right", class = "arrow arrow-dr fa-3x", lib = "font-awesome")
}

