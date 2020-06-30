library(shiny)
library(shinyBS) #modals
library(shinyWidgets)


diagram_box <- function(title, image_name, button_name, description)
{
  div(class = "diagram-box",

      div(title, class = "box-title"),
      div(description, class = "diagram-box-description"),


      #div(class = "image-border", style = "border-style: solid; border-width: 0px 1px 1px 1px; border-color: #AAAAAA;",
      div(class = "diagram-image",  style= paste0("background-image: url(", image_name, ".png);
          background-size: auto 130%; background-position: center; background-repeat: no-repeat;")),
      actionButton(button_name, NULL, class="diagram-button")



  )
}

diagram_elipse <- function(title, image)
{

  div(class = "diagram-elipse",
      div(title, class = "elipse")

      )
}


arrow_d <- function() {
  div(
    icon(name = "arrow-right", class = "arrow arrow-d fa-3x", lib = "font-awesome"),
  )
}

arrow_r <- function() {
  div(
    icon(name = "arrow-right", class = "arrow arrow-r fa-3x", lib = "font-awesome"),
  )
}

arrow_l <- function() {
  div(
    icon(name = "arrow-left", class = "arrow arrow-l fa-3x", lib = "font-awesome"),
  )
}

arrow_dl <- function(){
  icon(name = "arrow-right", class = "arrow arrow-dl fa-3x", lib = "font-awesome")
}

arrow_dr <- function(){
  icon(name = "arrow-right", class = "arrow arrow-dr fa-3x", lib = "font-awesome")
}


