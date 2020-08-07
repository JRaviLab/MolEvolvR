library(shiny)
library(shinyBS) #modals
library(shinyWidgets)


diagram_box <- function(title, image_name, button_name, description)
{
  div(class = "diagram-box",

      div(strong(title), class = "box-title"),
      div(description, class = "diagram-box-description"),


      #div(class = "image-border", style = "border-style: solid; border-width: 0px 1px 1px 1px; border-color: #AAAAAA;",
      div(class = "diagram-image",  style= paste0("background-image: url(", image_name, ".png);
          background-size: auto 80%; background-position: center; background-repeat: no-repeat;")),
      actionButton(button_name, NULL, class="diagram-button")



  )
}


basic_box <- function(title, image_name, description, height = '200px')
{
  div(class = "basic-box", #style = paste0("height:", height, ";") ,
      div(strong(title), class = "box-title")
      #div(description, class = "diagram-box-description", style = paste0("height:", height, ";")),
      #div(class = "diagram-image",  style= paste0("background-image: url(", image_name, ".png);
      #    background-size: auto 100%; background-position: center; background-repeat: no-repeat;"))
  )

}

sub_box <- function(title, image_name, description, button_name )
{
  div(class = "sub-box",
      h4(title),
      p(description),
      div(class = "diagram-image",  style= paste0("background-image: url(", image_name, ".png);
          background-size: 80% 80%; background-position: center; background-repeat: no-repeat;")),
      actionButton(button_name, NULL, class="sub-diagram-button")
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



construction_box <- function(text = "This page is under construction and some features
                             may not work")
{

  div( class = "construction-box", style = "color: black;",
       column(width = 12,
              icon("wrench", class= "inline-icon fa-2x")
              ,p(text, class = "inline-text")
       )
  )
}





