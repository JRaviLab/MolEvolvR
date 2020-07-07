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


basic_box <- function(title, image_name, description, height = '200px')
{
  div(class = "basic-box", style = paste0("height:", height, ";") ,
      div(title, class = "box-title"),
      div(description, class = "diagram-box-description"),
      div(class = "diagram-image",  style= paste0("background-image: url(", image_name, ".png);
          background-size: auto 130%; background-position: center; background-repeat: no-repeat;"))
  )

}

sub_box <- function(title, image_name, description, button_name )
{
  div(class = "sub-box",
      h4(title),
      p(description),
      div(class = "diagram-image",  style= paste0("background-image: url(", image_name, ".png);
          background-size: auto 130%; background-position: center; background-repeat: no-repeat;")),
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



#### UI components for the upload page ####

full_data_ui <- tagList(

  tabsetPanel( id = "fullDataTabs",
               tabPanel(title = "Full Data Table", value = "fullDataTab",
                        fluidRow(
                          column(width = 12, offset = 0,
                                 fileInput(inputId = "fileUpload",
                                           label = "Choose File"
                                 ),
                                 radioButtons(inputId = "fileType", label = "File Type:",
                                              choices = c("tsv", "csv")),
                                 actionButton(inputId = "loadExample", "Example Data"),



                                 textInput(inputId = "queryInput", label = "Queries:", value = ""),

                                 tags$div( class = "textOutput-box",
                                           textOutput("queryVect")
                                 ),

                                 tableOutput("uploadedContent")
                          )
                        )
               ),
               tabPanel(title = "FASTA", value = "fastaTab",
                        construction_box(),
                        fluidRow(
                          selectInput("fastaRepresentativeType", label = "Choose Representative Observations:",
                                      choices = c("One per Lineage", "One per Species")
                          ),
                          actionButton("fullDF2Fasta", "Generate FASTA"),

                          column(width = 12, offset = 0,
                                 div(class = "text-area-output",
                                     verbatimTextOutput(outputId = "DF2Fasta")
                                 )
                          )

                        )
               )
  )
)


fasta_input_ui <- tagList(


  column(width = 12, offset = 0,construction_box() ),

  fluidRow(
    column(width = 4, offset = 0,
           fileInput(inputId = "fastaFileUpload", label = "Choose Fasta File")
    )
  ),
  fluidRow(
    column(width = 4, offset = 2,
           strong("Or", style = "font-size: 20px; text-align:left;")
    )
  ),
  fluidRow(
    column( width = 12, offset = 0,
            textAreaInput(inputId = "fastaTextInput", label = "Enter Fasta Sequence")
    )
  )
)

accNum_input_ui <- tagList(

  construction_box(),

  fluidRow(
    column(
      width = 12, offset = 0,
      textAreaInput(inputId = "accNumTextInput", label = "Enter Protein Accession Number(s)"),
    ),
    column(
      width = 12, offset = 0,
      actionButton(inputId = "exampleAccNums", "Load Example"),
    ),
    column(
      width = 12, offset = 0,
      div(class = "text-area-output",
          verbatimTextOutput(outputId = "generatedFasta")
      )
    )
  )
)

blast_input_ui <- tagList(
  construction_box(),

  fluidRow(
    column(width = 4, offset = 0,
           fileInput(inputId = "blastFileUpload", label = "Choose Blast File")
    )
  )
)

interpro_input_ui <- tagList(
  construction_box(),

  fluidRow(
    column(width = 4, offset = 0,
           fileInput(inputId = "interproFileUpload", label = "Choose Interproscan File")
    )
  )
)





