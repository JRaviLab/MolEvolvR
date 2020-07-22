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
          background-size: auto 100%; background-position: center; background-repeat: no-repeat;")),
      actionButton(button_name, NULL, class="diagram-button")



  )
}


basic_box <- function(title, image_name, description, height = '200px')
{
  div(class = "basic-box", #style = paste0("height:", height, ";") ,
      div(title, class = "box-title")
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
                          column( width = 4,
                                  selectInput("fastaRepresentativeType", label = "Choose Representative Observations:",
                                              choices = c("One per Lineage", "One per Species")
                                  )
                          )
                          ,
                          column(width = 12,
                                 actionButton("fullDF2Fasta", "Generate FASTA")
                          ),

                          column(width = 12, offset = 0,
                                 strong("FASTA Sequences"),
                                 div(class = "text-area-output",
                                     verbatimTextOutput(outputId = "DF2Fasta")
                                 )
                          ),

                          column(width = 5,
                                 selectInput("DFAlignmentTool", label = "Choose Tool for MSA:",
                                             choices = c("Muscle", "ClustalOmega", "ClustalW")
                                 )
                          ),
                          column(width = 12, offset = 0,
                                 actionButton("DF2msa", "Generate MSA")
                          ),
                          column(width = 12, offset = 0,
                                 strong("Aligned FASTA"),
                                 div(class = "text-area-output",
                                     verbatimTextOutput(outputId = "DF2AlignedFasta")
                                 )
                          )

                        )
               )
  )
)


# fasta_input_ui <- tagList(
#
#
#   column(width = 12, offset = 0,construction_box() ),
#
#   fluidRow(
#     column(width = 4, offset = 0,
#            fileInput(inputId = "fastaFileUpload", label = "Choose Fasta File")
#     )
#   ),
#   fluidRow(
#     column(width = 4, offset = 2,
#            strong("Or", style = "font-size: 20px; text-align:left;")
#     )
#   ),
#   fluidRow(
#     column( width = 12, offset = 0,
#             textAreaInput(inputId = "fastaTextInput", label = "Enter Fasta Sequence")
#     )
#   ),
#   fluidRow(
#     column(width = 3,
#            selectInput("FastaAlignmentTool", label = "Choose Tool for MSA:",
#                        choices = c("Muscle", "ClustalOmega", "ClustalW")
#            )
#     ),
#     column(width = 3, offset = 0,
#            actionButton("fasta2msaBtn", "Generate MSA")
#     )
#   ),
#   fluidRow(
#     column(width = 12, offset = 0,
#            div(class = "text-area-output",
#                verbatimTextOutput(outputId = "fasta2msa")
#            )
#     )
#   )
# )
#
# accNum_input_ui <- tagList(
#
#   construction_box(),
#
#   fluidRow(
#     column(
#       width = 12, offset = 0,
#       textAreaInput(inputId = "accNumTextInput", label = "Enter Protein Accession Number(s)"),
#     ),
#     column(
#       width = 12, offset = 0,
#       actionButton(inputId = "exampleAccNums", "Load Example"),
#     ),
#     column(
#       width = 12, offset = 0,
#       actionButton( "accnum2Fasta","Generate FASTA")
#     ),
#     column(
#       width = 12, offset = 0,
#       div(class = "text-area-output",
#           verbatimTextOutput(outputId = "generatedFasta")
#       )
#     ),
#     column(
#       width = 12, offset = 0,
#       selectInput( "accnumAlignmentTool", label = "Choose Tool for MSA:",
#                    choices = c("Muscle", "ClustalOmega", "ClustalW")
#       ),
#       actionButton(inputId = "accnum2msa", label = "Generate MSA"),
#       div(class = "text-area-output",
#           verbatimTextOutput(outputId = "accnumMSA")
#       )
#
#     )
#   )
# )

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
           fileInput(inputId = "interproFileUpload", label = "Choose Interproscan File"),
           radioButtons(inputId = "fileTypeIPRScan", label = "File Type",choices = "tsv")
    )
  ),
  fluidRow(
    column(width = 12, offset = 0,
           DT::dataTableOutput("IPRScanData")
    )
  )
)


acc_fasta_ui <- tagList(
  bsCollapse(id = "accCollapse",
             multiple = TRUE,
             bsCollapsePanel(title = "Accession Numbers", value = "accnum", style = "primary",
                             fluidRow(
                               column( width = 4, offset = 0,
                                       fileInput(inputId = "accNumUpload", label = "Choose AccNum File")
                               ),
                               column(
                                 width = 12, offset = 0,
                                 textAreaInput(inputId = "accNumTextInput", label = "Enter Protein Accession Number(s)",
                                               width = "100%", height = "100%"),
                               ),
                               column(
                                 width = 12, offset = 0,
                                 actionButton(inputId = "exampleAccNums", "Load Example"),
                               ),
                               column(
                                 width = 12, offset = 0,
                                 actionButton( "accnum2Fasta","Generate FASTA")
                               )
                             )

             ),
             bsCollapsePanel(title = "FASTA Sequences", value = "fasta", style = "primary",
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
                                       textAreaInput(inputId = "fastaTextInput", label = "Enter Fasta Sequence",
                                                     width = "100%", height = "100%")
                               )
                             ),
                             bsAlert(anchorId = "invalidPin"),
                             textInput(inputId = "pinName", "Pin Name:"),
                             actionButton("pinFasta", "Pin FASTA"),

                             fluidRow(
                               column(width = 3,
                                      selectInput("FastaAlignmentTool", label = "Choose Tool for MSA:",
                                                  choices = c("Muscle", "ClustalOmega", "ClustalW")
                                      )
                               ),
                               column(width = 3, offset = 0,
                                      actionButton("fasta2msaBtn", "Generate MSA")
                               )
                             )


             ),
             bsCollapsePanel("FASTA Sequence Alignment", value = "msa", style = "primary",
                             fluidRow(
                               column(width = 4, offset = 0,
                                      fileInput(inputId = "msaFileUpload", label = "Choose MSA File")
                               )
                             ),

                             fluidRow(
                               column(width = 12, offset = 0,

                                      textAreaInput(inputId = "msaText", label = "Paste  Aligned FASTA Sequence",
                                                    width = "100%", height = "100%")

                               )

                             )
             )


  )


)


#### DATA UI Component ####


noTableComponent <- tagList(
  column(width = 12,

         actionLink(inputId = "dataTable2Upload", label = "Upload a data table ", class = "inline-text"),
         p(" to view your data here.", class = "inline-text")

  )
)


tableComponent <- tagList( column(width = 2,
                                  #Dropdown to select protein for viewing
                                  selectInput(inputId =  "mainSelect", label = "Protein",
                                              choices = c("All"), selected = "All")

),

#Buttons to select which file type to download
column( width = 3, offset= 1,
        #Radiobuttons to select what to download data table as: tab separated or comma seperated
        radioButtons(inputId = "downloadType", label = "Download Type:",
                     choices= c("tsv", "csv"), selected = "tsv" ),
        #Output download button
        downloadButton(outputId = "downloadData", label = "Download")),
#Create mainpanel where dataTable is displayed
column(
  DT::dataTableOutput(outputId = "mainTable"), width = 12)
)

noFastaComponent <- tagList(
  p("Upload or generate a FASTA sequence to view the sequence here")
)

noMsaComponent <- tagList(
  p("Upload or generate an aligned fasta sequence for viewing here")
)


fastaComponent <- tagList(
           strong("FASTA"),
           div(class = "text-area-output",
               verbatimTextOutput( outputId = "fastaDataText")
           )
)

msaComponent <- tagList(
  strong("FASTA Sequence Alignment"),
  div(class = "text-area-output",
      verbatimTextOutput(outputId = "msaDataText")
  )
)



