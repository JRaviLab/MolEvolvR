library(shiny)
library(shinyBS) #modals
library(shinyWidgets)

#****************************#
#### Upload UI Components ####
#****************************#

#### Full Data ####
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

#### BLAST Upload ####
blast_input_ui <- tagList(
  construction_box(),

  fluidRow(
    column(width = 4, offset = 0,
           fileInput(inputId = "blastFileUpload", label = "Choose Blast File")
    )
  )
)

#### IPRScan Upload ####
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

#### AccFASTA Upload ####
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

#*************************#
#### DATA UI Component ####
#*************************#

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


#**********************#
#### Result Summary ####
#**********************#

# Requires DomArch or DomArch.repeats column
rs_DomArch_component <-   tagList(
  fluidRow(
    column(width = 2, offset = 4,
           actionLink("rs2DomArch", label = h3("Domain Architecture")),
           p("Proximity Network")
    ),
    column(width = 6, offset = 0,
           visNetworkOutput(outputId = "rs_network")
    )
  )
)

# Requires GenContext column
rs_GenContext_component <- tagList(
  fluidRow(
    column(width = 6, offset = 0,
           plotOutput(outputId = "rs_gcHeatmap")
    ),
    column(width = 2, offset = 0,
           actionLink("rs2GenContext", label = h3("Genomic Context")),
           p("Lineage Heatmap")
    )
  )
)

# Requires Lineage Column
rs_sunburst_component <- tagList(
  fluidRow(
    column(width = 2, offset = 4,
           actionLink("rs2Phylogeny", label = h3("Phylogeny")),
           p("Sunburst")
    ),
    column(width = 6, offset = 0,
           sunburstOutput(outputId = "rs_sunburst")
    )
  )
)

# Requires FASTA File: use if fasta file present and no lin col?
rs_tree_component <- tagList(
  fluidRow(
    column(width = 2, offset = 4,
           actionLink("rs2PhylogenyTree", label = h3("Phylogeny")),
           p("Phylogenetic Tree")
    ),
    column(width = 6, offset = 0,
           plotOutput("rs_tree")
    )
  )
)






