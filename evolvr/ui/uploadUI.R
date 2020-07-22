


tabPanel(
  title = "Upload",
  icon = icon("upload", lib = "font-awesome"),
  value = "upload",
  fluidPage(
    # construction_box(),



    div( class = "upload-header",
    fluidRow(

           column(width = 3, offset = 0,
                  selectInput(inputId = "inputType", label = "Input Type:",
                              choices = c(
                                # "Protein Accession Numbers",
                                # "Fasta Sequence(s)",
                                "Blast Results",
                                "Interproscan Results",
                                "Full Data",
                                "AccNum/FASTA"
                              ),
                              selected = "Full Data"
                  )
           ),
           column(width = 3, offset = 6,
                  actionButton(inputId = "upload2RS", label = "View Results")
                  )
      )

    ),
    # Dynamically generate ui components depending on the selected input type
    uiOutput(outputId = "uploadComponents")

  )
)