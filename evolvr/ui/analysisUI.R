
tabPanel("Analysis",
         value = "analysis",
         fluidPage(
           column( width = 4, offset = 4,
                   fileInput(inputId = "accNumBLASTUpload", label = "Choose AccNum File")
           ),
           column(
             width = 4, offset = 4,
             textAreaInput(inputId = "accNumBLASTTextInput", label = "Enter Protein Accession Number(s)",
                           width = "100%", height = "100%"),
           ),
           column(
             width = 4, offset = 4,
             # textInput(inputId = "pinName", "Pin Name:"),
             actionButton("pinAccBlast", "BLAST"),
             bsAlert("successAccNumBlastAlert")
           )
         )

)