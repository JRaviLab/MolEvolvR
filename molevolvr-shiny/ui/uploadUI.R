


tabPanel(
  title = "Upload",
  icon = icon("upload", lib = "font-awesome"),
  value = "upload",
  fluidPage(
    # construction_box(),

    actionButton(inputId = "AddLinsBttn", label = "Add Lineage"),

    div( class = "upload-header",
         fluidRow(
           # tabset panel instead:
           tabsetPanel(id = "uploadTabs",
                       tabPanel(title = "AccNum/FASTA",value = "AccFastaTab",

                                tabsetPanel(id = "AccType", type = "pills",
                                            tabPanel("Homologous", value = "AccHomo",

                                                     tagList(
                                                       bsCollapse(id = "accCollapse",
                                                                  multiple = TRUE,
                                                                  bsCollapsePanel(title = "Homolog Accession Numbers", value = "accnum", style = "primary",
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
                                                                                    ),
                                                                                    column(
                                                                                      width = 12, offset = 0,
                                                                                      actionButton(inputId = "exampleFASTA", "Load Example"),
                                                                                    )
                                                                                  ),
                                                                                  # bsAlert(anchorId = "invalidPin"),
                                                                                  # textInput(inputId = "pinName", "Pin Name:"),
                                                                                  # actionButton("pinFasta", "Pin FASTA"),

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

                                                                                    ),
                                                                                    column(
                                                                                      width = 12, offset = 0,
                                                                                      actionButton(inputId = "exampleMSA", "Load Example"),
                                                                                    )

                                                                                  )
                                                                  )
                                                       )
                                                     )),
                                            tabPanel(title = "NonHomologous", value = "AccNoHomo",

                                                     fluidPage(
                                                       column( width = 4, offset = 4,
                                                               fileInput(inputId = "accNumBLASTUpload", label = "Choose AccNum File")
                                                       ),
                                                       column(
                                                         width = 4, offset = 4,
                                                         textAreaInput(inputId = "accNumBLASTTextInput", label = "Enter Protein Accession Number(s)",
                                                                       width = "100%", height = "100%"),
                                                       ),

                                                       div(style = "width: 80%; padding: 20px; margin: auto;",
                                                           # popify(
                                                           #   checkboxInput(inputId = "testnormCB",label = "Testing Normal"),title = "Normal CB"
                                                           #   ,placement = "right", content = "This is more information about what this check box does"
                                                           # ),

                                                           # Left Column
                                                           div(style = "width:50%; float: left;",
                                                               h4("Homology Search"),
                                                               prettyCheckbox(inputId = "Acc_deltaCB", label = "DeltaBLAST",
                                                                              inline = F,shape = "curve", # animation = "pulse",
                                                                              status = "primary"),
                                                               bsTooltip(id = "Acc_deltaCB", title = "DeltaBLAST will be run and ...",
                                                                         placement = "right"),

                                                               # doesn't work with shiny widgets
                                                               popify(
                                                                 prettyCheckbox(inputId = "Acc_acc2faCB", label = "Acc2FASTA",
                                                                                inline = F,shape = "curve",
                                                                                status = "primary"), title = "Acc2FASTA", content = "Selecting this option
                             will convert the accession numbers produced by the DeltaBLAST run into
                             FASTA files.",trigger = "click"
                                                               ),

                                                               prettyCheckbox(inputId = "Acc_blastclustCB", label = "BLASTCLUST",
                                                                              inline = F,shape = "curve",
                                                                              status = "primary")
                                                           ),
                                                           # Right Column
                                                           div(style = "width:50%; float: left;",
                                                               h4("Domain Architectures"),
                                                               prettyCheckbox(inputId = "Acc_iprscanCB", label = "InterproScan",
                                                                              inline = F,shape = "curve",
                                                                              status = "primary"),
                                                               prettyCheckbox(inputId = "Acc_rpsblastCB", label = "RPSBLAST",
                                                                              inline = F,shape = "curve",
                                                                              status = "primary"),
                                                               prettyCheckbox(inputId = "Acc_ipr2daCB", label = "InterproScan2DA",
                                                                              inline = F,shape = "curve",
                                                                              status = "primary"),
                                                               prettyCheckbox(inputId = "Acc_rps2daCB", label = "RPSBLAST2DA",
                                                                              inline = F,shape = "curve",
                                                                              status = "primary")
                                                           )
                                                       ),
                                                       column(width = 4, offset = 4,
                                                              textInput("acc2blastEmail","Email")
                                                              ),

                                                       column(
                                                         width = 4, offset = 4,
                                                         # textInput(inputId = "pinName", "Pin Name:"),
                                                         actionButton("pinAccBlast", "BLAST"),
                                                         bsAlert("successAccNumBlastAlert")
                                                       )
                                                     )

                                            )
                                )


                       ),
                       tabPanel(title = "BLAST Output", value = "BlastOutTab",
                                fluidRow(
                                  column(width = 4, offset = 0,
                                         fileInput(inputId = "blastFileUpload", label = "Choose Blast File")
                                  )
                                ),
                                column(width = 12, offset = 0,
                                       DT::dataTableOutput("BLASTData")
                                       )
                       ),
                       tabPanel(title = "Interproscan Output", value = "IprOutTab",
                                fluidRow(
                                  column(width = 4, offset = 0,
                                         fileInput(inputId = "interproFileUpload", label = "Choose Interproscan File"),
                                         actionButton(inputId = "loadIPRExample", "Example Data"),
                                         radioButtons(inputId = "fileTypeIPRScan", label = "File Type",choices = "tsv")
                                  )
                                ),
                                fluidRow(
                                  column(width = 12, offset = 0,
                                         DT::dataTableOutput("IPRScanData")
                                  )
                                )

                       ),
                       tabPanel(title = "Full Data", value = "FullDataTab",
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
                                                      # construction_box(),
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
                       ),
                       tabPanel(title = "Analysis", value = "AnalysisTab",
                                div(style= "text-align:center;",#  margin-left: auto;margin-right: auto;",

                                    column(
                                      p("Input the code for your data:"),
                                      width = 4, offset = 4, textInput(inputId = "analysisCode", label = ""),
                                      shinyWidgets::actionBttn(inputId = "FetchAnalysisBtn", label = "Fetch", style = "pill"),

                                      bsAlert("FetchAlert")
                                    )

                                ),

                                column(width = 12,
                                       tabsetPanel( id = "FetchedDataTabs",type = "pills",
                                                    tabPanel("BLAST Results", value = "FetchedBlastTab",
                                                             DT::dataTableOutput(outputId = "fetchedBlastOut")
                                                    ),
                                                    tabPanel("InterproScan Results", value = "FetchedIprTab",
                                                             DT::dataTableOutput(outputId = "fetchedIprOut")
                                                    )
                                       )
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
)