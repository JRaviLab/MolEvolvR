library(shiny)
library(shinyBS)
library(shinyWidgets)
library(tidyverse)

ui <- navbarPage(

  # tabPanel("Component 1", ... =  p("This is the c1")),
  # tabPanel("Component 2"),
  # navbarMenu("More",
  #            tabPanel("Sub-Component A"),
  #            tabPanel("Sub-Component B"))


tabsetPanel(
  tabPanel("AccNums/Fasta",
           column(
             width = 4, offset = 4,
             textAreaInput(inputId = "accNumBLASTTextInput", label = "Enter Protein Accession Number(s)",
                           width = "100%", height = "100%"),
           ),
           div(class = "note-box", style = "width: 80%; padding: 20px; margin: auto;",

               ## Don't do groups: I can make tool tips for each of the checkboxes then
               # prettyCheckboxGroup(inputId = 'testgroup',label = h4("Options"),
               #                     choices = c("DeltaBLAST",
               #                                 "Acc2FASTA",
               #                                 "BLASTClust"),
               #                     status = "primary",
               #                     shape = "curve"
               #                     ),
               # bsTooltip(id = "testgroup", title = "daflj"),
               popify(
               checkboxInput(inputId = "testnormCB",label = "Testing Normal"),title = "Normal CB"
              ,placement = "right", content = "This is more information about what this check box does"
               ),

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


           )

  )
)
)

server <- shinyServer(function(input, output, session) {

  # BLAST CLUST
  observeEvent(input$Acc_blastclustCB,
               {
                 # Selecting blast clust requires DeltaBlast
                 if(input$Acc_blastclustCB)
                 {
                    updatePrettyCheckbox(session, inputId = "Acc_deltaCB", label = "DeltaBLAST", value = T)
                    updatePrettyCheckbox(session, inputId = "Acc_acc2faCB", label = "Acc2FASTA", value = T)
                 }
                 # Unselecting it means that rps blast and ipr scan can't be run
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_rpsblastCB", label = "RPSBLAST", value = F)
                   updatePrettyCheckbox(session, inputId = "Acc_iprscanCB", label = "InterproScan", value = F)
                 }
               })

  # DeltaBLAST
  observeEvent(input$Acc_deltaCB,
               {
                 if(!input$Acc_deltaCB)
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = F)
                   updatePrettyCheckbox(session, inputId = "Acc_acc2faCB", label = "Acc2FASTA", value = F)
                 }

               })


  observeEvent(input$Acc_acc2faCB,
               {
                 if(input$Acc_acc2faCB)
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_deltaCB", label = "DeltaBLAST", value = T)
                 }
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = F)
                 }
               }
               )
  # Iprscan
  observeEvent(input$Acc_iprscanCB,
  {
    if(input$Acc_iprscanCB){
      updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = T)
    }
    else
    {
      updatePrettyCheckbox(session, inputId = "Acc_ipr2daCB", label = "InterproScan2DA", value = F)
    }
  })

  # rpsblast
  observeEvent(input$Acc_rpsblastCB,
               {
                 if(input$Acc_rpsblastCB){
                   updatePrettyCheckbox(session, inputId = "Acc_blastclustCB", label = "BLASTClust", value = T)
                 }
                 else
                 {
                   updatePrettyCheckbox(session, inputId = "Acc_rps2daCB", label = "RPSBLAST2DA", value = F)
                 }
               })
  # rps2da
  observeEvent(input$Acc_rps2daCB,
               {
                 if(input$Acc_rps2daCB){
                   updatePrettyCheckbox(session, inputId = "Acc_rpsblastCB", label = "RPSBLAST", value = T)
                 }
               })

  # ipr2da
  observeEvent(input$Acc_ipr2daCB,
               {
                 if(input$Acc_ipr2daCB){
                   updatePrettyCheckbox(session, inputId = "Acc_iprscanCB", label = "InterproScan", value = T)
                 }
               })

})
shinyApp(ui, server)