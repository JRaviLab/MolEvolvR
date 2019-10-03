#For Temp App, Commented out phylogeny, Removed DUFs

library(shiny)
library(tidyverse)
library(DT)
library(rmarkdown)
library(shinydashboard)
library(shinyjqui)
library(wordcloud)
library(shinyauthr)
library(svgPanZoom)
conflicted::conflict_prefer("intersect", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("strsplit", "base")
conflicted::conflict_prefer("count", "dplyr")
setwd("..")
source("shiny/PSP_Web_Data.R")
source("R/plotting.R")
#source("R/cleanup.R")
#source("R/reverse_operons.R")
source("shiny/shinyfunctions.R")

#########
##Users##
#########
user_base <- data.frame(
  user = c("pspevolution"),
  password = c("cpathogeno2019"),
  permissions = c("admin"),
  name = c("User One"),
  stringsAsFactors = FALSE
)



#########
#Sidebar#
#########
sidebar<- dashboardSidebar(
  sidebarMenu(id = "mainTabs",
              menuItem("Data Table", tabName = "datatable"),
              menuItem("Lineage Plots", tabName = "lineagePlots")
              ,
              menuItem("Phylogeny", tabName = "phylogeny")
              ,menuItem("Usage", tabName="usage")
              ,menuItem("About",tabName="about")
  )
)
######
#Body#
######
body <- dashboardBody(
  # must turn shinyjs on
  shinyjs::useShinyjs(),
  # add logout button UI
  div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
  # add login panel UI function
  shinyauthr::loginUI(id = "login"),

  uiOutput("user_table"),
  uiOutput("testUI"),
  tabItems(
    #Datatable tab contains all protein datatables
    tabItem("datatable",
            fluidPage(
              column(width = 2,
                     #Dropdown to select protein for viewing
                     selectInput(inputId =  "proSelec", label = "Protein",
                                 choices = c( "All","DUF1700", "DUF1707",
                                              "PspA", "PspB", "PspC", "PspM", "PspN","LiaI-LiaF-TM","LiaG","Toast-rack" )
                                 , selected = "All")
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
                DT::dataTableOutput(outputId = "proTable"), width = 12)))
    ,
    #lineagePlots contains a heatmap, datatable, and upset plot
    tabItem("lineagePlots",
            fluidPage(
              sidebarLayout(
                sidebarPanel(
                  #dropdown to select protein
                  selectInput(inputId =  "linSelec", label = "Protein",
                              choices = c( "PspA", "PspB", "PspC","PspN", "LiaF","Toast-rack")
                              , selected = "PspA"),
                  #Radiobuttons to selext domain architecture and genomic context
                  radioButtons(inputId = "DA_GC", label = "Lineage by:"
                               , choices= c("Domain Architecture", "Genomic Context")
                               , selected = "Domain Architecture"),
                  #Slider input to determine cutoff value for totalcounts
                  sliderInput(inputId = "cutoff", label = "Total Count Cutoff:", min = 0, max = 500, value = 30)
                  #sliderInput(inputId = "nr", label="Rows",min=0,max=10,value=0)
                ),
                #mainpanel dictates what is displayed on screen depending on which tabset panel is selected
                mainPanel(
                  tabsetPanel(
                    id= 'lin_data',
                    tabPanel("Heatmap",plotOutput(outputId = "LinPlot", height = '500px' )),
                    tabPanel("Table", DT::dataTableOutput(outputId = "LinTable"),
                             column(downloadButton(outputId = "downloadCounts", label = "Download"),radioButtons(inputId = "countDownloadType", label = "Download Type:",
                                                                                                                 choices= c("tsv", "csv"), selected = "tsv" ),width = 10)),
                    tabPanel("Upset Plot",plotOutput(outputId = "upsetP")),
                    tabPanel("Wordcloud", plotOutput(outputId = "wordcloud"))
                  )
                )


              )
            )
    ),
    tabItem("phylogeny",
            fluidPage(
              column(width = 2,
                     selectInput(inputId =  "alignSelec", label = "Protein",
                                 choices = c( "PspA")
                                 , selected = "PspA")),
              column(width = 12,
                     tabsetPanel(
                       id= "phylo",
                       tabPanel("Tree", value="Tree",
                                #,radioButtons(inputId = "plottype", label = "Plot Type",
                                #                   choices = c("Tree1", "Tree2", "Tree3"),selected = "Tree1"),

                                tags$head(tags$script(src = "http://www.elevateweb.co.uk/wp-content/themes/radial/jquery.elevatezoom.min.js")),
                                actionButton("myBtn", "Press Me for zoom!"),
                                p("If this button does not work, check if your browser is blocking this script from running"),
                                htmlOutput(outputId = "msaTree" ),
                                singleton(
                                  tags$head(tags$script('Shiny.addCustomMessageHandler("testmessage",
  function(message) {
    var image = $("#msaTree img");
    var zoomConfig = {scrollZoom : true};
    if(message.value == "ZoomOn"){
        //$("#msaTree img").elevateZoom({scrollZoom : true});
        image.elevateZoom(zoomConfig);
    }
    else{
       $.removeData(image, "elevateZoom");//remove zoom instance from image;

       $(".zoomContainer").remove();// remove zoom container from DOM;
    }
  }
);'))
                                )
                       ),
                       tabPanel("MSA", value="MSA",
                                htmlOutput(outputId="msaPlot")),
                       tabPanel("Paralog Table", value="Paralog",
                                DT::dataTableOutput(outputId = "ParalogTable",width = 1000)
                       )
                     )
              )
            )
    ),
    tabItem("about",
            fluidRow(column(10,
                            h2('Citation:')
                            #put citations here


                            )),
            fluidRow(column(10,
                            h2('Abstract:'),

                            p("The phage shock protein (Psp) stress-response system protects bacteria from envelope stress and
                            stabilizes the cell membrane. Despite the prevalence of the key effector, PspA, and the functional
                            Psp system, the various genomic contexts of Psp proteins, as well as their evolution across the kingdoms
                            of life, have not yet been characterized. Recent work from our group suggests that the psp systems have
                            evolved independently in distinct Gram-positive and Gram-negative bacterial clades to effect similar stress
                            response functions. We developed a computational pipeline for comparative genomics and protein
                            sequence-structure-function analyses to identify sequence homologs, phyletic patterns, domain architectures,
                            gene neighborhoods, sequence conservation and evolution of the candidates across the tree of life. Using
                            contextual information from conserved gene neighborhoods and their domain architectures, we delineated the
                            phyletic patterns of all the Psp members. Next, we systematically identified all possible 'flavors' and
                            genomic neighborhoods of the Psp systems. Finally, we have traced their evolution leading us to several
                            interesting observations as to their occurrence and co-migration, suggesting their function and role in
                            stress-response systems that are often lineage-specific. Conservation of the Psp systems across bacterial
                            phyla emphasizes the established importance of this stress response system in prokaryotes, while the
                            modularity in various lineages is indicative of adaptation to bacteria-specific cell-envelope structures,
                            lifestyles, and adaptation strategies.", style = "font-size:120%")

                            )),
            fluidRow(column(10,
                            h2("Links:"),
                            p("GitHub", style = "font-size:120%"),
                            tags$br(),
                            p("Paper", style = "font-size:120%"),
                            tags$br(),
                            p("link to lab", style = "font-size:120%")



                            # tags$dl(
                            #   tags$dt("Links:"),
                            #   tags$dd("GitHub"),
                            #   tags$dd("Paper"),
                            #   tags$dd("link to lab"))


                            )),
            fluidRow(column(10,
                            h2("Contacts:")


            ))

                            #put link to paper, lab, github, other readings

            ),
    tabItem("usage",
            fluidRow(
              column(12,
                     tabsetPanel(
                       id="manual_panel",
                       tabPanel("Data Table"),
                       tabPanel("Lineage Plots"),
                       tabPanel("Phylogeny")


                     )
                     )

            )
            #Use tabs to have different aspects of the page:
            #datatable(main)
            #Lineage tab, with descriptions for the various tools
            #Phylogeny, again with descriptions for the various tools


            )
  ),
  HTML('<div data-iframe-height></div>')
)
####
#UI#
####
ui <- dashboardPage(
  skin = "green",
                    #App title
                    dashboardHeader(title = "PSP Data"),
                    #Create Sidebar with inputs and download button output
                    sidebar,
                    body
)

########
#Server#
########
server <- function(input, output,session){
  # call the logout module with reactive trigger to hide/show
  logout_init <- callModule(shinyauthr::logout,
                            id = "logout",
                            active = reactive(credentials()$user_auth))

  # call login module supplying data frame, user and password cols
  # and reactive trigger
  credentials <- callModule(shinyauthr::login,
                            id = "login",
                            data = user_base,
                            user_col = user,
                            pwd_col = password,
                            log_out = reactive(logout_init()))

  # pulls out the user information returned from login module
  user_data <- reactive({credentials()$info})

  observe({
    if(credentials()$user_auth) {
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    } else {
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
  })


  #Reactive expression used to determine which data table should be displayed
  #based on selected input
  pspTable<- reactive({
    req(credentials()$user_auth)
    switch(input$proSelec,
           "All" = all,
           "DUF1700" = all %>% filter(Query=="^DUF1700"),
           "DUF1707" = all %>% filter(Query=="^DUF1707-SHOCT"),
           "PspA" = all%>% filter(Query=="pspa"),
           "PspB" = all%>% filter(Query=="pspb"),
           "PspC" = all%>% filter(Query=="pspc"),
           "PspM" = pspm_data,
           "PspN" = all%>% filter(Query=="pspn"),
           "LiaI-LiaF-TM" = all%>% filter(Query=="LiaI-LiaF-TM"),
           "Toast-rack" = all%>%filter(Query=="Toast-rack"),
           "LiaG" = liag_data)
  })
  #Render the Data table for selected protein
  output$proTable <- DT::renderDT({
    req(credentials()$user_auth)
    paged_table(pspTable())}, extensions = c('FixedColumns',"FixedHeader"),
                                  options = list(pageLength = 10,
                                                 #The below line seems to disable other pages and the search bar
                                                 #dom = 't',
                                                 scrollX = TRUE,
                                                 paging=TRUE,
                                                 fixedHeader=TRUE,
                                                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))

  observe({
    if(input$DA_GC== "Domain Architecture"){
    switch(input$linSelec,
        "PspA"= {hs <-heatmap_slider(pspa_DA_lin,"da")},
        "PspB"= {hs <-heatmap_slider(pspb_DA_lin,"da")},
        "PspC"= {hs <-heatmap_slider(pspc_DA_lin,"da")},
        "PspN"= {hs <-heatmap_slider(pspn_DA_lin,"da")},
        "LiaF"= {hs <-heatmap_slider(liai_liaf_DA_lin,"da")},
        #"LiaG"= {hs <-heatmap_slider(liag_DA_lin)},
        "Toast-rack"= {hs <-heatmap_slider(toast_rack_DA_lin,"da")}
    )}
    else{
      switch(input$linSelec,
             "PspA"= {hs <-heatmap_slider(pspa_GC_lin,"gc")},
             "PspB"= {hs <-heatmap_slider(pspb_GC_lin,"gc")},
             "PspC"= {hs <-heatmap_slider(pspc_GC_lin,"gc")},
             "PspN"= {hs <-heatmap_slider(pspn_GC_lin,"gc")},
             "LiaF"= {hs <-heatmap_slider(liai_liaf_GC_lin,"gc")},
             "Toast-rack"= {hs <-heatmap_slider(toast_rack_GC_lin,"gc")}
      )}
    updateSliderInput(session,inputId = "cutoff",min=0, max=hs[,"max"], value=hs[,"cutoff_init"])
  })

#   prot_da_lin <- reactive({switch(input$linSelec,
#                         "PspA"= pspa_DA_Lin,
#                         "PspB"= pspb_DA_Lin,
#                         "PspC"= pspc_DA_Lin,
#                         "LiaG"= liag_DA_Lin,
#                         "LiaF"= liaf_DA_Lin,
#                         "LiaI"= liai_DA_Lin
#   )})
#     prot<- reactive({
#       prot_da_lin%>% group_by(DomArch.norep)%>%
#         summarise(totalcount = sum(count)) %>%
#         filter(totalcount >1) %>% arrange(desc(totalcount))
#     })
#
#     prot_filt <- reactive({
#       filter(prot,totalcount >= input$cutoff)
#     })
# current_row <- reactive({
#   length(prot_filt$totalcount)
#   }
# )
# current_cutoff <- reactive({
#   prot[current_row,"totalcount"]
# })
#   observeEvent(input$nr,{
#     updateSliderInput(session,inputId = "cutoff",value = current_cutoff)
#   })
#
#   observeEvent(input$cutoff,  {
#     updateSliderInput(session = session, inputId = "nr", value = current_row)
#   })


#   observe({
#     current_row <- as.numeric(input$nr)
#     prot_da_lin <- switch(input$linSelec,
#            "PspA"= pspa_DA_Lin,
#            "PspB"= pspb_DA_Lin,
#            "PspC"= pspc_DA_Lin,
#            "LiaG"= liag_DA_Lin,
#            "LiaF"= liaf_DA_Lin,
#            "LiaI"= liai_DA_Lin
#            )
#
#
#     prot <- prot_da_lin %>% group_by(DomArch.norep)%>%
#       summarise(totalcount = sum(count)) %>%
#       filter(totalcount >1) %>% arrange(desc(totalcount))
#     cutoff_val <- prot[current_row,"totalcount"]
#   })
#
#   prot<- reactive({
#     prot_da_lin <- switch(input$linSelec,
#                           "PspA"= pspa_DA_Lin,
#                           "PspB"= pspb_DA_Lin,
#                           "PspC"= pspc_DA_Lin,
#                           "LiaG"= liag_DA_Lin,
#                           "LiaF"= liaf_DA_Lin,
#                           "LiaI"= liai_DA_Lin
#     )
#     prot <- prot_da_lin%>% group_by(DomArch.norep)%>%
#       summarise(totalcount = sum(count)) %>%
#       filter(totalcount >1) %>% arrange(desc(totalcount))
#   })
#   cutoff_val <- reactive({
#     prot[current_row,"totalcount"]}
#   )
#
#   observeEvent(input$cutoff, {
#     v$x <- input$cutoff
#   })
#   observeEvent(input$nr, {
#     v$x <- prot[input$nr,"totalcount"]
#   })
#
# observeEvent(input$cutoff,{
#   if(v$x != input$cutoff) updateSliderInput(session, "cutoff",value=cutoff_val)
#
#   if(v$x != input$nr) updateSliderInput(session,"nr",value =length(filter(prot,totalcount>=cutoff_val)))
#   }
#)




  #Renders the heatmap
  output$LinPlot <- renderPlot({
    req(credentials()$user_auth)
    if(input$DA_GC == "Domain Architecture"){
      switch(input$linSelec,
             "PspA" = lineage.DA.plot(pspa, filter(pspa_DA.cummulative,totalcount >= input$cutoff),"DomArch.norep", ""),
             "PspB" = lineage.DA.plot(pspb, filter(pspb_DA.cummulative,totalcount >= input$cutoff),"DomArch.norep", ""),
             "PspC" = lineage.DA.plot(pspc, filter(pspc_DA.cummulative,totalcount>= input$cutoff),"DomArch.norep", ""),
             "PspN" = lineage.DA.plot(pspc, filter(pspn_DA.cummulative,totalcount>= input$cutoff),"DomArch.norep", ""),
             "LiaF" = lineage.DA.plot(liai_liaf, filter(liai_liaf_DA.cummulative,totalcount>=input$cutoff),"DomArch.norep", ""),
             "Toast-rack" = lineage.DA.plot(toast_rack, filter(toast_rack_DA.cummulative,totalcount>=input$cutoff),"DomArch.norep", ""))
             #"LiaI" = lineage.DA.plot(pspa, pspa_DA.cummulative,"DomArch.norep", ""))
    }
    else{
      switch(input$linSelec,
             "PspA" = lineage.DA.plot(pspa, filter(pspa_GC.cummulative,totalcount>=input$cutoff),"GenContext.norep", ""),
             "PspB" = lineage.DA.plot(pspb, filter(pspb_GC.cummulative,totalcount>=input$cutoff),"GenContext.norep", ""),
             "PspC" = lineage.DA.plot(pspc, filter(pspc_GC.cummulative,totalcount>=input$cutoff),"GenContext.norep", ""),
             "PspN" = lineage.DA.plot(pspc, filter(pspn_GC.cummulative,totalcount>=input$cutoff),"GenContext.norep", ""),
             "LiaF" = lineage.DA.plot(liai_liaf, filter(liai_liaf_GC.cummulative,totalcount>=input$cutoff),"GenContext.norep", ""),
             "Toast-rack" = lineage.DA.plot(toast_rack, filter(toast_rack_GC.cummulative,totalcount>=input$cutoff),"GenContext.norep", ""))
    }
  }, height = 500)

  #Renders the datatable for the lineage counts
  output$LinTable <- DT::renderDT({
    req(credentials()$user_auth)
    paged_table(
    if(input$DA_GC == "Domain Architecture"){
      switch(input$linSelec,
             "PspA" = filter(pspa_DA.cummulative,totalcount >= input$cutoff),
             "PspB" = filter(pspb_DA.cummulative,totalcount >= input$cutoff),
             "PspC" = filter(pspc_DA.cummulative,totalcount >= input$cutoff),
             "PspN" = filter(pspn_DA.cummulative,totalcount >= input$cutoff),
             "LiaF" = filter(liai_liaf_DA.cummulative, totalcount >= input$cutoff),
             "Toast-rack" = filter(toast_rack_DA.cummulative, totalcount >= input$cutoff)
      )}
    else{
      switch(input$linSelec,
             "PspA" = filter(pspa_GC.cummulative,totalcount >= input$cutoff),
             "PspB" = filter(pspb_GC.cummulative,totalcount >= input$cutoff),
             "PspC" = filter(pspc_GC.cummulative,totalcount >= input$cutoff),
             "PspN" = filter(pspn_GC.cummulative,totalcount >= input$cutoff),
             "LiaF" = filter(liai_liaf_GC.cummulative,totalcount >= input$cutoff),
             "Toast-rack" = filter(toast_rack_GC.cummulative,totalcount >= input$cutoff)
      )
    }
  )
  }, extensions = c('FixedColumns',"FixedHeader"),
  options = list(pageLength = 15,
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))


  upset_vals <- reactiveValues(GC.DA.wc=pspa$GC.DA %>% words2wc(), DA.doms.wc=pspa$DA.doms %>% words2wc())
  observeEvent(input$linSelec,{
      upset_vals$GC.DA.wc <- switch(input$linSelec,
                              "PspA" = pspa$GC.DA %>% words2wc(),
                              "PspB" = pspb$GC.DA %>% words2wc(),
                              "PspC" = pspc$GC.DA %>% words2wc(),
                              "PspN" = pspn$GC.DA %>% words2wc(),
                              "LiaF" = liai_liaf$GC.DA %>% words2wc(),
                              "Toast-rack" = toast_rack$GC.DA %>% words2wc(),
                              "Liag" = liag$GC.DA %>% words2wc())
      upset_vals$DA.doms.wc <- switch(input$linSelec,
                                "PspA" = pspa$DA.doms %>% words2wc(),
                                "PspB" = pspb$DA.doms %>% words2wc(),
                                "PspC" = pspc$DA.doms %>% words2wc(),
                                "PspN" = pspn$DA.doms %>% words2wc(),
                                "LiaF" = liai_liaf$DA.doms %>% words2wc(),
                                "Toast-rack" = toast_rack$DA.doms %>% words2wc(),
                                "Liag" = liag$DA.doms %>% words2wc())
  }
  )

  #Renders the upsetPlot
  output$upsetP <- renderPlot({
    req(credentials()$user_auth)
    selected <- input$linSelec
    DA_or_GC <- input$DA_GC
    if(DA_or_GC == "Domain Architecture"){
      switch(selected,
             "PspA"= lineage.upset(pspa,input$cutoff, "da2doms",upset_vals$DA.doms.wc),
             "PspB"= lineage.upset(pspb,input$cutoff, "da2doms",upset_vals$DA.doms.wc),
             "PspC"= lineage.upset(pspc,input$cutoff, "da2doms",upset_vals$DA.doms.wc),
             "PspN"= lineage.upset(pspn,input$cutoff, "da2doms",pspn$DA.doms %>% words2wc()),
             "LiaF"= lineage.upset(liai_liaf,input$cutoff, "da2doms",upset_vals$DA.doms.wc),
             "Toast-rack"= lineage.upset(toast_rack,input$cutoff, "da2doms",upset_vals$DA.doms.wc)
      )
    }
    else{
      switch(selected,
             "PspA"= lineage.upset(pspa, input$cutoff, "gc2da",upset_vals$GC.DA.wc),
             "PspB"= lineage.upset(pspb, input$cutoff, "gc2da",upset_vals$GC.DA.wc),
             "PspC"= lineage.upset(pspc, input$cutoff, "gc2da",upset_vals$GC.DA.wc),
             "PspN"= lineage.upset(pspn, input$cutoff, "gc2da",pspn$GC.DA %>% words2wc()),
             "LiaF"= lineage.upset(liai_liaf, input$cutoff, "gc2da",upset_vals$GC.DA.wc),
             "Toast-rack"= lineage.upset(toast_rack, input$cutoff, "gc2da",upset_vals$GC.DA.wc)
      )
    }
  }, height = 550)


  #Render Wordcloud
  #reverse_operons probably not necessary anymore
  output$wordcloud <- renderPlot({
    req(credentials()$user_auth)
    if(input$DA_GC == "Genomic Context"){
      wordcloud(upset_vals$GC.DA.wc$words, upset_vals$GC.DA.wc$freq, min.freq = input$cutoff,colors = brewer.pal(8, "Spectral"))
    }
    else{
      wordcloud(upset_vals$DA.doms.wc$words, upset_vals$DA.doms.wc$freq, min.freq = input$cutoff,colors = brewer.pal(8, "Spectral"))
    }

  }, height = 550)



  output$msaTree <- renderUI({
    print(input$phylo)
    req(credentials()$user_auth)
    img(src="rseqtree.png", "data-zoom-image" ="rseqtree.png", height=1024,width=800)
  })

  output$msaPlot <- renderUI({
    req(credentials()$user_auth)
    tags$iframe(style="height:600px; width:100%", src="pspa_reduced.fasta.pdf", seamless=T)
    # mytest <-tags$iframe(src="www.rstudio.com",scrolling="yes",seamless=T, height=600, width=535)
    # print(mytest)
    # mytest
  })

  vals <- reactiveValues(btn = 0, tab = "home")


  observeEvent(input$myBtn,{
    if( input$phylo =="Tree" ){
      vals$btn <- 1
      vals$tab <- input$phylo
    }
  }
  )
  observeEvent(input$phylo,{
    if( input$phylo !="Tree"){
      vals$tab <- "notphylo"
      vals$btn <- 0
    }
  }
  )


  observe({
    if(vals$btn == 1 && input$mainTabs == "phylogeny"){
      session$sendCustomMessage(type = 'testmessage'
                                ,message = list(value="ZoomOn")
      )
    }
    else{
      session$sendCustomMessage(type = 'testmessage'
                                ,message = list(value="ZoomOff")
      )
    }
  })


  output$ParalogTable <- DT::renderDataTable({
    req(credentials()$user_auth)
    switch(input$alignSelec,
           "PspA"= find_paralogs(all%>% filter(Query=="pspa")))
  },extensions = c('FixedColumns'),
  options = list(pageLength = 10,
                 #The below line seems to disable other pages and the search bar
                 #dom = 't',
                 scrollX = TRUE,
                 paging=TRUE,
                 fixedHeader=TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 0)))






  #Reactive expresion to change file name depending on which protein is selected
  fileNam <- reactive({
    if(input$downloadType == "tsv"){
      paste(input$proSelec, ".txt", sep = "")
    }
    else  paste(input$proSelec, ".csv", sep = "")
  })

  #Downloads the data from datatable
  output$downloadData <- downloadHandler(
    filename = function() {
      fileNam()
    },
    content = function(file) {

      if(input$downloadType == "tsv"){
        write_tsv(pspTable(), file)
      }
      else  write.csv(pspTable(), file, row.names = FALSE)
    }
  )

  #Download Data for Lin Table
  output$downloadCounts <- downloadHandler(
    req(credentials()$user_auth),
    filename = function(){
      if(input$DA_GC == "Domain Architecture"){extension <- "-da_lin_counts"}
      else if(input$DA_GC == "Genomic Context"){ extension <- "-gc_lin_counts"}
      if(input$downloadType == "tsv"){
        paste(input$linSelec,extension,".txt", sep = "")
      }
      else{  paste(input$linSelec,extension , ".csv", sep = "")}
    },
    content = function(file){
      if(input$DA_GC == "Domain Architecture"){
        selected <- switch(input$linSelec,
                           "PspA" = pspa_totalC,
                           "PspB" = pspb_totalC,
                           "PspC" = pspc_totalc,
                           "LiaF" = liaf_totalC,
                           "LiaG" = liag_totalC,
                           "LiaI" = liai_totalC)
      }
      else if(input$DA_GC == "Genomic Context"){
        selected <- switch(input$linSelec,
                           "PspA" = pspa_cum,
                           "PspB" = pspb_cum,
                           "PspC" = pspc_cum,
                           "LiaF" = liaf_cum,
                           "LiaG" = liag_cum,
                           "LiaI" = liai_cum)
      }
      if(input$countDownloadType == "tsv"){
        write_tsv(selected, file )
      }
      else if(input$countDownloadType == "csv"){
        write.csv(selected,file)
      }
    }
  )
}

#Call to shiny app
shinyApp(ui = ui, server = server)


