tabItem("phylogeny",
        fluidPage(
          column(width = 2,
                 selectInput(inputId =  "alignSelec", label = "Protein",
                             choices = c(
                               "LiaI-LiaF-TM.allfa50",
                               "LiaI-LiaF-TM_LiaFN.2",
                               "LiaI-LiaF-TM_LiaI.1",
                               "LiaI-LiaF-TM_PspC.3",
                               "PspA Only",
                               "PspA Snf7 Gismo",
                               "PspA Snf7",
                               "PspB Gismo",
                               "PspC Gismo",
                               "Snf7 Only",
                               "Toast-rack DUF2154-LiaF",
                               "Toast-rack DUF2807",
                               "Toast-rack DUF4097",
                               "Toast-rack PspC-Cterm"
                               )
                             , selected = "PspA Snf7")),
          column(width = 12,
                 tabsetPanel(
                   id= "phylo",
                   tabPanel("Tree", value="Tree",

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
)