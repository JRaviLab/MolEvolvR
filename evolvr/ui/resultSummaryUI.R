
tabPanel(
        title = "Result Summary",
        value = "resultSummary",
        icon = icon("poll"),
        fluidPage(
                # uiOutput(outputId = "resultSummaryUI")
                # fluidRow(
                #   div(
                #     class = "cardRow", style = "background-color: blue; width: 100%; height: 500px; position: relative; margin: 0; padding: 0;",
                #     strong("Result Summary")
                #
                #
                # ),

                uiOutput(outputId = "rs_DomArch_ui"),

                uiOutput(outputId = "rs_GenContext_ui"),

                uiOutput(outputId = "rs_Phylogeny_ui")


                # fluidRow(
                #         column(width = 2, offset = 4,
                #                actionLink("rs2DomArch", label = h3("Domain Architecture")),
                #                p("Proximity Network")
                #         ),
                #         column(width = 6, offset = 0,
                #                visNetworkOutput(outputId = "rs_network")
                #         )
                # ),
                # fluidRow(
                #         column(width = 6, offset = 0,
                #                plotOutput(outputId = "rs_gcHeatmap")
                #         ),
                #         column(width = 2, offset = 0,
                #                actionLink("rs2GenContext", label = h3("Genomic Context")),
                #                p("Lineage Heatmap")
                #         )
                # ),
                # fluidRow(
                #         column(width = 2, offset = 4,
                #                actionLink("rs2Phylogeny", label = h3("Phylogeny")),
                #                p("Sunburst")
                #         ),
                #         column(width = 6, offset = 0,
                #                sunburstOutput(outputId = "rs_sunburst")
                #         )
                # )

        )
)

# full_data_results =