
query_LinPlot_txt <- div(class = "legend",
                         tags$span(
                           strong("The phyletic spread of query proteins."),
                           p(" The color gradient in the heatmap represents the number of homologs within each lineage. Rows: Query proteins/domains. Columns: Key lineages within the three kingdoms of life.")
                         )
                         )

DA_Network_txt <- div(class = "legend",
                      strong("The Domain Proximity Network"),
                      p("Domains (nodes) that co-occur within a protein/domain architecture are connected (edges). Node size and edge weight correspond to the frequency of occurrence.
                        ")
                      )

DA_Wordcloud_txt <- div(class = "legend",
                        strong("Domain Wordcloud"),
                        p("The commonly occurring domain partners are shown by frequency.")
                        )

DA_LinTable_txt <- div( class = "legend",
  p("Summary of domain architectures across lineages.")
)

DA_upset_txt <-
  tagList( strong("Frequencies of Domains and Domain Architectures"),
                    p("Blue histogram: Distribution of the constituent domains underlying all homologs. Combination matrix: Various configurations of the constituent domains that come together in the domain architectures. Red histogram: Frequency of occurrences of the indicated domain architectures.")
)

DA_LineageHeatmap_txt <- div(
  strong("The phyletic spread of the domain architectures."),
  p("The color gradient in the heatmap
                                                                 represents the number of homologs identified within each lineage. Rows: Predominant domain
                                                                 architectures. Columns: Key lineages within the three kingdoms of life."))

GC_LinTable_txt <- div( class = "legend",
  p("Summary of genomic contexts across lineages.")
)

GC_LineageHeatmap_txt <- div(
  class = "legend",
  tags$span(
    strong("The phyletic spread of genomic contexts. "),
    p("The color gradient in the heatmap represents the number of homologs with the indicated genomic context within each lineage. Rows: Predominant genomic contexts. Columns: Key lineages within the three kingdoms of life.")
  )
)


GC_upset_txt <- div(class = "legend",strong("Frequencies of Domain Architectures and Genomic Contexts."),
                    p("Blue histogram: Distribution of the constituent domain architectures. Combination matrix: Various configurations of the domain architectures that come together in the genomic context. Red histogram: Frequency of occurrences of the indicated genomic contexts.")
)

GC_Wordcloud_txt <- div(class = "legend",
                        strong("Genomic Neighborhood Wordcloud"),
                        p("The common genomic neighbors of query proteins are shown by frequency.")
)

tree_txt <- div(class = "legend",
                tags$span(strong("Phylogenetic Tree. "), p("Constructed based on a multiple sequence alignment performed using representative homologs across all the major kingdoms/phyla. Colors: key phyla. Tree leaves: labeled by lineage, species (three-letter abbreviation), and accession numbers."))
                )
sunburst_txt <- div(class = "legend",
                    tags$span(strong("Sunburst plots for lineages."), p(" The interactive plot shows the phyletic spread of the query protein across bacteria, archaea, and eukaryota. Toggle checkbox to see the legend."))
                    )
msa_txt <- div(class = "legend",
               style = "display: inline;",
               tags$span(strong("Multiple sequence alignment"), p("of representative homologs across distinct lineages."))
               )

paralog_txt <- div(class = "legend",
                   style = "display: inline;",
                   strong("Paralog Table"),
                   p("Species containing multiple paralogs are summarized along with detailed protein-level data.")
                   )
