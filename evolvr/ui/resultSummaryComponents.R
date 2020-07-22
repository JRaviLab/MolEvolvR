
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
