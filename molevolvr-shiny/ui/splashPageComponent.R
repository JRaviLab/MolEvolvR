

splashUIComponent <- tagList(
  fluidRow(
    column(
    width = 4, offset = 4,
    basic_box("Query", image = "none", height = '20px'),
    # sub_box("Accession Number(s)", image_name = "none", button_name = "dAccNumBtn", color = "#9de19a"),
    # sub_box("FASTA", image_name = "Fasta_sequence", button_name = 'dFastaBtn', color = "#a4c5ea"),
    # sub_box("BLAST Output", image_name = "protein_blast_logo", button_name = "dBlastBtn", color = "#bca9e1"),
    # sub_box("InterproScan Output", image_name = "interproscan_logo", button_name = "dIprScanBtn", color = "#e7eca3"),
    sub_box("Accession Number(s)", image_name = "none", button_name = "dAccNumBtn", color = "#9de19a"),
    sub_box("FASTA", image_name = "none", button_name = 'dFastaBtn', color = "#a4c5ea"),
    sub_box("BLAST Output", image_name = "none", button_name = "dBlastBtn", color = "#bca9e1"),
    sub_box("InterproScan Output", image_name = "none", button_name = "dIprScanBtn", color = "#e7eca3"),

    sub_box("Full Data", image_name = "none", button_name = "dFullBtn", color = "#98a7f2")
    )
  ),
  fluidRow(
      column(width = 2, offset = 3, arrow_dl()),
      column(width = 2, offset = 0, arrow_d()),
      column(width = 2, offset = 0, arrow_dr())
  ),
  fluidRow(
    column(width = 4, offset = 0,
           diagram_box("Domain Architecture", image_name = "DomainArchitectureCartoon", button_name = "dDomArchBtn", description = "")
           ),
    column(width = 1, offset = 0, arrow_l()),
    column(width = 2, offset = 0,
           diagram_box("Homologs", image_name = "HomologyHeatmap", button_name = "dHomologBtn", description = "")
           ),
    column(width =1, offset = 0, arrow_r()),
    column(width = 4, offset = 0,
           diagram_box("Genomic Context", image_name = "GenomicContextCartoon", button_name = "dGenContextBtn", description = "")
           )
  ),
  fluidRow(
    column(width = 2, offset = 3, arrow_dr()),
    column(width = 2, offset = 0, arrow_d()),
    column(width = 2, offset = 0, arrow_dl())
  ),
  fluidRow(
    column(width = 4, offset = 4,
           diagram_box("Phylogeny", image_name = "phylo_tree", button_name = "dPhyloBtn", description = "")
           )
  )


)



# fluidRow(
#   column(width = 3, offset = 1,
#          basic_box("Query Protein", image = "none", description = "", height = '0px'),
#          sub_box("Protein Accession Number", image_name = "none", description = "", button_name = "dAccNumBtn"),
#          sub_box("Fasta", image_name = "none", description = "", button_name = "dFastaBtn")
#   ),
#   column(width = 1, offset = 0,
#          arrow_r()
#   ),
#   column(width = 3, offset = 0,
#          diagram_box("Molecular Characterization of Target Protein", image = "mol_characterization", button_name = "dMolTargetBtn", description = "")
#   ),
#   column(width = 1, offset = 0,
#          arrow_r()
#   ),
#   column(width = 3, offset = 0,
#          diagram_box("Identify Homologs", image = "none", button_name = "dHomoBtn", description = "")
#   )
# ),
# fluidRow(
#   column(width = 3, offset = 7, arrow_dl())
# ),
# fluidRow(
#   column(width = 3, offset = 5,
#          diagram_box("Molecular Characterization of Homologs", image = "mol_characterization", button_name = "dMolHomoBtn", description = "")
#   )
#   # column(width = 3, offset = 5,
#   #        diagram_elipse("Analysis", image = "none") ),
#   # column(width = 1, offset = 0, arrow_l()),
#   # column(width = 3, offset = 0,
#   #        diagram_box("Molecular Characterization of Homologs", image = "mol_characterization", button_name = "dMolHomoBtn", description = "Description") )
# ),
# fluidRow(
#   column(width = 2, offset = 4, arrow_dl()),
#   column(width = 1, offset = 0, arrow_d()),
#   column(width = 2, offset = 0, arrow_dr())
# ),
# fluidRow(
#
#   column(width = 3, offset = 2,
#          diagram_box("Phyletic Spread", image = "heatmap", button_name = "dPhyleticBtn", description = "")
#   ),
#   column(width = 3, offset = 0,
#          diagram_box("Proximity Network", image = "domain_network", button_name = "dNetworkBtn", description = "")
#   ),
#   column(width = 3, offset = 0,
#          diagram_box("Structure Based Sequence Alignment", image = "phylo_tree", button_name = "dAlignmentBtn", description = "")
#   )
# )