discrete_palettes <- c( # arc>bac>euk>vir
  RColorBrewer::brewer.pal(8, "Set2"),                         # 1-8
  'darksalmon', 'rosybrown', 'navy', 'firebrick3',             # 9-12
  'yellowgreen', 'deepskyblue3', 'orange', 'gold3',            # 13-16
  'royalblue1', 'cyan4', 'darkmagenta', 'khaki',               # 18-20
  'salmon', 'darkolivegreen3', 'powderblue', 'mediumpurple1',  # 21-24
  RColorBrewer::brewer.pal(7, "Dark2")                         # 25-31
)


lineages <- c("Vir", "Euk", "Bacteria", "B.thermot",
             "B.thermob", "B.teneri", "B.spiro", "B.PVC",
             "B.proteo", "B.nitro", "B.igna", "B.gemma",
             "B.fuso", "B.firmi", "B.FCB", "B.elusi",
             "B.dictyo", "B.deino.therm", "B.cyano" , "B.chloroflexi",
             "B.chlorobi", "B.caldit", "B.caldis", "B.arma",
             "B.aqui", "B.actino", "B.acido", "Archaea",
             "A.TACK", "A.eury", "A.asgard")



par(mar = rep(0, 4))
pie(rep(1, length(discrete_palettes)), angle=0,
    col=discrete_palettes,
    labels=lineages)
