library(tidyverse)

deltaBlast <- function(deltaBlast_path, db_search_path, db, query, evalue = "1e-5", out, num_alignments, num_threads = 1)
{
  start = Sys.time()

  system(paste0("export BLASTDB=/",db_search_path ))

  system2(command = deltaBlast_path,
          args = c("-db", db,
                   "-query", query,
                   "-evalue", evalue,
                   "-out", out,
                   "-num_threads", num_threads,
                   "-num_alignments", num_alignments
                   #   ,"-outfmt", outfmt
          )
  )
  print(Sys.time()-start)
}


rpsBlast <- function(rpsBlast_path,db_search_path,  db, query, evalue = "1e-5", out, num_threads = 1)
{
  start = Sys.time()
  system(paste0("export BLASTDB=/",db_search_path))
  system2(command = rpsBlast_path,
          args = c("-db", db,
                   "-query", query,
                   "-evalue", evalue,
                   "-out", out,
                   "-num_threads", num_threads
                   #                  , "-outfmt", outfmt
          )
  )
  print(Sys.time()-start)
}
