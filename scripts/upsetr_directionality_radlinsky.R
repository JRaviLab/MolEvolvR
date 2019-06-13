library(UpSetR)
test <- as.data.frame(matrix(rnorm(1:12), nrow = 4))
colnames(test) <- c("A", "B", "C")
test_up <- test
test_up[test >= 0] <- 1
test_up[test_up != 1] <- 0
test_up$Gene <- paste0("Gene", 1:nrow(test), "_UP")
test_down <- test
test_down[test < 0] <- 1
test_down[test_down != 1] <- 0
test_down$Gene <- paste0("Gene", 1:nrow(test), "_DOWN")
test <- rbind(test_up, test_down)


upordown <- function(row, direction) {
  gene <- row["Gene"]
  if (grepl(x = gene, pattern = direction)) {
    newData <- T
  }
  else {
    newData <- F
  }
}


metadata <- data.frame(
  c("A", "B", "C"),
  as.numeric(apply(test_up[, 1:3], 2, sum))
)
colnames(metadata) <- c(
  "sets",
  "NumberUP"
)


upset(test,
      sets = c("A", "B", "C"), set.metadata = list(
        data = metadata,
        plots = list(
          list(type = "hist", 
               column = "NumberUP", 
               assign = 20, # defines width of the meta-data histogram
               colors = "red")
        )
      ),
      queries = list(list(query = upordown, params = list("_UP"), color = "red", active = TRUE))
)