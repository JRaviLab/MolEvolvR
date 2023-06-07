library(pins)
library(tidyverse)

GeneratePinName <- function(length = 6, postfix = "", board = "github")
{
  name = ""
  # Technically can infinitely loop, but only if a lot of pins generated at length 1 or 2
  while(name == "" || nrow(pin_find(name = name, board = board)) != 0)
  {
    # alpha num has 62 options, so there are 62^length possibilities for names!
    alphaNum <- c(letters, LETTERS, 0:9 )
    # Don't start or end with hyphen/underscore
    name <- alphaNum[runif(n = length, min = 1,max =  length(alphaNum))] %>% paste0(collapse = "")
    name <- paste0(name, postfix)
  }
  return(name)
}

