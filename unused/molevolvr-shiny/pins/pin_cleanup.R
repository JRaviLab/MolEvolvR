library(lubridate)
library(pins)

pin_date <- function(name, board)
{
  # Retrieve the date that the pin was last modified
  pin_recent_date <- pin_versions(name = name, board = board)$created[1]
  return(as_datetime(pin_recent_date))
}

clean_pins <- function(board, expiration_time, timetype = "hours")
{
  #' Remove pins that have been pinned for longer than the expiration time
  #'@author Samuel Chen
  #'@param board Pin Board to clean
  #'@param expiration_time The time limit that pins are allowed to remain
  #'@param timetype String of the unit of the expiration time. Can be either
  #'"hours", "minutes", "seconds", or "days". Defaults to hours
  #'@return Integer value of the number of pins removed

  # Convert expiration time to seconds
  expiration_time <- switch(timetype,
                            "hours" = expiration_time*3600,
                            "minutes" = expiration_time*60,
                            "days" = expiration_time * 24*3600,
                            "seconds" = expiration_time
  )

  expiration_time <- make_difftime(expiration_time, units = "hours")

  pin_names <- pin_find(board = "github")$name
  # pin_get(board = "github")

  current_time = as_datetime(Sys.time())


  removed = 0
  for(pin in pin_names)
  {
    if(current_time - pin_date(pin, board) > expiration_time )
    {
      pin_remove(name = pin, board = board)
      removed = removed + 1
    }

  }
  return(removed)

}

AnalyzePins <- function(board, timeBound = 5, timetype = "minutes" )
{
  # Convert expiration time to hours
  timeBound <- switch(timetype,
                            "hours" = timeBound*3600,
                            "minutes" = timeBound*60,
                            "days" = timeBound * 24*3600,
                            "seconds" = timeBound
  )

  timeBound <- make_difftime(timeBound, units = "hours")

  pin_names <- pin_find(board = board)$name
  # pin_get(board = "github")

  current_time = as_datetime(Sys.time())

  for(pin in pin_names)
  {
    postfix = substr(pin, 8, stop = nchar(pin))

    if(postfix == "AccNum2BLAST" & current_time - pin_date(pin, board) < timeBound )
    {
      print(pin)
    }

  }

}

