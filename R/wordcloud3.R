wordcloud3 <- function (data, size = 1, minSize = 0, gridSize = 0, fontFamily = "Segoe UI",
                        fontWeight = "bold", color = "random-dark", backgroundColor = "white",
                        minRotation = -pi/4, maxRotation = pi/4, shuffle = TRUE,
                        rotateRatio = 0.4, shape = "circle", ellipticity = 0.65,
                        widgetsize = NULL, figPath = NULL, hoverFunction = NULL) {
  if ("table" %in% class(data)) {
    dataOut = data.frame(name = names(data), freq = as.vector(data))
  }
  else {
    data = as.data.frame(data)
    dataOut = data[, 1:3]
    names(dataOut) = c("name", "freq", "label")
  }
  if (!is.null(figPath)) {
    if (!file.exists(figPath)) {
      stop("cannot find fig in the figPath")
    }
    spPath = strsplit(figPath, "\\.")[[1]]
    len = length(spPath)
    figClass = spPath[len]
    if (!figClass %in% c("jpeg", "jpg", "png", "bmp", "gif")) {
      stop("file should be a jpeg, jpg, png, bmp or gif file!")
    }
    base64 = base64enc::base64encode(figPath)
    base64 = paste0("data:image/", figClass, ";base64,",
                    base64)
  }
  else {
    base64 = NULL
  }
  weightFactor = size * 180/max(dataOut$freq)
  settings <- list(word = dataOut$name, freq = dataOut$freq, label = dataOut$label,
                   fontFamily = fontFamily, fontWeight = fontWeight, color = color,
                   minSize = minSize, weightFactor = weightFactor, backgroundColor = backgroundColor,
                   gridSize = gridSize, minRotation = minRotation, maxRotation = maxRotation,
                   shuffle = shuffle, rotateRatio = rotateRatio, shape = shape,
                   ellipticity = ellipticity, figBase64 = base64, hover = htmlwidgets::JS(hoverFunction))
  chart = htmlwidgets::createWidget("wordcloud2", settings,
                                    width = widgetsize[1], height = widgetsize[2], sizingPolicy = htmlwidgets::sizingPolicy(viewer.padding = 0,
                                                                                                                            browser.padding = 0, browser.fill = TRUE))
  htmlwidgets::onRender(chart, "function(el,x){\n                        console.log(123);\n                        if(!iii){\n                          window.location.reload();\n                          iii = False;\n\n                        }\n  }")
}


#'
#'
#'
#' getClickedWord <- function(WordInputId) {
#'   #OUPUT
#'   #       - referencing input in server will return a string of form word:freq (same as hover info shown in wordcloud; ie 'super:32')
#'   shiny::tags$script(shiny::HTML(
#'     "$(document).on('click', '#canvas', function(evt) {",
#'     'var id = evt.target.nextElementSibling.firstChild.id;',
#'     'word = document.getElementById(id).innerHTML;',
#'     'console.log(id);',
#'     "Shiny.onInputChange(id.replace('wcSpan','_clicked_word'), word);",
#'     #sprintf("Shiny.onInputChange(id'_%s', word);",WordInputId),
#'     "});"
#'   ))
#' }
#'
#'
#'
#' #' @rdname wordcloud2-shiny
#' #' @export
#' wordcloud2Output <- function(outputId, width = "100%", height = "400px") {
#'   widget_out <- htmlwidgets::shinyWidgetOutput(outputId, "wordcloud2", width, height, package = "wordcloud2")
#'
#'   shiny::div(getClickedWord(), widget_out)
#' }
#'
#' #' @rdname wordcloud2-shiny
#' #' @export
#' renderWordcloud2 <- function(expr, env = parent.frame(), quoted = FALSE) {
#'   if (!quoted) { expr <- substitute(expr) } # force quoted
#'   htmlwidgets::shinyRenderWidget(expr, wordcloud2Output, env, quoted = TRUE)
#' }