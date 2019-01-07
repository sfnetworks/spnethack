#' @examples
#' library(sf)
#' sfn = sfnet(promenade)
#' plot(sfn)
#' plot(attr(sfn, "graph"))
sfnet = function(x) {
  xrnet = stplanr::SpatialLinesNetwork(x)
  attr(x, "graph") = xrnet@g
  class(x) = c("sfnet", class(x))
  x
}
