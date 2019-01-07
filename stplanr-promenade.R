rstp = SpatialLinesNetwork(sl = promenade_min)
rstp@sl$b = igraph::edge.betweenness(rstp@g)
lwd = rstp@sl$b / mean(rstp@sl$b) * 5
plot(rstp, lwd = lwd)
