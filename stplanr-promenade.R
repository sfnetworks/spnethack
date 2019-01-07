rstp = SpatialLinesNetwork(sl = promenade_min)
rstp@sl$b = igraph::edge.betweenness(rstp@g)
lwd = rstp@sl$b / mean(rstp@sl$b) * 5
plot(rstp, lwd = lwd)

rtg = as_tbl_graph(x = rstp@g)

rtg_edges = activate(rtg,edges) %>% mutate(geometry = promenade_min$geometry)

st_as_sf(as.tibble(rtg_edges)) 
