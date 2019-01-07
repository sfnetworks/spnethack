## From dodgr
muenster = dodgr_streetnet('Muenster, DE')
promenade_dodgr = muenster %>% filter(name == 'Promenade')

promenade_graph = weight_streetnet(promenade_dodgr, wt_profile = 'bicycle')
# promenade_vert = dodgr_vertices(promenade_graph)

## Sample points along the route to create flows
from = sample(promenade_graph$from_id, size = 1000)
to = sample(promenade_graph$to_id, size = 1000)
flows = matrix(10*runif(
  length(from)*length(to)),
  nrow = length (from)
)

graph_f = dodgr_flows_aggregate(promenade_graph, from = from, to = to, flows = flows)
graph_undir = merge_directed_flows(graph_f)

dodgr_flowmap(graph_f, linescale = 5)
