The R community can lean on very powerful packages for spatial vector
data analysis ([sf](https://github.com/r-spatial/sf)) and network
analysis ([tidygraph](https://github.com/thomasp85/tidygraph)), that
both support the popular 'tidy' approach to data science. In sf, spatial
vector data are stored as objects of class `sf`, which are flat tables
with a list column that contains the geometry of the features. Tidygraph
is build on top of the widely-used
[igraph](https://github.com/igraph/igraph) package, and stores networks
in objects of class `tbl_graph`. A `tbl_graph` is an igraph-object, but
enables the user to manipulate both the edges and nodes elements as
being flat tables.

Despite the existence of sf and tidygraph, R seems to lack a general,
modern way to store networks whose nodes are embedded in space, i.e.
spatial networks. The [stplanr](https://github.com/ropensci/stplanr)
package contains the `SpatialLinesNetwork` class. This, however, is
based on [sp](https://github.com/edzer/sp/), a package for spatial data
analysis launched in 2005, that is used less since sf entered the stage.
The same yields for [spnetwork](https://github.com/edzer/spnetwork), a
package that combined sp and igraph. More recently,
[dodgr](https://github.com/ATFutures/dodgr) was created. This package
provides very fast analytical tools for spatial networks, but deals
solely with dual-weighted directed graphs, mainly used for calculating
shortest paths in street networks.

This blogpost presents a general approach to store spatial networks in a
tidy way by combining sf and tidygraph. At the same time, improvement
points to make this approach more convenient, are discussed.

Example 1: Fietstelweek
-----------------------

Fietstelweek is, translated from Dutch, the National Bicycle Count week.
It is a crowdsourced initiative to collect data on what routes cyclists
take to get from their origins to their destinations. The data were
collected between September 14th and 20th in 2015, and between September
19th and 25th in 2016, through a mobile app which makes use of the
cellphone GPS to track the trips. Around 50 thousand people
participated. The data are available on
[bikeprint.nl](http://www.bikeprint.nl/fietstelweek/) for download, and
include seperate spatial files for edges and nodes.

To illustrate how such a spatial network could be stored in R, we use
the Fietstelweek data for the municipality of Groningen as an example.

    library(piggyback)
    library(sf)

    pb_download("groningen_edges.gpkg")
    pb_download("groningen_nodes.gpkg")

    edges <- st_read("groningen_edges.gpkg")
    nodes <- st_read("groningen_nodes.gpkg")

Each node has a `POINT` geometry and several attributes: a unique ID
(`KNOOPNUMME`), the average waiting time at the node (`TIJD`), and the
number of routes that passed the node (`AANTAL`).

    nodes

    ## Simple feature collection with 7615 features and 7 fields
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: 720027 ymin: 7016407 xmax: 746520 ymax: 7032136
    ## epsg (SRID):    NA
    ## proj4string:    +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs
    ## First 10 features:
    ##    KNOOPNUMME       TIJD AANTAL  TIJD_OP  id code gemeentena
    ## 1    48260253  3.6481000     20 20.05680 142 0014  Groningen
    ## 2   972989406 19.5109247    146 39.41951 142 0014  Groningen
    ## 3    48260241  0.5067500      4 12.55825 142 0014  Groningen
    ## 4    48289812  1.1117500      4 12.12375 142 0014  Groningen
    ## 5   956932338 21.8750619    113 37.56719 142 0014  Groningen
    ## 6  1140905063  1.8569245     53 14.24666 142 0014  Groningen
    ## 7    48205957  4.3340000     82 19.82232 142 0014  Groningen
    ## 8  1307020232  0.7828000     10 15.51460 142 0014  Groningen
    ## 9    48281675  3.6616154     65 19.81608 142 0014  Groningen
    ## 10  571125378  0.9427255     51 14.24422 142 0014  Groningen
    ##                      geom
    ## 1  POINT (741183 7027692)
    ## 2  POINT (730615 7018718)
    ## 3  POINT (741396 7027677)
    ## 4  POINT (745784 7030454)
    ## 5  POINT (738479 7026407)
    ## 6  POINT (737646 7018928)
    ## 7  POINT (722372 7022711)
    ## 8  POINT (721538 7024870)
    ## 9  POINT (726748 7029664)
    ## 10 POINT (737877 7018747)

Each edge has a `MULTILINESTRING` geometry and several attributes: a
unique ID (`LINKNUMMER`), the ID's of the source and target nodes, the
functional type of the edge (`HIGHWAY`), the average speed, and the
number of routes using that edge (`INTENSITEI`).

    edges

    ## Simple feature collection with 13243 features and 11 fields
    ## geometry type:  MULTILINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 719848.4 ymin: 7016376 xmax: 747384.2 ymax: 7032180
    ## epsg (SRID):    NA
    ## proj4string:    +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs
    ## First 10 features:
    ##    LINKNUMMER      SOURCE     TARGET      HIGHWAY   SPEED INTENSITEI
    ## 1     1434443   918318147  429213139     cycleway 20.7575         60
    ## 2      424855  1104862685 1082653501  residential 15.0106         10
    ## 3      410048    48230404   48231375 unclassified 15.4530          2
    ## 4     2169959 -1880429983 -205866091      footway  0.0000          0
    ## 5      411767    48258309   48260439  residential  0.0000          0
    ## 6      410304    48198378   48200514  residential 14.7470         12
    ## 7      413144  1178147202   48174051  residential 19.4029        127
    ## 8     2524954   319036944  319036947      footway 20.2268          1
    ## 9     1717229  1183135174 1141319605     cycleway 17.2862        222
    ## 10    2460064  -812818258 -812818255      footway 16.7918          1
    ##    INTENSI_01 SNELHEID_R  id code gemeentena
    ## 1          25   0.932773 142 0014  Groningen
    ## 2          12   0.830246 142 0014  Groningen
    ## 3           4   0.956290 142 0014  Groningen
    ## 4          53   0.000000 142 0014  Groningen
    ## 5           8   0.000000 142 0014  Groningen
    ## 6          20   0.788605 142 0014  Groningen
    ## 7         102   0.916165 142 0014  Groningen
    ## 8           0   1.000000 142 0014  Groningen
    ## 9         184   0.817026 142 0014  Groningen
    ## 10          0   0.843472 142 0014  Groningen
    ##                              geom
    ## 1  MULTILINESTRING ((727715.4 ...
    ## 2  MULTILINESTRING ((726348.4 ...
    ## 3  MULTILINESTRING ((735131.7 ...
    ## 4  MULTILINESTRING ((733808.5 ...
    ## 5  MULTILINESTRING ((732116 70...
    ## 6  MULTILINESTRING ((732363.1 ...
    ## 7  MULTILINESTRING ((730488.9 ...
    ## 8  MULTILINESTRING ((729630.1 ...
    ## 9  MULTILINESTRING ((728993.5 ...
    ## 10 MULTILINESTRING ((728296.7 ...

The `sf` class subclasses a `data.frame`, and the tidygraph package has
a specific function to create a `tbl_graph` object from two data frames,
where one contains the nodes and the other the edges. That is, we can
easily construct a `tbl_graph` from the two `sf` objects that we created
above. Regarding the edges data frame, it is important that the columns
that specify the ID's of the source an target node are either the first
two columns of the data frame, or are named 'to' and 'from',
respectively. Both may be inconvenient and even unwanted, but for now,
it unfortunately is the only way.

    # Reorder the columns in the edges data frame
    (edges <- edges[, c(2, 3, 1, c(4:ncol(edges)))])

    ## Simple feature collection with 13243 features and 11 fields
    ## geometry type:  MULTILINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 719848.4 ymin: 7016376 xmax: 747384.2 ymax: 7032180
    ## epsg (SRID):    NA
    ## proj4string:    +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs
    ## First 10 features:
    ##         SOURCE     TARGET LINKNUMMER      HIGHWAY   SPEED INTENSITEI
    ## 1    918318147  429213139    1434443     cycleway 20.7575         60
    ## 2   1104862685 1082653501     424855  residential 15.0106         10
    ## 3     48230404   48231375     410048 unclassified 15.4530          2
    ## 4  -1880429983 -205866091    2169959      footway  0.0000          0
    ## 5     48258309   48260439     411767  residential  0.0000          0
    ## 6     48198378   48200514     410304  residential 14.7470         12
    ## 7   1178147202   48174051     413144  residential 19.4029        127
    ## 8    319036944  319036947    2524954      footway 20.2268          1
    ## 9   1183135174 1141319605    1717229     cycleway 17.2862        222
    ## 10  -812818258 -812818255    2460064      footway 16.7918          1
    ##    INTENSI_01 SNELHEID_R  id code gemeentena
    ## 1          25   0.932773 142 0014  Groningen
    ## 2          12   0.830246 142 0014  Groningen
    ## 3           4   0.956290 142 0014  Groningen
    ## 4          53   0.000000 142 0014  Groningen
    ## 5           8   0.000000 142 0014  Groningen
    ## 6          20   0.788605 142 0014  Groningen
    ## 7         102   0.916165 142 0014  Groningen
    ## 8           0   1.000000 142 0014  Groningen
    ## 9         184   0.817026 142 0014  Groningen
    ## 10          0   0.843472 142 0014  Groningen
    ##                              geom
    ## 1  MULTILINESTRING ((727715.4 ...
    ## 2  MULTILINESTRING ((726348.4 ...
    ## 3  MULTILINESTRING ((735131.7 ...
    ## 4  MULTILINESTRING ((733808.5 ...
    ## 5  MULTILINESTRING ((732116 70...
    ## 6  MULTILINESTRING ((732363.1 ...
    ## 7  MULTILINESTRING ((730488.9 ...
    ## 8  MULTILINESTRING ((729630.1 ...
    ## 9  MULTILINESTRING ((728993.5 ...
    ## 10 MULTILINESTRING ((728296.7 ...

    # Construct a tbl_graph
    graph <- tbl_graph(nodes = nodes, edges = edges)

Oops, something seems to go wrong!
