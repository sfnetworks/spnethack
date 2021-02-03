Spatial networks in R with sf and tidygraph
================
Lucas van der Meer, Robin Lovelace & Lorena Abad
September 25, 2019

## Introduction

Street networks, shipping routes, telecommunication lines, river
bassins. All examples of spatial networks: organized systems of nodes
and edges embedded in space. For most of them, these nodes and edges can
be associated with geographical coordinates. That is, the nodes are
geographical points, and the edges geographical lines.

Such spatial networks can be analyzed using graph theory. Not for
nothing, Leonhard Eulers famous work on the
<a href="https://www.mathsisfun.com/activity/seven-bridges-konigsberg.html" target="_blank">Seven Bridges of Köningsberg</a>,
which laid the foundations of graph theory and network analysis, was in
essence a spatial problem.

In R, there are advanced, modern tools for both the analysis of spatial
data and networks. Furthermore, several packages have been developed
that cover (parts of) spatial network analysis. As stated in this
<a href="https://github.com/r-spatial/sf/issues/966" target="_blank">github issue</a>
and this
<a href="https://twitter.com/zevross/status/1089908839816794118" target="_blank">tweet</a>,
concensus on how best to represent and analyse spatial networks has
proved elusive.

This blogpost demonstrates an approach to spatial networks that starts
with a set of geographic lines, and leads to an object ready to be used
for network analysis. Along the way, we will see that there are still
steps that need to be taken before the process of analyzing spatial
networks in R is user friendly and efficient.

## Existing R packages for spatial networks

Although R was originally designed as a language for statistical
computing, an active ‘R-spatial’ ecosystem has evolved. Powerful and
high performance packages for spatial data analysis have been developed,
thanks largely to interfaces to mature C/C++ libraries such as GDAL,
GEOS and PROJ, notably in the package
<a href="https://github.com/r-spatial/sf" target="_blank">sf</a> (see
<a href="https://geocompr.robinlovelace.net/intro.html#the-history-of-r-spatial" target="_blank">section 1.5</a>
of Geocomputation with R for a brief history). Likewise, a number of
packages for graph representation and analysis have been developed,
notably
<a href="https://github.com/thomasp85/tidygraph" target="_blank">tidygraph</a>,
which is based on
<a href="https://igraph.org/" target="_blank">igraph</a>.

Both sf and tidygraph support the `tibble` class and the broader ‘tidy’
approach to data science, which involves data processing pipelines, type
stability and a convention of representing everything as a data frame
(well a `tibble`, which is a data frame with user friendly default
settings). In sf, this means storing spatial vector data as objects of
class `sf`, which are essentially the same as a regular data frame (or
tibble), but with an additional ‘sticky’ list column containing a
geometry for each feature (row), and attributes such as bounding box and
CRS. Tidygraph stores networks in objects of class `tbl_graph`. A
`tbl_graph` is an `igraph` object, but enables the user to manipulate
both the edges and nodes elements as if they were data frames also.

Both sf and tidygraph are relatively new packages (first released on
CRAN in 2016 and 2017, respectively). It is unsurprising, therefore,
that they have yet to be combined to allow a hybrid, tibble-based
representation of spatial networks.

Nevertheless, a number of approaches have been developed for
representing spatial networks, and some of these are in packages that
have been published on CRAN.
<a href="https://github.com/ropensci/stplanr" target="_blank">stplanr</a>,
for instance, contains the `SpatialLinesNetwork` class, which works with
both the <a href="https://github.com/edzer/sp/" target="_blank">sp</a>
(a package for spatial data analysis launched in 2005) and sf packages.
<a href="https://github.com/ATFutures/dodgr" target="_blank">dodgr</a>
is a more recent package that provides analytical tools for street
networks, with a focus on directed graphs (that can have
direction-dependent weights, e.g. representing a one-way street). Other
packages seeking to implement spatial networks in R include
<a href="https://github.com/edzer/spnetwork" target="_blank">spnetwork</a>,
a package that defined a class system combining sp and igraph, and
<a href="https://cran.r-project.org/web/packages/shp2graph/index.html" target="_blank">shp2graph</a>,
which provides tools to switch between sp and igraph objects.

Each package has its merits that deserve to be explored in more detail
(possibly in a subsequent blog post). The remainder of this post
outlines an approach that combines `sf` and `igraph` objects in a
`tidygraph` object.

## Set-up

The following code chunk will install the packages used in this post:

``` r
# We'll use remotes to install packages, install it if needs be:
if(!"remotes" %in% installed.packages()) {
  install.packages("remotes")
}

cran_pkgs = c(
  "sf",
  "tidygraph",
  "sfnetworks",
  "osmdata",
  "tidyverse",
  "ggplot2",
  "units",
  "tmap",
)

remotes::install_cran(cran_pkgs)
```

``` r
library(sf)
library(tidygraph)
library(sfnetworks)
library(tidyverse)
library(ggplot2)
library(units)
library(tmap)
library(osmdata)
```

## Getting the data

As an example, we use the street network of the city center of Münster,
Germany. We will get the data from OpenStreetMap. Packages like `dodgr`
have optimized their code for such data, however considering that we
want to showcase this workflow for any source of data, we will generate
an object of class `sf` containing only `LINESTRING` geometries.
However, streets that form loops, are returned by `osmdata` as polygons,
rather than lines. These, we will convert to lines, using the
`osm_poly2line` function. One additional variable, the type of street,
is added to show that the same steps can be used for `sf` objects that
contain any number of additional variables.

``` r
muenster <- opq(bbox =  c(7.61, 51.954, 7.636, 51.968)) %>% 
  add_osm_feature(key = 'highway') %>% 
  osmdata_sf() %>% 
  osm_poly2line()

muenster_center <- muenster$osm_lines %>% 
  select(highway)
```

``` r
muenster_center
```

    ## Simple feature collection with 2365 features and 1 field
    ## geometry type:  LINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ## geographic CRS: WGS 84
    ## First 10 features:
    ##             highway                       geometry
    ## 4064462   secondary LINESTRING (7.624554 51.955...
    ## 4064463   secondary LINESTRING (7.626498 51.956...
    ## 4064467 residential LINESTRING (7.630898 51.955...
    ## 4064474   secondary LINESTRING (7.61972 51.9554...
    ## 4064476   secondary LINESTRING (7.619844 51.954...
    ## 4064482    tertiary LINESTRING (7.616402 51.957...
    ## 4064485     service LINESTRING (7.63275 51.9603...
    ## 4984982   secondary LINESTRING (7.614167 51.967...
    ## 4985138    cycleway LINESTRING (7.61525 51.9673...
    ## 4985140 residential LINESTRING (7.616774 51.968...

``` r
ggplot(data = muenster_center) + geom_sf()
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## From sf to sfnetwork: so much has changed!

If you read our previous blogpost, you will realize we went through a
long list of steps to come up with an object of class `tbl_graph`. Since
that blogpost, the `sfnetworks` package has made this whole process
significantly easier. Let’s look now at our simplified stepwise
approach!

### Step 1: Construct the network

``` r
muenster_center_comb = muenster_center %>% 
  st_combine() %>% 
  st_node() %>% 
  st_cast('LINESTRING') 
```

    ## Warning in st_node.sfc(.): st_node may not give correct results for longitude/
    ## latitude data

``` r
muenster_center_comb2 = muenster_center_comb %>% 
  st_as_sf() %>% 
  st_join(muenster_center, join = st_covered_by)
```

    ## although coordinates are longitude/latitude, st_covered_by assumes that they are planar

``` r
graph = as_sfnetwork(muenster_center_comb2, directed = F)
graph
```

    ## # A sfnetwork with 3480 nodes and 4861 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An undirected multigraph with 40 components with spatially explicit edges
    ## #
    ## # Node Data:     3,480 x 1 (active)
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##              geometry
    ##           <POINT [°]>
    ## 1 (7.624554 51.95561)
    ## 2 (7.625074 51.95564)
    ## 3 (7.625275 51.95566)
    ## 4 (7.626498 51.95617)
    ## 5  (7.626103 51.9561)
    ## 6 (7.625926 51.95595)
    ## # ... with 3,474 more rows
    ## #
    ## # Edge Data:     4,861 x 4
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##    from    to highway                                                   geometry
    ##   <int> <int> <chr>                                             <LINESTRING [°]>
    ## 1     1     2 secondary (7.624554 51.95561, 7.624838 51.95562, 7.624961 51.9556~
    ## 2     2     3 secondary                   (7.625074 51.95564, 7.625275 51.95566)
    ## 3     4     5 secondary (7.626498 51.95617, 7.626437 51.95617, 7.626377 51.9561~
    ## # ... with 4,858 more rows

### Step 2: Clean the network

To perform network analysis, we need a network with a clean topology. In
theory, the best way to clean up the network topology is by manual
editing, but this can be very labour intensive and time consuming,
mainly for large networks.

In sfnetworks, we have implemented several cleaning functions as spatial
morphers. To learn more about spatial morphers and what they are, see
the [dedicated
vignette](https://luukvdmeer.github.io/sfnetworks/articles/morphers.html)
for that.

``` r
(graph = graph %>%
  convert(to_spatial_subdivision, .clean = T))
```

    ## Warning: to_spatial_subdivision assumes attributes are constant over geometries

    ## # A sfnetwork with 3480 nodes and 4861 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An undirected multigraph with 40 components with spatially explicit edges
    ## #
    ## # Node Data:     3,480 x 1 (active)
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##              geometry
    ##           <POINT [°]>
    ## 1 (7.624554 51.95561)
    ## 2 (7.625074 51.95564)
    ## 3 (7.625275 51.95566)
    ## 4 (7.626498 51.95617)
    ## 5  (7.626103 51.9561)
    ## 6 (7.625926 51.95595)
    ## # ... with 3,474 more rows
    ## #
    ## # Edge Data:     4,861 x 4
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##    from    to highway                                                   geometry
    ##   <int> <int> <chr>                                             <LINESTRING [°]>
    ## 1     1     2 secondary (7.624554 51.95561, 7.624838 51.95562, 7.624961 51.9556~
    ## 2     2     3 secondary                   (7.625074 51.95564, 7.625275 51.95566)
    ## 3     4     5 secondary (7.626498 51.95617, 7.626437 51.95617, 7.626377 51.9561~
    ## # ... with 4,858 more rows

``` r
(graph = graph %>% 
  convert(to_spatial_smooth, .clean = T))
```

    ## # A sfnetwork with 2744 nodes and 4125 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An undirected multigraph with 20 components with spatially explicit edges
    ## #
    ## # Node Data:     2,744 x 1 (active)
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##              geometry
    ##           <POINT [°]>
    ## 1 (7.624554 51.95561)
    ## 2 (7.625074 51.95564)
    ## 3 (7.626498 51.95617)
    ## 4  (7.626103 51.9561)
    ## 5 (7.625926 51.95595)
    ## 6 (7.625906 51.95573)
    ## # ... with 2,738 more rows
    ## #
    ## # Edge Data:     4,125 x 4
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##    from    to highway                                                   geometry
    ##   <int> <int> <chr>                                             <LINESTRING [°]>
    ## 1     1     2 secondary (7.624554 51.95561, 7.624838 51.95562, 7.624961 51.9556~
    ## 2     3     4 secondary (7.626498 51.95617, 7.626437 51.95617, 7.626377 51.9561~
    ## 3     4     5 secondary (7.626103 51.9561, 7.626058 51.95607, 7.626017 51.95605~
    ## # ... with 4,122 more rows

``` r
(graph = graph %>% 
  convert(to_spatial_simple, .clean = T))
```

    ## # A sfnetwork with 2744 nodes and 4002 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An undirected simple graph with 20 components with spatially explicit edges
    ## #
    ## # Node Data:     2,744 x 1 (active)
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##              geometry
    ##           <POINT [°]>
    ## 1 (7.624554 51.95561)
    ## 2 (7.625074 51.95564)
    ## 3 (7.626498 51.95617)
    ## 4  (7.626103 51.9561)
    ## 5 (7.625926 51.95595)
    ## 6 (7.625906 51.95573)
    ## # ... with 2,738 more rows
    ## #
    ## # Edge Data:     4,002 x 4
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##    from    to highway                                                   geometry
    ##   <int> <int> <chr>                                             <LINESTRING [°]>
    ## 1     1     2 secondary  (7.624554 51.95561, 7.624838 51.95562, 7.624961 51.955~
    ## 2     1   925 residenti~                  (7.624559 51.95567, 7.624554 51.95561)
    ## 3     1   953 secondary                   (7.623884 51.95563, 7.624554 51.95561)
    ## # ... with 3,999 more rows

## Combining the best of both worlds

Having the network stored in the tbl\_graph structure, with a geometry
list column for both the edges and nodes, enables us to combine the wide
range of functionalities in sf and tidygraph, in a way that fits neatly
into the tidyverse.

With the `activate()` verb, we specify if we want to manipulate the
edges or the nodes. Then, most dplyr verbs can be used in the familiar
way, also when directly applied to the geometry list column. For
example, we can add a variable describing the length of each edge,
which, later, we use as a weight for the edges.

``` r
graph <- graph %>%
  activate(edges) %>%
  mutate(length = edge_length())

graph
```

    ## # A sfnetwork with 2744 nodes and 4002 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An undirected simple graph with 20 components with spatially explicit edges
    ## #
    ## # Edge Data:     4,002 x 5 (active)
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##    from    to highway                                          geometry   length
    ##   <int> <int> <chr>                                    <LINESTRING [°]>      [m]
    ## 1     1     2 secondary (7.624554 51.95561, 7.624838 51.95562, 7.62496~ 35.8537~
    ## 2     1   925 resident~          (7.624559 51.95567, 7.624554 51.95561)  5.8968~
    ## 3     1   953 secondary          (7.623884 51.95563, 7.624554 51.95561) 46.0666~
    ## 4     1  1713 resident~          (7.624554 51.95561, 7.624553 51.95557)  4.5289~
    ## 5     2   935 service             (7.625074 51.95564, 7.625095 51.9556)  4.5754~
    ## 6     2  2309 <NA>      (7.625074 51.95564, 7.625275 51.95566, 7.62564~ 39.8855~
    ## # ... with 3,996 more rows
    ## #
    ## # Node Data:     2,744 x 1
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##              geometry
    ##           <POINT [°]>
    ## 1 (7.624554 51.95561)
    ## 2 (7.625074 51.95564)
    ## 3 (7.626498 51.95617)
    ## # ... with 2,741 more rows

With one flow of pipes, we can ‘escape’ the graph structure, turn either
the edges or nodes back into real `sf` objects, and, for example,
summarise the data based on a specific variable.

``` r
graph %>%
  activate(edges) %>%
  st_as_sf() %>%
  group_by(highway) %>%
  summarise(length = sum(length))
```

    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar
    ## although coordinates are longitude/latitude, st_union assumes that they are planar

    ## Simple feature collection with 14 features and 2 fields
    ## geometry type:  MULTILINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ## geographic CRS: WGS 84
    ## # A tibble: 14 x 3
    ##    highway         length                                               geometry
    ##    <chr>              [m]                                  <MULTILINESTRING [°]>
    ##  1 cycleway     5193.300~ ((7.61525 51.96736, 7.615165 51.96744, 7.615148 51.96~
    ##  2 footway     35269.647~ ((7.616968 51.95749, 7.616983 51.95756), (7.617861 51~
    ##  3 path        19689.760~ ((7.625926 51.95595, 7.62582 51.95588, 7.625795 51.95~
    ##  4 pedestrian   7232.768~ ((7.626354 51.96529, 7.626341 51.96521), (7.634264 51~
    ##  5 residential 19170.934~ ((7.624559 51.95567, 7.624554 51.95561), (7.624554 51~
    ##  6 secondary    4509.889~ ((7.624554 51.95561, 7.624838 51.95562, 7.624961 51.9~
    ##  7 secondary_~   344.213~ ((7.614357 51.9673, 7.614268 51.96727, 7.614167 51.96~
    ##  8 service     20179.297~ ((7.625074 51.95564, 7.625095 51.9556), (7.630898 51.~
    ##  9 steps         256.964~ ((7.616523 51.95744, 7.616562 51.95755), (7.618848 51~
    ## 10 tertiary     2907.179~ ((7.619573 51.95543, 7.619585 51.95543, 7.61972 51.95~
    ## 11 tertiary_l~    43.856~ ((7.623719 51.96629, 7.623699 51.96623, 7.623664 51.9~
    ## 12 track         356.181~ ((7.613738 51.96494, 7.61359 51.96493, 7.613383 51.96~
    ## 13 unclassifi~   553.101~ ((7.631481 51.95742, 7.631599 51.95738), (7.631658 51~
    ## 14 <NA>        28597.482~ ((7.625074 51.95564, 7.625275 51.95566, 7.625643 51.9~

Switching back to `sf` objects is useful as well when plotting the
network, in a way that preserves its spatial properties.

``` r
ggplot() +
  geom_sf(data = st_as_sf(graph, "edges")) + 
  geom_sf(data = st_as_sf(graph, "nodes"), size = 0.5)
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Or, alternatively, in only a few lines of code, plot the network as an
interactive map. On this page, the interactive map might show as an
image, but
<a href="https://luukvdmeer.github.io/spnethack/imap.html" target="_blank">here</a>,
you should be able to really interact with it!

``` r
tmap_mode('view')

tm_shape(st_as_sf(graph, "edges")) +
  tm_lines() +
tm_shape(st_as_sf(graph, "nodes")) +
  tm_dots() +
tmap_options(basemaps = 'OpenStreetMap')
```

<!--html_preserve-->

<div id="htmlwidget-bdc34cf187cc7cd6d354" class="leaflet html-widget"
style="width:672px;height:480px;">

</div>

<!--/html_preserve-->

All nice and well, but these are not things that we necessarily need the
graph representation for. The added value of tidygraph, is that it opens
the door to the functions of the igraph library, all specifically
designed for network analysis, and enables us to use them inside a
‘tidy’ workflow. To cover them all, we would need to write a book, but
let’s at least show a few examples below.

### Centrality measures

Centraltity measures describe the importances of nodes in the network.
The simplest of those measures is the degree centrality: the number of
edges connected to a node. Another example is the betweenness
centrality, which, simply stated, is the number of shortest paths that
pass through a node. In tidygraph, we can calculate these and many other
centrality measures, and simply add them as a variable to the nodes.

The betweenness centrality can also be calculated for edges. In that
case, it specifies the number of shortest paths that pass through an
edge.

``` r
graph <- graph %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree()) %>%
  mutate(betweenness = centrality_betweenness(weights = length)) %>%
  activate(edges) %>%
  mutate(betweenness = centrality_edge_betweenness(weights = length))

graph
```

    ## # A sfnetwork with 2744 nodes and 4002 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An undirected simple graph with 20 components with spatially explicit edges
    ## #
    ## # Edge Data:     4,002 x 6 (active)
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##    from    to highway                              geometry   length betweenness
    ##   <int> <int> <chr>                        <LINESTRING [°]>      [m]       <dbl>
    ## 1     1     2 seconda~ (7.624554 51.95561, 7.624838 51.955~ 35.8537~        7908
    ## 2     1   925 residen~ (7.624559 51.95567, 7.624554 51.955~  5.8968~       62916
    ## 3     1   953 seconda~ (7.623884 51.95563, 7.624554 51.955~ 46.0666~        6709
    ## 4     1  1713 residen~ (7.624554 51.95561, 7.624553 51.955~  4.5289~       50040
    ## 5     2   935 service  (7.625074 51.95564, 7.625095 51.955~  4.5754~        7111
    ## 6     2  2309 <NA>     (7.625074 51.95564, 7.625275 51.955~ 39.8855~        1878
    ## # ... with 3,996 more rows
    ## #
    ## # Node Data:     2,744 x 3
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.603762 ymin: 51.94823 xmax: 7.645245 ymax: 51.97241
    ##              geometry degree betweenness
    ##           <POINT [°]>  <dbl>       <dbl>
    ## 1 (7.624554 51.95561)      4       62444
    ## 2 (7.625074 51.95564)      3        7106
    ## 3 (7.626498 51.95617)      3      142779
    ## # ... with 2,741 more rows

``` r
ggplot() +
  geom_sf(data = st_as_sf(graph, "edges"), col = 'grey50') + 
  geom_sf(data = st_as_sf(graph, "nodes"), aes(col = betweenness, size = betweenness)) +
  scale_colour_viridis_c(option = 'inferno') +
  scale_size_continuous(range = c(0,4))
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggplot() +
  geom_sf(data = st_as_sf(graph, "edges"), 
          aes(col = betweenness, size = betweenness)) +
  scale_colour_viridis_c(option = 'inferno') +
  scale_size_continuous(range = c(0,4))
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

### Shortest paths

A core part of spatial network analysis is generally finding the path
between two nodes that minimizes either the travel distance or travel
time. In igraph, there are several functions that can be used for this
purpose, and since a `tbl_graph` is just a subclass of an `igraph`
object, we can directly input it into every function in the igraph
package.

The function `distances`, for example, returns a numeric matrix
containing the distances of the shortest paths between every possible
combination of nodes. It will automatically choose a suitable algorithm
to calculate these shortest paths.

``` r
distances <- st_network_cost(graph, weights = "length")

distances[1:5, 1:5]
```

    ##           [,1]      [,2]      [,3]      [,4]      [,5]
    ## [1,]   0.00000  35.85371 159.67782 130.62567 109.90021
    ## [2,]  35.85371   0.00000 126.03936  96.98722  76.26176
    ## [3,] 159.67782 126.03936   0.00000  29.05214  49.77761
    ## [4,] 130.62567  96.98722  29.05214   0.00000  20.72546
    ## [5,] 109.90021  76.26176  49.77761  20.72546   0.00000

The function ‘shortest\_paths’ not only returns distances, but also the
indices of the nodes and edges that make up the path. When we relate
them to their corresponding geometry columns, we get the spatial
representation of the shortest paths. Instead of doing this for all
possible combinations of nodes, we can specify from and to which nodes
we want to calculate the shortest paths. Here, we will show an example
of a shortest path from one node to another, but it is just as well
possible to do the same for one to many, many to one, or many to many
nodes. Whenever the graph is weighted, the Dijkstra algoritm will be
used under the hood. Note here that we have to define the desired output
beforehand: `vpath` means that only the nodes (called vertices in
igraph) are returned, `epath` means that only the edges are returned,
and `both` returns them both.

``` r
path = st_network_paths(graph, from = 18, to = 1823, weights = "length")

path %>% 
  pull(node_paths) %>%
  unlist()
```

    ##  [1]   18   96   95   94  149  170  171 1491    3    4    5 2446  925    1 1713
    ## [16] 1623 2691 1828 1818 2450 1823

``` r
path %>% 
  pull(edge_paths) %>%
  unlist()
```

    ##  [1]   46  243  239  240  378  430  431    8    7   10   13 1940    2    4 3016
    ## [16] 3018 3232 3219 3220 3229

A subgraph of the original graph can now be created, containing only
those nodes and edges that are part of the calculated path. Note here
that something inconvenient happens: `igraph` seems to have it’s own
underlying structure for attaching indices to nodes, which makes that
the `from` and `to` columns in the egdes data don’t match with our
`nodeID` column in the nodes data anymore, but instead simply refer to
the row numbers of the nodes data.

``` r
path_graph <- graph %>%
    convert(to_spatial_shortest_paths, 18, 1823, weights = "length")

path_graph
```

    ## # A sfnetwork with 21 nodes and 20 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An unrooted tree with spatially explicit edges
    ## #
    ## # Edge Data:     20 x 7 (active)
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.620014 ymin: 51.9548 xmax: 7.63235 ymax: 51.96061
    ##    from    to highway                  geometry   length betweenness
    ##   <int> <int> <chr>            <LINESTRING [°]>      [m]       <dbl>
    ## 1     1    12 reside~ (7.624559 51.95567, 7.62~  5.8968~       62916
    ## 2     1    15 reside~ (7.624554 51.95561, 7.62~  4.5289~       50040
    ## 3     2     3 second~ (7.626498 51.95617, 7.62~ 29.0521~       84946
    ## 4     2    13 <NA>    (7.626552 51.95648, 7.62~ 33.7254~      119677
    ## 5     3     4 second~ (7.626103 51.9561, 7.626~ 20.7254~       78046
    ## 6     4    19 path    (7.625926 51.95595, 7.62~ 28.7160~       44828
    ## # ... with 14 more rows, and 1 more variable: .tidygraph_edge_index <int>
    ## #
    ## # Node Data:     21 x 4
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.620014 ymin: 51.9548 xmax: 7.63235 ymax: 51.96061
    ##              geometry degree betweenness .tidygraph_node_index
    ##           <POINT [°]>  <dbl>       <dbl>                 <int>
    ## 1 (7.624554 51.95561)      4       62444                     1
    ## 2 (7.626498 51.95617)      3      142779                     3
    ## 3  (7.626103 51.9561)      3      102529                     4
    ## # ... with 18 more rows

This subgraph can now be analyzed, for example by calculating the total
length of the path, …

``` r
path_graph %>%
  activate(edges) %>%
  as_tibble() %>%
  summarise(length = sum(length))
```

    ## although coordinates are longitude/latitude, st_union assumes that they are planar

    ## Simple feature collection with 1 feature and 1 field
    ## geometry type:  MULTILINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.620014 ymin: 51.9548 xmax: 7.63235 ymax: 51.96061
    ## geographic CRS: WGS 84
    ## # A tibble: 1 x 2
    ##     length                                                              geometry
    ##        [m]                                                 <MULTILINESTRING [°]>
    ## 1 1218.149 ((7.624559 51.95567, 7.624554 51.95561), (7.624554 51.95561, 7.62455~

… and plotted.

``` r
ggplot() +
  geom_sf(data = st_as_sf(graph, "edges"), col = 'darkgrey') +
  geom_sf(data = st_as_sf(graph, "nodes"), col = 'darkgrey', size = 0.5) +
  geom_sf(data = st_as_sf(path_graph, "edges"), lwd = 1, col = 'firebrick') +
  geom_sf(data = st_as_sf(path_graph, "nodes"), size = 2)
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

However, often we will be interested in shortest paths between
geographical points that are not necessarily nodes in the network. For
example, we might want to calculate the shortest path from the railway
station of Münster to the cathedral.

``` r
muenster_station <- st_point(c(7.6349, 51.9566)) %>% 
  st_sfc(crs = 4326)

muenster_cathedral <- st_point(c(7.626, 51.962)) %>%
  st_sfc(crs = 4326)

ggplot() +
  geom_sf(data = st_as_sf(graph, "edges"), col = 'darkgrey') +
  geom_sf(data = st_as_sf(graph, "nodes"), col = 'darkgrey', size = 0.5) +
  geom_sf(data = muenster_station, size = 2, col = 'firebrick') +
  geom_sf(data = muenster_cathedral, size = 2, col = 'firebrick') +
  geom_sf_label(data = muenster_station, aes(label = 'station'), nudge_x = 0.004) +
  geom_sf_label(data = muenster_cathedral, aes(label = 'cathedral'), nudge_x = 0.005)
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

To find the route on the network, we must first identify the nearest
points on the network. The `nabor` package has a well performing
function to do so. It does, however, require the coordinates of the
origin and destination nodes to be given in a matrix.

Like before, we use the ID to calculate the shortest path, and plot it:

``` r
path_graph <- graph %>%
    convert(
      to_spatial_shortest_paths, 
      from = muenster_station, 
      to = muenster_cathedral,
      weights = "length"
    )
```

    ## although coordinates are longitude/latitude, st_nearest_points assumes that they are planar
    ## although coordinates are longitude/latitude, st_nearest_points assumes that they are planar

``` r
path_graph
```

    ## # A sfnetwork with 52 nodes and 51 edges
    ## #
    ## # CRS:  EPSG:4326 
    ## #
    ## # An unrooted tree with spatially explicit edges
    ## #
    ## # Edge Data:     51 x 7 (active)
    ## # Geometry type: LINESTRING
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.625887 ymin: 51.95635 xmax: 7.634871 ymax: 51.96205
    ##    from    to highway                  geometry   length betweenness
    ##   <int> <int> <chr>            <LINESTRING [°]>      [m]       <dbl>
    ## 1     1     2 pedest~ (7.634264 51.95679, 7.63~ 10.3675~      197702
    ## 2     1    28 pedest~ (7.634264 51.95679, 7.63~  3.8174~      204688
    ## 3     2     3 pedest~ (7.634123 51.95682, 7.63~  9.3208~      201209
    ## 4     3    39 <NA>    (7.633995 51.95685, 7.63~  2.9228~      205349
    ## 5     4     5 unclas~ (7.631481 51.95742, 7.63~  9.0307~      285308
    ## 6     4    48 <NA>    (7.631481 51.95742, 7.63~ 74.0336~      177395
    ## # ... with 45 more rows, and 1 more variable: .tidygraph_edge_index <int>
    ## #
    ## # Node Data:     52 x 4
    ## # Geometry type: POINT
    ## # Dimension:     XY
    ## # Bounding box:  xmin: 7.625949 ymin: 51.95635 xmax: 7.634871 ymax: 51.96205
    ##              geometry degree betweenness .tidygraph_node_index
    ##           <POINT [°]>  <dbl>       <dbl>                 <int>
    ## 1 (7.634264 51.95679)      4      208480                    77
    ## 2 (7.634123 51.95682)      4      204467                    78
    ## 3 (7.633995 51.95685)      4      212715                    79
    ## # ... with 49 more rows

``` r
ggplot() +
  geom_sf(data = st_as_sf(graph, "edges"), col = 'darkgrey') +
  geom_sf(data = st_as_sf(graph, "nodes"), 
          col = 'darkgrey', size = 0.5) +
  geom_sf(data = st_as_sf(path_graph, "edges"), 
          lwd = 1, col = 'firebrick') +
  geom_sf(data = muenster_station, size = 2) +
  geom_sf(data = muenster_cathedral, size = 2)  +
  geom_sf_label(data = muenster_station, 
                aes(label = 'station'), nudge_x = 0.004) +
  geom_sf_label(data = muenster_cathedral, 
                aes(label = 'cathedral'), nudge_x = 0.005)
```

![](blogpost_2_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

It worked! We calculated a path from the rail station to the centre, a
common trip taken by tourists visiting Muenster. Clearly this is not a
finished piece of work but the post has demonstrated what is possible.
Future functionality should look to make spatial networks more user
friendly, including provision of ‘weighting profiles’, batch routing and
functions that reduce the number of steps needed to work with spatial
network data in R.

For alternative approaches and further reading, the following resources
are recommended:

-   sfnetworks, a GitHub package that implements some of the ideas in
    this post: <https://github.com/luukvdmeer/sfnetworks>
-   stplanr, a package for transport planning:
    <https://github.com/ropensci/stplanr>
-   dodgr, distances on directed graphs:
    <https://github.com/ATFutures/dodgr>
-   cppRouting, a package for routing in C++:
    <https://github.com/vlarmet/cppRouting>
-   Chapter 10 of Geocomputation with R, which provides context and
    demonstrates a transport planning workflow in R:
    <https://geocompr.robinlovelace.net/transport.html>