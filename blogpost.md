Spatial networks in R with sf and tidygraph
================

## Introduction

Street networks, shipping routes, telecommunication lines, river
bassins. All examples of spatial networks: organized systems of nodes
and edges embedded in space. For most of them, these nodes and edges can
be associated with geographical coordinates. That is, the nodes are
geographical points, and the edges geographical lines.

Such spatial networks can be analyzed using graph theory. Not for
nothing, Leonhard Eulers famous work on the [Seven Bridges of
Köningsberg](https://www.mathsisfun.com/activity/seven-bridges-konigsberg.html),
which laid the foundations of graph theory and network analysis, was in
essence a spatial problem.

In R, there are advanced, modern tools for both the analysis of spatial
data and networks. Furthermore, several packages have been developed
that cover (parts of) spatial network analysis. As stated in this
[github issue](https://github.com/r-spatial/sf/issues/966) and this
[tweet](https://twitter.com/zevross/status/1089908839816794118),
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
[sf](https://github.com/r-spatial/sf) (see
[section 1.5](https://geocompr.robinlovelace.net/intro.html#the-history-of-r-spatial)
of Geocomputation with R for a brief history). Likewise, a number of
packages for graph representation and analysis have been developed,
notably [tidygraph](https://github.com/thomasp85/tidygraph), which is
based on [igraph](https://igraph.org/).

Both sf and tidygraph support the the `tibble` class and the broader
‘tidy’ approach to data science, which involves data processing
pipelines, type stability and a convention of representing everything as
a data frame (well a `tibble`, which is a data frame with user friendly
default settings). In sf, this means storing spatial vector data as
objects of class `sf`, which are essentially the same as a regular data
frame (or tibble), but with an additional ‘sticky’ list column
containing a geometry for each feature (row), and attributes such as
bounding box and CRS. Tidygraph stores networks in objects of class
`tbl_graph`. A `tbl_graph` is an `igraph` object, but enables the user
to manipulate both the edges and nodes elements as if they were data
frames also.

Both sf and tidygraph are relatively new packages (first released on
CRAN in 2016 and 2017, respectively). It is unsurprising, therefore,
that they have yet to be combined to allow a hybrid, tibble-based
representation of spatial networks.

Nevertheless, a number of approaches have been developed for
representing spatial networks, and some of these are in packages that
have been published on CRAN.
[stplanr](https://github.com/ropensci/stplanr), for instance, contains
the `SpatialLinesNetwork` class, which works with both the
[sp](https://github.com/edzer/sp/) (a package for spatial data analysis
launched in 2005) and sf packages.
[dodgr](https://github.com/ATFutures/dodgr) is a more recent package
that provides analytical tools for street networks, with a focus on
directed graphs (that can have direction-dependent weights,
e.g. representing a one-way street). Other packages seeking to
implement spatial networks in R include
[spnetwork](https://github.com/edzer/spnetwork), a package that defined
a class system combining sp and igraph, and
[shp2graph](https://cran.r-project.org/web/packages/shp2graph/index.html),
which provides tools to switch between sp and igraph objects.

Each package has its merits that deserve to be explored in more detail
(possibly in a subsequent blog post). The remainder of this post
outlines an approach that combines `sf` and `igraph` objects in a
`tidygraph` object.

## Set-up

The following code chunk will install the packages used in this post:

``` r
# we'll use remotes to install packages, install it if needs be:
if(!"remotes" %in% installed.packages()) {
  install.packages("remotes")
}

cran_pkgs = c(
  "sf",
  "tidygraph",
  "igraph",
  "osmdata",
  "dodgr",
  "stplanr",
  "dplyr",
  "tibble",
  "ggplot2",
  "units",
  "piggyback",
  "tmap",
  "rgrass7",
  "link2GI"
)

remotes::install_cran(cran_pkgs)
```

``` r
library(sf)
library(tidygraph)
library(igraph)
library(dplyr)
library(tibble)
library(ggplot2)
library(units)
```

## Getting the data

As an example, we use the street network of the city center of Münster,
Germany. We will get the data from OpenStreetMap. Packages like `dodgr`
have optimized their code for such data, however considering that we
want to showcase this workflow for any source of data, we will generate
an object of class `sf` containing only `LINESTRING` geometries. One
additional variable, the type of street, is added to show that the same
steps can be used for `sf` objects that contain any number of additional
variables.

``` r
library(osmdata)
```

    ## Data (c) OpenStreetMap contributors, ODbL 1.0. http://www.openstreetmap.org/copyright

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

    ## Simple feature collection with 2191 features and 1 field
    ## geometry type:  LINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.601942 ymin: 51.94823 xmax: 7.645597 ymax: 51.97241
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## First 10 features:
    ##             highway                       geometry
    ## 4064462     primary LINESTRING (7.624554 51.955...
    ## 4064463     primary LINESTRING (7.626498 51.956...
    ## 4064467 residential LINESTRING (7.630898 51.955...
    ## 4064474     primary LINESTRING (7.61972 51.9554...
    ## 4064476     primary LINESTRING (7.619844 51.954...
    ## 4064482    tertiary LINESTRING (7.616395 51.957...
    ## 4064485     service LINESTRING (7.63275 51.9603...
    ## 4984982   secondary LINESTRING (7.614156 51.967...
    ## 4985138    cycleway LINESTRING (7.61525 51.9673...
    ## 4985140 residential LINESTRING (7.616774 51.968...

``` r
ggplot(data = muenster_center) + geom_sf()
```

![](blogpost_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## From sf to tbl\_graph: a step wise approach

### Step 1: Clean the network

To perform network analysis, we need a network with a clean topology. In
theory, the best way to clean up the network topology is by manual
editing, but this can be very labour intensive and time consuming,
mainly for large networks. The
[v.clean](https://grass.osgeo.org/grass77/manuals/v.clean.html) toolset
from the GRASS GIS software provides automated functionalities for this
task, and is therefore a popular instrument within the field of spatial
network analysis. As far as we know, there is no R equivalent for this
toolset, but fortunately, the
[rgrass7](https://cran.r-project.org/web/packages/rgrass7/index.html)
and [link2GI](https://github.com/r-spatial/link2GI) packages enable us
to easily ‘bridge’ to GRASS GIS. Obviously, this requires to have GRASS
GIS installed on your computer. For an in depth description of combining
R with open source GIS software, see
[Chapter 9](https://geocompr.robinlovelace.net/gis.html) of
Geocomputation with R. Take into account that the linking process may
take up some time, especially on windows operating systems.

Here, we will clean the network topology by breaking lines at
intersections and also breaking lines that form a collapsed loop. This
will be followed by a removal of duplicated geometry features. Once
done, we will read the data back into R, and convert again into an `sf`
object with `LINESTRING` geometry.

``` r
# Link tot GRASS GIS
link2GI::linkGRASS7(muenster_center, ver_select = TRUE)
```

    ## Registered S3 method overwritten by 'xts':
    ##   method     from
    ##   as.zoo.xts zoo

    ## Registered S3 method overwritten by 'R.oo':
    ##   method        from       
    ##   throw.default R.methodsS3

``` r
library(rgrass7)

# Add data to GRASS spatial database
writeVECT(
  SDF = as(muenster_center, 'Spatial'), 
  vname = 'muenster_center', 
  v.in.ogr_flags = 'overwrite'
)

# Execute the v.clean tool
execGRASS(
  cmd = 'v.clean', 
  input = 'muenster_center', 
  output = 'muenster_cleaned',        
  tool = 'break', 
  flags = c('overwrite', 'c')
)

# Read back into R
muenster_center <- readVECT('muenster_cleaned') %>% 
  st_as_sf() %>%
  st_cast('LINESTRING') %>%
  select(-cat)
```

``` r
muenster_center
```

    ## Simple feature collection with 4673 features and 1 field
    ## geometry type:  LINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.601942 ymin: 51.94823 xmax: 7.645597 ymax: 51.97241
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## First 10 features:
    ##         highway                       geometry
    ## 1       primary LINESTRING (7.625074 51.955...
    ## 1.1     primary LINESTRING (7.624554 51.955...
    ## 2       primary LINESTRING (7.626498 51.956...
    ## 2.1     primary LINESTRING (7.625926 51.955...
    ## 2.2     primary LINESTRING (7.626103 51.956...
    ## 3   residential LINESTRING (7.630898 51.955...
    ## 3.1 residential LINESTRING (7.63084 51.9556...
    ## 3.2 residential LINESTRING (7.63048 51.9552...
    ## 3.3 residential LINESTRING (7.630307 51.955...
    ## 4       primary LINESTRING (7.61972 51.9554...

### Step 2: Give each edge an unique index

The edges of the network, are simply the linestrings in the data. Each
of them gets an unique index, which can be later be related to their
start and end node.

``` r
edges <- muenster_center %>%
  mutate(edgeID = c(1:n()))

edges
```

    ## Simple feature collection with 4673 features and 2 fields
    ## geometry type:  LINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.601942 ymin: 51.94823 xmax: 7.645597 ymax: 51.97241
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## First 10 features:
    ##        highway                       geometry edgeID
    ## 1      primary LINESTRING (7.625074 51.955...      1
    ## 2      primary LINESTRING (7.624554 51.955...      2
    ## 3      primary LINESTRING (7.626498 51.956...      3
    ## 4      primary LINESTRING (7.625926 51.955...      4
    ## 5      primary LINESTRING (7.626103 51.956...      5
    ## 6  residential LINESTRING (7.630898 51.955...      6
    ## 7  residential LINESTRING (7.63084 51.9556...      7
    ## 8  residential LINESTRING (7.63048 51.9552...      8
    ## 9  residential LINESTRING (7.630307 51.955...      9
    ## 10     primary LINESTRING (7.61972 51.9554...     10

### Step 3: Create nodes at the start and end point of each edge

The nodes of the network, are the start and end points of the edges. The
locations of these points can be derived by using the `st_coordinates`
function in sf. When given a set of linestrings, this function breaks
down each of them into the points they are build up. It returns a matrix
with the X and Y coordinates of those points, and additionally an
integer indicator L1 specifying to which line a point belongs. These
integer indicators correspond to the edge indices defined in step 1.
That is, if we convert the matrix into a `data.frame` or `tibble`, group
the features by the edge index, and only keep the first and last feature
of each group, we have the start and end points of the linestrings.

``` r
nodes <- edges %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(edgeID = L1) %>%
  group_by(edgeID) %>%
  slice(c(1, n())) %>%
  ungroup() %>%
  mutate(start_end = rep(c('start', 'end'), times = n()/2))

nodes
```

    ## # A tibble: 9,346 x 4
    ##        X     Y edgeID start_end
    ##    <dbl> <dbl>  <dbl> <chr>    
    ##  1  7.63  52.0      1 start    
    ##  2  7.63  52.0      1 end      
    ##  3  7.62  52.0      2 start    
    ##  4  7.63  52.0      2 end      
    ##  5  7.63  52.0      3 start    
    ##  6  7.63  52.0      3 end      
    ##  7  7.63  52.0      4 start    
    ##  8  7.63  52.0      4 end      
    ##  9  7.63  52.0      5 start    
    ## 10  7.63  52.0      5 end      
    ## # … with 9,336 more rows

### Step 4: Give each node an unique index

Each of the nodes in the network needs to get an unique index, such that
they can be related to the edges. However, we need to take into account
that edges can share either startpoints and/or endpoints. Such
duplicated points, that have the same X and Y coordinate, are one single
node, and should therefore get the same index. Note that the coordinate
values as displayed in the tibble are rounded, and may look the same for
several rows, even when they are not. We can use the `group_indices`
function in dplyr to give each group of unique X,Y-combinations an
unique index.

``` r
nodes <- nodes %>%
  mutate(xy = paste(.$X, .$Y)) %>% 
  mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
  select(-xy)

nodes
```

    ## # A tibble: 9,346 x 5
    ##        X     Y edgeID start_end nodeID
    ##    <dbl> <dbl>  <dbl> <chr>      <int>
    ##  1  7.63  52.0      1 start          1
    ##  2  7.63  52.0      1 end            2
    ##  3  7.62  52.0      2 start          3
    ##  4  7.63  52.0      2 end            1
    ##  5  7.63  52.0      3 start          4
    ##  6  7.63  52.0      3 end            5
    ##  7  7.63  52.0      4 start          6
    ##  8  7.63  52.0      4 end            7
    ##  9  7.63  52.0      5 start          5
    ## 10  7.63  52.0      5 end            6
    ## # … with 9,336 more rows

### Step 5: Combine the node indices with the edges

Now each of the start and endpoints from step has been assigned a node
ID in step 4, we can add the node indices to the edges. In other words,
we can specify for each edge, in which node it starts, and in which node
it ends.

``` r
source_nodes <- nodes %>%
  filter(start_end == 'start') %>%
  pull(nodeID)

target_nodes <- nodes %>%
  filter(start_end == 'end') %>%
  pull(nodeID)

edges = edges %>%
  mutate(from = source_nodes, to = target_nodes)

edges
```

    ## Simple feature collection with 4673 features and 4 fields
    ## geometry type:  LINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.601942 ymin: 51.94823 xmax: 7.645597 ymax: 51.97241
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## First 10 features:
    ##        highway                       geometry edgeID from to
    ## 1      primary LINESTRING (7.625074 51.955...      1    1  2
    ## 2      primary LINESTRING (7.624554 51.955...      2    3  1
    ## 3      primary LINESTRING (7.626498 51.956...      3    4  5
    ## 4      primary LINESTRING (7.625926 51.955...      4    6  7
    ## 5      primary LINESTRING (7.626103 51.956...      5    5  6
    ## 6  residential LINESTRING (7.630898 51.955...      6    8  9
    ## 7  residential LINESTRING (7.63084 51.9556...      7    9 10
    ## 8  residential LINESTRING (7.63048 51.9552...      8   10 11
    ## 9  residential LINESTRING (7.630307 51.955...      9   11 12
    ## 10     primary LINESTRING (7.61972 51.9554...     10   13 14

### Step 6: Remove duplicate nodes

Having added the unique node ID’s to the edges data, we don’t need the
duplicated start and endpoints anymore. After removing them, we end up
with a `tibble` in which each row represents a unique, single node. This
tibble can be converted into an `sf` object, with `POINT` geometries.

``` r
nodes <- nodes %>%
  distinct(nodeID, .keep_all = TRUE) %>%
  select(-c(edgeID, start_end)) %>%
  st_as_sf(coords = c('X', 'Y')) %>%
  st_set_crs(st_crs(edges))

nodes
```

    ## Simple feature collection with 3322 features and 1 field
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: 7.601942 ymin: 51.94823 xmax: 7.645597 ymax: 51.97241
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## # A tibble: 3,322 x 2
    ##    nodeID            geometry
    ##     <int>         <POINT [°]>
    ##  1      1 (7.625074 51.95564)
    ##  2      2 (7.625275 51.95566)
    ##  3      3 (7.624554 51.95561)
    ##  4      4 (7.626498 51.95617)
    ##  5      5  (7.626103 51.9561)
    ##  6      6 (7.625926 51.95595)
    ##  7      7 (7.625906 51.95573)
    ##  8      8 (7.630898 51.95573)
    ##  9      9  (7.63084 51.95567)
    ## 10     10  (7.63048 51.95524)
    ## # … with 3,312 more rows

### Step 7: Convert to tbl\_graph

The first six steps led to one `sf` object with `LINESTRING` geometries,
representing the edges of the network, and one `sf` object with `POINT`
geometries, representing the nodes of the network. The `tbl_graph`
function allows us to convert these two into a `tbl_graph` object. There
are two tricky parts in this step that need to be highlighted. One, is
that the columns containing the indices of the source and target nodes
should either be the first two columns of the `sf` object, or be named
‘to’ and ‘from’, respectively. Secondly, inside the `tbl_graph`
function, these columns are converted into a two-column matrix. However,
an `sf` object has a so-called ‘sticky geometry’, which means that the
geometry column sticks to the attributes whenever specific columns are
selected. Therefore, the matrix created inside `tbl_graph` has three
columns instead of two, and that causes an error. Therefore, we first
need to convert the `sf` object to a regular `data.frame` or `tibble`,
before we can construct a `tbl_graph`. In the end, this doesn’t matter,
since both the nodes and edges will be ‘integrated’ into an `igraph`
structure, and loose their specific `sf`
characteristics.

``` r
graph = tbl_graph(nodes = nodes, edges = as_tibble(edges), directed = FALSE)

graph
```

    ## # A tbl_graph: 3322 nodes and 4673 edges
    ## #
    ## # An undirected multigraph with 37 components
    ## #
    ## # Node Data: 3,322 x 2 (active)
    ##   nodeID            geometry
    ##    <int>         <POINT [°]>
    ## 1      1 (7.625074 51.95564)
    ## 2      2 (7.625275 51.95566)
    ## 3      3 (7.624554 51.95561)
    ## 4      4 (7.626498 51.95617)
    ## 5      5  (7.626103 51.9561)
    ## 6      6 (7.625926 51.95595)
    ## # … with 3,316 more rows
    ## #
    ## # Edge Data: 4,673 x 5
    ##    from    to highway                                       geometry edgeID
    ##   <int> <int> <fct>                                 <LINESTRING [°]>  <int>
    ## 1     1     2 primary         (7.625074 51.95564, 7.625275 51.95566)      1
    ## 2     1     3 primary (7.624554 51.95561, 7.624838 51.95562, 7.6249…      2
    ## 3     4     5 primary (7.626498 51.95617, 7.626437 51.95617, 7.6263…      3
    ## # … with 4,670 more rows

### Step 8: Putting it together

For make the approach more convenient, we can combine all steps above
into a single function, that takes a cleaned `sf` object with
`LINESTRING` geometries as input, and returns a spatial `tbl_graph`.

``` r
sf_to_tidygraph = function(x, directed = TRUE) {
  
  edges <- x %>%
    mutate(edgeID = c(1:n()))
  
  nodes <- edges %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(edgeID = L1) %>%
    group_by(edgeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
    mutate(xy = paste(.$X, .$Y)) %>% 
    mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
    select(-xy)
  
  source_nodes <- nodes %>%
    filter(start_end == 'start') %>%
    pull(nodeID)

  target_nodes <- nodes %>%
    filter(start_end == 'end') %>%
    pull(nodeID)

  edges = edges %>%
    mutate(from = source_nodes, to = target_nodes)
  
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(edgeID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(edges))
  
  tbl_graph(nodes = nodes, edges = as_tibble(edges), directed = directed)
  
}

sf_to_tidygraph(muenster_center, directed = FALSE)
```

    ## # A tbl_graph: 3322 nodes and 4673 edges
    ## #
    ## # An undirected multigraph with 37 components
    ## #
    ## # Node Data: 3,322 x 2 (active)
    ##   nodeID            geometry
    ##    <int>         <POINT [°]>
    ## 1      1 (7.625074 51.95564)
    ## 2      2 (7.625275 51.95566)
    ## 3      3 (7.624554 51.95561)
    ## 4      4 (7.626498 51.95617)
    ## 5      5  (7.626103 51.9561)
    ## 6      6 (7.625926 51.95595)
    ## # … with 3,316 more rows
    ## #
    ## # Edge Data: 4,673 x 5
    ##    from    to highway                                       geometry edgeID
    ##   <int> <int> <fct>                                 <LINESTRING [°]>  <int>
    ## 1     1     2 primary         (7.625074 51.95564, 7.625275 51.95566)      1
    ## 2     1     3 primary (7.624554 51.95561, 7.624838 51.95562, 7.6249…      2
    ## 3     4     5 primary (7.626498 51.95617, 7.626437 51.95617, 7.6263…      3
    ## # … with 4,670 more rows

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
  mutate(length = st_length(geometry))

graph
```

    ## # A tbl_graph: 3322 nodes and 4673 edges
    ## #
    ## # An undirected multigraph with 37 components
    ## #
    ## # Edge Data: 4,673 x 6 (active)
    ##    from    to highway                              geometry edgeID   length
    ##   <int> <int> <fct>                        <LINESTRING [°]>  <int>      [m]
    ## 1     1     2 primary  (7.625074 51.95564, 7.625275 51.955…      1 14.0982…
    ## 2     1     3 primary  (7.624554 51.95561, 7.624838 51.955…      2 35.8537…
    ## 3     4     5 primary  (7.626498 51.95617, 7.626437 51.956…      3 29.0521…
    ## 4     6     7 primary  (7.625926 51.95595, 7.625907 51.955…      4 24.8519…
    ## 5     5     6 primary  (7.626103 51.9561, 7.626058 51.9560…      5 20.7254…
    ## 6     8     9 residen… (7.630898 51.95573, 7.63084 51.9556…      6  8.0905…
    ## # … with 4,667 more rows
    ## #
    ## # Node Data: 3,322 x 2
    ##   nodeID            geometry
    ##    <int>         <POINT [°]>
    ## 1      1 (7.625074 51.95564)
    ## 2      2 (7.625275 51.95566)
    ## 3      3 (7.624554 51.95561)
    ## # … with 3,319 more rows

With one flow of pipes, we can ‘escape’ the graph structure, turn either
the edges or nodes back into real `sf` objects, and, for example,
summarise the data based on a specific variable.

``` r
graph %>%
  activate(edges) %>%
  as_tibble() %>%
  st_as_sf() %>%
  group_by(highway) %>%
  summarise(length = sum(length))
```

    ## Warning: Factor `highway` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    ## Simple feature collection with 17 features and 2 fields
    ## geometry type:  MULTILINESTRING
    ## dimension:      XY
    ## bbox:           xmin: 7.601942 ymin: 51.94823 xmax: 7.645597 ymax: 51.97241
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## # A tibble: 17 x 3
    ##    highway         length                                          geometry
    ##  * <fct>              [m]                             <MULTILINESTRING [°]>
    ##  1 corridor        9.377… ((7.620506 51.96262, 7.620542 51.96266), (7.6205…
    ##  2 cycleway    22929.349… ((7.619683 51.95395, 7.619641 51.95378, 7.619559…
    ##  3 footway     42782.207… ((7.640529 51.95325, 7.640528 51.95323), (7.6405…
    ##  4 path         7642.403… ((7.624007 51.95379, 7.624223 51.95378, 7.624253…
    ##  5 pedestrian  11439.762… ((7.620362 51.95471, 7.620477 51.9547), (7.62012…
    ##  6 primary      3539.164… ((7.625556 51.95272, 7.625594 51.95284, 7.625714…
    ##  7 primary_li…   184.385… ((7.617285 51.96609, 7.617286 51.96624, 7.617295…
    ##  8 residential 22714.998… ((7.614509 51.95351, 7.614554 51.95346), (7.6326…
    ##  9 secondary    4471.930… ((7.631252 51.95402, 7.631405 51.95399), (7.6311…
    ## 10 secondary_…   160.708… ((7.635309 51.95946, 7.635705 51.95948), (7.6349…
    ## 11 service     26954.215… ((7.624803 51.95393, 7.625072 51.95393), (7.6158…
    ## 12 steps        1293.775… ((7.634423 51.9546, 7.634438 51.95462), (7.61430…
    ## 13 tertiary     4305.225… ((7.607112 51.94991, 7.607126 51.94992, 7.607183…
    ## 14 tertiary_l…    43.856… ((7.623592 51.96612, 7.623568 51.96611, 7.623468…
    ## 15 track         389.866… ((7.610671 51.95778, 7.610571 51.95759, 7.610585…
    ## 16 unclassifi…   610.488… ((7.634492 51.95613, 7.634689 51.95611), (7.6343…
    ## 17 <NA>         3162.396… ((7.634374 51.95579, 7.634545 51.95575, 7.634662…

Switching back to `sf` objects is useful as well when plotting the
network, in a way that preserves its spatial properties.

``` r
ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf()) + 
  geom_sf(data = graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), size = 0.5)
```

![](blogpost_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Or, alternatively, in only a few lines of code, plot the network as an
interactive
    map.

``` r
library(tmap)
```

    ## Warning in fun(libname, pkgname): rgeos: versions of GEOS runtime 3.7.1-CAPI-1.11.1
    ## and GEOS at installation 3.7.0-CAPI-1.11.0differ

``` r
tmap_mode('view')
```

    ## tmap mode set to interactive viewing

``` r
tm_shape(graph %>% activate(edges) %>% as_tibble() %>% st_as_sf()) +
  tm_lines() +
tm_shape(graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf()) +
  tm_dots() +
tmap_options(basemaps = 'OpenStreetMap')
```

![](blogpost_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

All nice and well, but these are not things that we necessarily need the
graph representation for. The added value of tidygraph, is that it opens
the door to the functions of the igraph library, all specifically
designed for network analysis, and enables us to use them inside a
‘tidy’ workflow. To cover them all, we would need to write a book,
but let’s at least show a few examples below.

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

    ## # A tbl_graph: 3322 nodes and 4673 edges
    ## #
    ## # An undirected multigraph with 37 components
    ## #
    ## # Edge Data: 4,673 x 7 (active)
    ##    from    to highway                  geometry edgeID   length betweenness
    ##   <int> <int> <fct>            <LINESTRING [°]>  <int>      [m]       <dbl>
    ## 1     1     2 primary (7.625074 51.95564, 7.62…      1 14.0982…       47471
    ## 2     1     3 primary (7.624554 51.95561, 7.62…      2 35.8537…       54910
    ## 3     4     5 primary (7.626498 51.95617, 7.62…      3 29.0521…      112002
    ## 4     6     7 primary (7.625926 51.95595, 7.62…      4 24.8519…       19822
    ## 5     5     6 primary (7.626103 51.9561, 7.626…      5 20.7254…       98127
    ## 6     8     9 reside… (7.630898 51.95573, 7.63…      6  8.0905…       26508
    ## # … with 4,667 more rows
    ## #
    ## # Node Data: 3,322 x 4
    ##   nodeID            geometry degree betweenness
    ##    <int>         <POINT [°]>  <dbl>       <dbl>
    ## 1      1 (7.625074 51.95564)      3       54368
    ## 2      2 (7.625275 51.95566)      2       45630
    ## 3      3 (7.624554 51.95561)      4      130103
    ## # … with 3,319 more rows

``` r
ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf()) + 
  geom_sf(data = graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), aes(col = betweenness), size = 2) +
  scale_colour_viridis_c(option = 'inferno')
```

![](blogpost_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), aes(col = betweenness), lwd = 1) +
  scale_colour_viridis_c(option = 'inferno')
```

![](blogpost_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

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
distances <- distances(
  graph = graph,
  weights = graph %>% activate(edges) %>% pull(length)
)

distances[1:5, 1:5]
```

    ##           [,1]      [,2]      [,3]      [,4]      [,5]
    ## [1,]   0.00000  14.09821  35.85371 126.03936  96.98722
    ## [2,]  14.09821   0.00000  49.95192 111.94116  82.88901
    ## [3,]  35.85371  49.95192   0.00000 159.67782 130.62567
    ## [4,] 126.03936 111.94116 159.67782   0.00000  29.05214
    ## [5,]  96.98722  82.88901 130.62567  29.05214   0.00000

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
from_node <- graph %>%
  activate(nodes) %>%
  filter(nodeID == 2404) %>%
  pull(nodeID)

to_node <- graph %>%
  activate(nodes) %>%
  filter(nodeID == 1679) %>%
  pull(nodeID)

path <- shortest_paths(
  graph = graph,
  from = from_node,
  to = to_node,
  output = 'both',
  weights = graph %>% activate(edges) %>% pull(length)
)

path$vpath
```

    ## [[1]]
    ## + 34/3322 vertices, from 297924e:
    ##  [1] 2404 2403 1230 1229 1228 1367 1368 1403 2868 1187 1188  148 1220 1221
    ## [15] 1222 1424 1423 1430 1431 1445  458  459 1614  506  498  499  500  501
    ## [29]  502  503 1721 1722 1678 1679

``` r
path$epath
```

    ## [[1]]
    ## + 33/4673 edges from 297924e:
    ##  [1] 2403--2404 1230--2403 1229--1230 1228--1229 1228--1367 1367--1368
    ##  [7] 1368--1403 1403--2868 1187--2868 1187--1188  148--1188  148--1220
    ## [13] 1220--1221 1221--1222 1222--1424 1423--1424 1423--1430 1430--1431
    ## [19] 1431--1445  458--1445  458-- 459  459--1614  506--1614  498-- 506
    ## [25]  498-- 499  499-- 500  500-- 501  501-- 502  502-- 503  503--1721
    ## [31] 1721--1722 1678--1722 1678--1679

``` r
path_graph <- graph %>%
  activate(edges) %>%
  slice(path$epath %>% unlist()) %>%
  activate(nodes) %>%
  slice(path$vpath %>% unlist())

path_graph
```

    ## # A tbl_graph: 34 nodes and 33 edges
    ## #
    ## # An undirected simple graph with 1 component
    ## #
    ## # Node Data: 34 x 4 (active)
    ##   nodeID            geometry degree betweenness
    ##    <int>         <POINT [°]>  <dbl>       <dbl>
    ## 1    148 (7.612122 51.95986)      3      116092
    ## 2    458 (7.612711 51.96182)      3       66173
    ## 3    459 (7.612722 51.96202)      2       48801
    ## 4    498 (7.612767 51.96253)      4       27074
    ## 5    499 (7.612659 51.96259)      3       14700
    ## 6    500 (7.611899 51.96294)      4       38902
    ## # … with 28 more rows
    ## #
    ## # Edge Data: 33 x 7
    ##    from    to highway                  geometry edgeID   length betweenness
    ##   <int> <int> <fct>            <LINESTRING [°]>  <int>      [m]       <dbl>
    ## 1     2     3 service (7.612711 51.96182, 7.61…    399 21.8535…       51010
    ## 2     4     5 service (7.612767 51.96253, 7.61…    433  9.6982…       16712
    ## 3     6     7 service (7.611899 51.96294, 7.61…    434 31.0929…       34059
    ## # … with 30 more rows

``` r
ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), col = 'darkgrey') + 
  geom_sf(data = graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), size = 0.5, col = 'darkgrey') +
  geom_sf(data = path_graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), lwd = 1, col = 'magenta') + 
  geom_sf(data = path_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), size = 1)
```

![](blogpost_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

However, often we will be interested in shortest paths between
geographical points that are not necessarily nodes in the network. For
example, we might want to calculate the shortest path from the railway
station of Münster to the
cathedral.

``` r
muenster_station <- st_point(c(7.6349, 51.9566)) %>% st_sfc(crs = 4326) %>% st_coordinates()

muenster_station
```

    ##        X       Y
    ## 1 7.6349 51.9566

``` r
muenster_cathedral <- st_point(c(7.626, 51.962)) %>%
  st_sfc(crs = 4326) %>%
  st_coordinates()

muenster_cathedral
```

    ##       X      Y
    ## 1 7.626 51.962

``` r
coords_o = matrix(muenster_station, ncol = 2)
coords_d = matrix(muenster_cathedral, ncol = 2)
```

To find the route on the network, we must first identify the nearest
points on the network.

``` r
nodes_coords = st_coordinates(nodes)
node_index_o = nabor::knn(data = nodes_coords, query = coords_o, k = 1)
node_index_d = nabor::knn(data = nodes_coords, query = coords_d, k = 1)
node_o = nodes[node_index_o$nn.idx, ]
node_d = nodes[node_index_d$nn.idx, ]
```

Like before, we use the ID to calculate the shortest path:

``` r
path <- shortest_paths(
  graph = graph,
  from = node_o$nodeID, # new origin
  to = node_d$nodeID,   # new destination
  output = 'both',
  weights = graph %>% activate(edges) %>% pull(length)
)

path$vpath
```

    ## [[1]]
    ## + 58/3322 vertices, from 297924e:
    ##  [1] 3275 3085 3084  317 1604  844  309 1606   98   99  100 2206 2209 2207
    ## [15] 1449  487 1160 2208 2501  102  101 2579 1920  144  145 1384  119 1593
    ## [29] 1586  850  851  848  849 1589 2578 2791 1460 1478 1477  336  337 2067
    ## [43] 2066 2069 2070 2072 2071 2074 2073 2075 2068 2065  200  199  831  832
    ## [57]  833  400

``` r
path$epath
```

    ## [[1]]
    ## + 57/4673 edges from 297924e:
    ##  [1] 3085--3275 3084--3085  317--3084  317--1604  844--1604  309-- 844
    ##  [7]  309--1606   98--1606   98--  99   99-- 100  100--2206 2206--2209
    ## [13] 2207--2209 1449--2207  487--1449  487--1160 1160--2208 2208--2501
    ## [19]  102--2501  101-- 102  101--2579 1920--2579  144--1920  144-- 145
    ## [25]  145--1384  119--1384  119--1593 1586--1593  850--1586  850-- 851
    ## [31]  848-- 851  848-- 849  849--1589 1589--2578 2578--2791 1460--2791
    ## [37] 1460--1478 1477--1478  336--1477  336-- 337  337--2067 2066--2067
    ## [43] 2066--2069 2069--2070 2070--2072 2071--2072 2071--2074 2073--2074
    ## [49] 2073--2075 2068--2075 2065--2068  200--2065  199-- 200  199-- 831
    ## [55]  831-- 832  832-- 833  400-- 833

``` r
path_graph <- graph %>%
  activate(edges) %>%
  slice(path$epath %>% unlist()) %>%
  activate(nodes) %>%
  slice(path$vpath %>% unlist())

path_sf = path_graph %>% 
  activate(edges) %>% 
  as_tibble() %>%
  st_as_sf()
mapview::mapview(path_sf)
```

![](blogpost_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

It worked\! We calculated a path from the rail station to the centre, a
common trip taken by tourists visiting Muenster. Clearly this is not a
finished piece of work but the post has demonstrated what is possible.
Future functionality should look to make spatial networks more user
friendly, including provision of ‘weighting profiles’, batch routing and
functions that reduce the number of steps needed to work with spatial
network data in R.

For alternative approaches and further reading, the following resources
are recommended:

  - sfnetworks, a GitHub package that implements some of the ideas in
    this post: <https://github.com/luukvdmeer/sfnetworks>

  - stplanr, a package for transport planning:
    <https://github.com/ropensci/stplanr>

  - dodgr, distances on directed graphs:
    <https://github.com/ATFutures/dodgr>

  - cppRouting, a package for routing in C++:
    <https://github.com/vlarmet/cppRouting>

  - Chapter 10 of Geocomputation with R, which provides context and
    demonstrates a transport planning workflow in R:
    <https://geocompr.robinlovelace.net/transport.html>

  - 

<!-- For this purpose, the dodgr package, designed with a focus on routing in street networks, comes in handy. If we provide it the network, and the coordinates of the desired from and two points as two-column matrices, it will find the network points closest to the given locations, and calculate the shortest paths. In computing time, dodgr even outperforms igraph.  -->

<!-- However, dodgr does not work with sf-like geometry list columns (yet?). Instead, it assumes the graph to be a single data.frame, or similar object, containing the longitudes and latitudes of the from and to nodes in distinct columns. Also, the distance column should be a numeric column called `dist`, and not a `units` class called `length`. Therefore, we need to take some pre-processing steps before we can use dodgr. Again, we will show an example of a shortest path from one location to another, but it is just as well possible to do the same for one to many, many to one, or many to many locations. -->

<!-- ```{r, results='asis'} -->

<!-- library(dodgr) -->

<!-- node_coordinates <- graph %>% -->

<!--   activate(nodes) %>% -->

<!--   as_tibble() %>% -->

<!--   st_as_sf() %>% -->

<!--   st_coordinates() %>% -->

<!--   as_tibble() %>% -->

<!--   mutate(nodeID = graph %>% activate(nodes) %>% pull(nodeID)) -->

<!-- dodgr_graph <- graph %>% -->

<!--   activate(edges) %>% -->

<!--   as_tibble() %>% -->

<!--   mutate(dist = as.numeric(length)) %>%  -->

<!--   rename(edge_id = edgeID) %>% -->

<!--   left_join(node_coordinates, by = c('from' = 'nodeID')) %>% -->

<!--   rename(from_id = from, from_lon = X, from_lat = Y) %>% -->

<!--   left_join(node_coordinates, by = c('to' = 'nodeID')) %>% -->

<!--   rename(to_id = to, to_lon = X, to_lat = Y) %>% -->

<!--   select(edge_id, from_id, from_lon, from_lat, to_id, to_lon, to_lat, dist) %>%  -->

<!--   as.data.frame() %>% # a tibble is not recognized by dodgr -->

<!--   dodgr_components() # a requirement to for a dodgr_streetnet  -->

<!-- dodgr_graph %>% head(10) %>% knitr::kable() -->

<!-- ``` -->

<!-- ```{r, error = TRUE} -->

<!-- muenster_station <- st_point(c(7.6349, 51.9566)) %>% st_sfc(crs = 4326) %>% st_coordinates() -->

<!-- muenster_station -->

<!-- muenster_cathedral <- st_point(c(7.62199751, 51.957829502)) %>% st_sfc(crs = 4326) %>% st_coordinates() -->

<!-- muenster_cathedral -->

<!-- path <- dodgr::dodgr_paths( -->

<!--   graph = dodgr_graph, -->

<!--   from = muenster_cathedral, -->

<!--   to = muenster_station -->

<!-- ) -->

<!-- path ## It works but just cannot find a path between the two points. Possibly has to do with the fact that the graph is undirected and dodgr only works with directed graphs.  -->

<!-- ``` -->
