
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spnethack

The goal of spnethack is to explore spatial networks in R.

# Requirements

We’ll use the following packages:

``` r
library(sf)
library(osmdata)
library(dodgr)
#> Warning: package 'dodgr' was built under R version 3.5.2
library(stplanr)
library(dplyr)
library(piggyback)
#> Warning: package 'piggyback' was built under R version 3.5.2
```

# Data

It makes sense to have some input data.

## From josm

Data from OSM was downloaded with the `josm` GUI. It can be read-in as
follows:

``` r
pb_download("promenade-all.geojson")
promenade_all = read_sf("promenade-all.geojson")
summary(factor(promenade_all$highway))
promenade_min = promenade_all %>% 
  filter(name == "Promenade")
promenade_way = promenade_all %>% 
  filter(!is.na(highway))
# write_sf(promenade_way, "promenade-way.geojson")
# write_sf(promenade_min, "promenade-min.geojson")
```

The minimum dataset can be read-in as follows:

``` r
pb_download("promenade-min.geojson")
#> Warning in get_token(): Using default public GITHUB_TOKEN.
#>                      Please set your own token
#> <U+2714> Setting active project to 'C:/Users/Lore/Documents/spnethack'
#> All files up-to-date already
promenade_min = read_sf("promenade-min.geojson")
plot(promenade_min$geometry)
```

![](README_files/figure-gfm/plot1-1.png)<!-- -->

A slightly larger dataset can be read-in and plotted as follows:

``` r
pb_download("promenade-way.geojson")
#> Warning in get_token(): Using default public GITHUB_TOKEN.
#>                      Please set your own token
#> All files up-to-date already
promenade_way = geojsonsf::geojson_sf("promenade-way.geojson")
plot(promenade_way$geometry)
```

![](README_files/figure-gfm/pway-1.png)<!-- -->

## From osmdata

``` r
promenade_osmdata = opq(bbox = 'Muenster, DE') %>% 
  add_osm_feature(key = 'name', value = 'Promenade') %>% 
  osmdata_sf %>% 
  unique_osmdata()
```

## From dodgr

``` r
muenster = dodgr_streetnet('Muenster, DE')
promenade_dodgr = muenster %>% filter(name == 'Promenade')
```

## Route networks with stplanr

The code to create route network data is in `stplanr-promenade.R`. It
can generate a ‘betweenness’ graph like this:

``` r
source(file = "stplanr-promenade.R")
```

![](README_files/figure-gfm/promenade-stplanr-1.png)<!-- -->

## Route networks with dodgr

An example of how to create route data from sample points along the
edges can be found `dodgr-promenade.R`, it can also generate a flow
aggregate which is quite similar to the betweenness.

``` r
source(file = "dodgr-promenade.R")

dodgr_flowmap(graph_f, linescale = 5)
```

To do this we sample 100 000 points on the `to` column and the same
number on the `from`, and it takes 0.6 sec.

## Route networks with sfnetworks

## Route networks with spnetwork

## Route networks with osmnx
