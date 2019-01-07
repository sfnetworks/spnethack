
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spnethack

The goal of spnethack is to explore spatial networks in R.

# Requirements

We’ll use the following packages:

``` r
library(sf)
library(osmdata)
library(dodgr)
library(stplanr)
```

# Data

It makes sense to have some input data.

## From josm

## From geojson

``` r
promenade = opq(bbox = "munster") %>% 
  add_osm_feature(key = "name", value = "Promenade")
```
