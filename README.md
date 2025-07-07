## *Grassland bird species decline with colonial-era landscape change in a tropical montane ecosystem*  

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) <!-- badges: end -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15831546.svg)](https://doi.org/10.5281/zenodo.15831546)  

This repository contains code and analysis for a manuscript that examines how a century of landscape changes relate to bird species populations.   

This manuscript is published in *Global Change Biology*. Please reach out to the repository owner if you have questions.  

![a) The Nilgiri hills of the Western Ghats-Sri Lanka biodiversity hotspot is the focus of our study. Shown here is the rare historical map of the Nilgiri hills, along with representative locations where land cover transitions occurred between historical (1848, colored in violet) and modern (2018, colored in red) time periods. Images of land cover were taken by Vijay Ramesh, Chandrasekar Das and Wikimedia Commons (captured by Jaseem Hamza). b) A timeline of key events that shaped the landscape of the Nilgiri hills alongside our data sources (colored in orange) is shown here.](figs/figure1.png)

### [Readable version](https://vjjan91.github.io/nilgiris-resurvey-project/)

A readable version of this analysis is available in bookdown format by clicking on the heading above.  

### Source code for the analyses

We describe what each script (`.Rmd`) of this repository is intended to achieve below. Please note that the first few scripts deal with processing climatic data. However, our manuscript excludes climate data from the final analysis owing to poor spatial resolution.    

-   *01_preparing-climate-data.Rmd:*. The aim of this script is to prepare historical climate data associated with specific locations, and to link the data to geographical covariates.  

-   *02_semi-spatial-gam-climate.Rmd:*. The overall aim here is to model historical climate variables --- mean temperature, standard deviation in temperature, and total precipitation --- as functions of geographical characteristics of the locations with which the data are associated.  

-   *03_climate-correction-layer-for-downscaling.Rmd*: In this script, we focus on quantifying the error over the entire study site between model-predicted climate and climate measures derived from remote sensing data.  

-   *04_climate-reconstruction.Rmd*: In this script, we reconstruct historical climate data using the correction layer prepared previously.  

-   *05_climate-change-analysis.Rmd*: In this script, we will process the reconstructed climate data and analyze trends in climate over time.  

-   *06_land-cover-classification.Rmd*: Here, we examine the digitized historical map of the Nilgiris and compare it to high-resolution satellite imagery to quantify landscape changes over time.  

-   *07_grassland-forest-comparisons.Rmd*: In this script, we compare changes in grassland cover and habitat over time in relation to changes in forest cover over time.  

-   *08_time-series-landscape-change.Rmd:* Using additional historical maps from 1910 and satellite imagery from 1973 and 1995, we verify patterns of grassland habitat loss.  

-   *09_occurrence-exploratory-analyses.Rmd:* Here, we explore and analyze historical bird observations, collated from natural history collections, diaries, journals and other historical sources.  

-   *10_species-relative-abundances.Rmd:* We apply statistical methods employed by Gotelli et al. (2021) to quantify differences in relative abundance over time.  

-   *11_relative-abundances-traits.Rmd:* We classify species as a grassland specialist, forest specialist or a generalist bird species and asked if their relative abundances varied over time.  

-   *12_relative-abundance-envtChange.Rmd:* Here, we test if change in grassland area is significantly associated with a change in grassland bird relative abundance over time.  

All [code](https://zenodo.org/records/15831546) and [data](https://doi.org/10.5281/zenodo.15831614) have been archived in Zenodo and can be accessed freely.  