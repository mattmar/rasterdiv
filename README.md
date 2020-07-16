# rasterdiv

[![cran version](http://www.r-pkg.org/badges/version/rasterdiv)](https://cran.r-project.org/package=rasterdiv)
[![downloads](http://cranlogs.r-pkg.org/badges/rasterdiv)](http://cranlogs.r-pkg.org/badges/rasterdiv)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/rasterdiv)](http://cranlogs.r-pkg.org/badges/grand-total/rasterdiv)

# How to install?

### Stable versions can be installed from the CRAN:

```{r}
install.packages("rasterdiv")
```

### Development versions stored in the Git repository can be installed with:

```{r}
library(devtools)
install_github("mattmar/rasterdiv")
```

## What is rasterdiv?

*rasterdiv* is a package for the R statistical software and environment. It aims to provide functions to apply Information Theory based diversity indeces on RasterLayer or numerical matrices, such as Shannon's entropy or Cumulative Residual Entropy (CRE).

## Rasterdiv basics. Derive indices of diversity from NDVI.

Here, we show hot to use **rasterdiv** to derive global series of indices of diversity based on Information Theory. The input dataset is the Copernicus Long-term (1999-2017) average Normalise Difference Vegetation Index for the 21st of June (copNDVI).

## Overview
A RasterLayer called copNDVI is loaded together with package **rasterdiv**. *copNDVI* is at 8-bit meaning that pixel values range from 0 to 255. You could *stretch* it to a more familiar (-1,1) range using `raster::stretch(copNDVI,minv=-1,maxv=1)` .  

## Reclassify NDVI 
Pixels with values 253, 254 and 255 (water) will be set as NA's.

```{r}
copNDVI <- reclassify(copNDVI, cbind(252, 255, NA), right=TRUE)
```

## Resample NDVI to coarser resolution 
To speed up the calculation, the RasterLayer will be "resampled" at a resolution 20 times coarser than original.

```{r}
#Resample using raster::aggregate and a linear factor of 10
copNDVIlr <- raster::aggregate(copNDVI, fact=20)
#Set float numbers as integers to further speed up the calculation
storage.mode(copNDVIlr[]) = "integer"
```

## Compare NDVI low and high resolution

```{r fig01}
levelplot(copNDVI,layout=c(0,1,1), main="NDVI 21st of June 1999-2017 - ~8km pixel resolution")
levelplot(copNDVIlr,layout=c(0,1,1), main="NDVI 21st of June 1999-2017 - ~150km pixel resolution")`
```

## Compute all indexes in rasterdiv
**rasterdiv** allows the computation of 8 diversity indexes based on information theory. In the following section, all these indexes will be computed for *copNDVIlr* using a moving window of 81 pixels (9 px side). Alpha values for the Hill, Renyi and parametric Rao indexes will be set from 0 to 2 every 0.5. In addition, we will set `na.tolerance=0.1`, meaning that all moving windows with more than 10% of pixels equal NA will be set to NA.

```{r echo = T, results = 'hide', warning=FALSE, message=FALSE}
#Shannon's Diversity
sha <- Shannon(copNDVIlr,window=9,na.tolerance=0.1,np=1)

#Pielou's Evenness
pie <- Pielou(copNDVIlr,window=9,na.tolerance=0.1,np=1)

#Berger-Parker's Index
ber <- BergerParker(copNDVIlr,window=9,na.tolerance=0.1,np=1)

#Rao's quadratic Entropy
rao <- Rao(copNDVIlr,window=9,na.tolerance=0.1,dist_m="euclidean",shannon=FALSE,np=1)

#Parametric Rao's quadratic entropy with alpha ranging from 1 to 5
prao <- paRao(copNDVIlr,window=9,alpha=1:5,na.tolerance=0.1,dist_m="euclidean",np=1)

#Cumulative Residual Entropy
cre <- CRE(copNDVIlr,window=9,na.tolerance=0.1,np=1)

#Hill's numbers
hil <- Hill(copNDVIlr,window=9,alpha=seq(0,2,0.5),na.tolerance=0.1,np=1)

#Renyi's Index
ren <- Renyi(copNDVIlr,window=9,alpha=seq(0,2,0.5),na.tolerance=0.1,np=1)
```

## Visualise RasterLayers

```{r fig02}
#Shannon's Diversity
levelplot(sha,main="Shannon's entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig03}
#Pielou's Evenness
levelplot(pie,main="Pielou's evenness from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig04}
#Berger-Parker' Index
levelplot(ber,main="Berger-Parker's index from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig05}
#Rao's quadratic Entropy
levelplot(rao,main="Rao's quadratic entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig06}
#Parametric Rao's quadratic Entropy
levelplot(stack(prao),main="Parametric Rao's quadratic entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,5,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig07}
#Cumulative Residual Entropy
levelplot(cre,main="Cumulative Residual Entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig08}
#Hill's numbers (alpha=0, 1, 1.5 and 2)
levelplot(stack(hil),main="Hill's numbers from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,5,1),names.attr=paste("alpha",seq(0,2,0.5),sep=" "), ylim=c(-60,75))
```

```{r fig09}
#Renyi' Index (alpha=0, 1, 1.5 and 2)
levelplot(stack(ren),main="Renyi's entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,5,1),names.attr=paste("alpha",seq(0,2,0.5),sep=" "), ylim=c(-60,75), margin = list(draw = FALSE))
```
