# What is rasterdiv?

*rasterdiv* is a package for the R statistical software and environment. It aims to provide functions to apply Information Theory based diversity indeces on RasterLayer or numerical matrices, such as Shannon's entropy or Cumulative Residual Entropy (CRE).

## Rasterdiv basics. Derive indices of diversity from NDVI.

Here, we show hot to use **rasterdiv** to derive global series of indices of diversity based on Information Theory. The input dataset is the Copernicus Long-term (1999-2017) average Normalise Difference Vegetation Index for the 21st of June (copNDVI).

## Overview
A RasterLayer called copNDVI is loaded together with package **rasterdiv**. *copNDVI* is at 8-bit meaning that pixel values range from 0 to 255. You could *stretch* it to a more familiar (-1,1) range using `raster::stretch(copNDVI,minv=-1,maxv=1)` .  

## Reclassify NDVI 
Pixels with values 253, 254 and 255 (water) will be set as NA's.

```{r}
copNDVI <- reclassify(copNDVI, cbind(253, 255, NA), right=TRUE)
```

## Resample NDVI to coarser resolution 
To speed up the calculation, the RasterLayer will be "resampled" at a resolution 10 times coarser than original.

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

## Derive Hill's and Renyi's indices 
Renyi and Hill's indices are calculated from copNDVI with a moving window of 81 pixels (9 px side) for alpha values from 0 to 2 every 0.5. In addition, we set `na.tolerance=0.5`, meaning that pixels whose 9x9 px window has more than 50% of NA values will be set to NA.

```{r echo = T, results = 'hide', warning=FALSE, message=FALSE}
#Shannon's Diversity
sha <- Shannon(copNDVIlr,window=9,np=1,na.tolerance=0.1)

#Pielou's Evenness
pie <- Pielou(copNDVIlr,window=9,np=1,na.tolerance=0.1)

#Berger-Parker
ber <- BergerParker(copNDVIlr,window=9,np=1,na.tolerance=0.1)

#Rao's quadratic Entropy
rao <- Rao(copNDVIlr,window=9,np=1,na.tolerance=0.1,dist_m="euclidean",shannon=FALSE)

#Cumulative Residual Entropy
cre <- CRE(copNDVIlr,window=9,np=1,na.tolerance=0.1)

#Hill's numbers
hil <- Hill(copNDVIlr,window=9,np=1,na.tolerance=0.1,alpha=seq(0,2,0.5))

#Renyi
ren <- Renyi(copNDVIlr,window=9,np=1,na.tolerance=0.1,alpha=seq(0,2,0.5))
```

## Transform output matrices to RasterLayer
All functions in *rasterdiv* output numerical matrices or list of numerical matrices. Therefore, they need to be placed in space to be transformed in RasterLayer. This is done by using the input NDVI RasterLayer as template.

```{r}
shara <- raster(sha,template=copNDVIlr)
piera <- raster(pie,template=copNDVIlr)
berra <- raster(ber,template=copNDVIlr)
raora <- raster(rao,template=copNDVIlr)
crera <- raster(cre,template=copNDVIlr)
hilra <- lapply(hil, function(x) raster(x,template=copNDVIlr))
renra <- lapply(ren, function(x) raster(x,template=copNDVIlr))
```

## Visualise RasterLayers

```{r fig02}
#Shannon's Diversity
levelplot(shara,main="Shannon's entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig03}
#Pielou's Evenness
levelplot(piera,main="Pielou's evenness from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig04}
#Berger-Parker' Index
levelplot(berra,main="Berger-Parker's index from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig05}
#Rao's quadratic Entropy
levelplot(crera,main="Rao's quadratic entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```


```{r fig06}
#Cumulative Residual Entropy
levelplot(crera,main="Cumulative Resiudal Entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,1,1), ylim=c(-60,75), margin = list(draw = TRUE))
```

```{r fig07}
#Hill's numbers (alpha=0, 1, 1.5 and 2)
levelplot(stack(hilra),main="Hill's numbers from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,5,1),names.attr=paste("alpha",seq(0,2,0.5),sep=" "), ylim=c(-60,75))
```

```{r fig08}
#Renyi' Index (alpha=0, 1, 1.5 and 2)
levelplot(stack(renra),main="Renyi's entropy from Copernicus NDVI 5 km (9 px-side moving window)",as.table = T,layout=c(0,5,1),names.attr=paste("alpha",seq(0,2,0.5),sep=" "), ylim=c(-60,75), margin = list(draw = FALSE))
```
