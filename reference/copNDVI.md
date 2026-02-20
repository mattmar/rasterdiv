# Copernicus Long Term (1999-2017) NDVI Overview (5km)

A \`SpatRaster\` (EPSG: 4326) of the global average NDVI value per pixel
for the 21st of June over the period 1999-2017.

## Format

A \`SpatRaster\` containing the following elements:

- NDVI:

  Normalised Difference Vegetation Index value (0-255) for each 5 km
  pixel. This index provides an indication of the presence of live green
  vegetation in the area.

## Source

<https://land.copernicus.eu/en/products/vegetation>

## Details

This dataset provides a long-term overview of the Normalised Difference
Vegetation Index (NDVI) across the globe. Each pixel represents a 5 km
area, with NDVI values ranging from 0 to 255.

## References

<https://land.copernicus.eu/en/products/vegetation>

## Examples

``` r
copNDVI <- load_copNDVI()
```
