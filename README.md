# RadoHydro

Radohydro is a small python tool which processes time series precipitation records from DWD Radolan data and among others map them on polygon shapefiles (e.g. representing subcatchments/catchments)

## Features
* Automatize download and unzip from dwd ftp server for given time interval
* Fully stream-based solution
* Automatic reprojection to destination EPSG given by user
* Equal-area approach to compute precipitation rates per basin per time
* Output 1) as csv files containing time series for each basin and <br/>2) as polygon shapefile having the precipitation rates as attribute tables

## Dependencies

* Python 3
* numpy
* geopandas
* shapely
* rasterio
* gdal

## Limitations

* until now supports only download of hourly datasets and full day downloads
* for more than 500 precipation records(~21days), output as shp is not supported (dbase problem)
* large amounts of polygons per shapefile (>10000) may lead to considerable performance decrease

```python
# Quickstart
import radohydro
radohydro.radohydro(start_date='20171230',
              end_date='20180102',
              shape_inpt='.\Examples\einzugsgebiet.shp',
              outpt_proj='epsg:25833',
              Output=True):
```

## Authors

* **Erik Nixdorf**
* **Marco Hannemann**



## Acknowledgments

* Thx to Nico Trauth for improving the code and the discussions


