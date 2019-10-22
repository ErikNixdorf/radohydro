#!/usr/bin/env python3.5
"""This script intersects a given shapefile of one or multiple basins with the Radolan ASCII data by DWD and returns a
   .csv with the rainfall intensity for each basin in the given period as output. Optionally the clipped shape file
   and the geoTIFFs can be requested
   Major Changes in V0.3:
       -Performance by vectorization
   Major Changes in V04:
       -Performance (factor 3) by having fully streambased solution
       -Script splitted into functions to increase flexibility
"""

__author__ = "Erik Nixdorf, Marco Hannemann"
__propertyof__ = "Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. "
__email__ = "erik.nixdorf@ufz.de, marco.hannemann@ufz.de"
__version__ = "0.4"
import time
import os
from ftplib import FTP
import tarfile
from datetime import datetime
from dateutil.relativedelta import relativedelta
from itertools import product
import re
import numpy as np
import pandas as pd
from osgeo import gdal
import geopandas as gp
from shapely.geometry import box
from io import BytesIO
from rasterio.io import MemoryFile


# generator to create list of dates
def daterange(start, end, delta):
    """
    generator to create list of dates
    """
    current = start
    while current <= end:
        yield current
        current += delta


# define radolan projection
def get_radoproj(earthradius=6370040):
    """
    Defines the native radoproj
    """
    # Native projection ("RADOLAN-Projektion"):
    R = earthradius
    # Erdradius (Kugel)
    R = "+a=" + str(R) + " +b=" + str(R)
    # zu 'R' zusammenfassen -> Kugel
    nat_proj = "+proj=stere +lon_0=10.0 +lat_0=90.0 +lat_ts=60.0 " + R + " +units=m"
    return nat_proj


# a function which clips through raster provides as numpy array, too


def buffered_raster_clipping(raster_inpt,
                             shape_inpt=None,
                             raster_transfrm=None,
                             raster_proj=None,
                             buffrcllsz=2,
                             na_value=-1):
    """
    This function clipped a raster dataset with shapefile, even if rasterdata
    is numpy array only
    raster_inpt could be either a georeferenced raster or a numpy array 
    (requires raster_transfrm and raster_crs)
    shape_inpt (currently either string or geodataframe)
    Output:
        r_clip_data: clipped_numpy array
        r_clip_transform: New transform in gdal transform style
        colums: Positions in original dataarray for begin/end clip (column)
        rows: Positions in original dataarray for begin/end clip (row)
    """
    # read the rasterfile
    if isinstance(raster_inpt, str):
        rds = gdal.Open(raster_inpt)
        r_proj = rds.GetProjection()
        r_transform = rds.GetGeoTransform()
        r_data = np.array(rds.GetRasterBand(1).ReadAsArray())
        r_na = rds.GetRasterBand(1).GetNoDataValue()
    else:
        r_data = raster_inpt
        r_proj = raster_proj
        r_transform = raster_transfrm
        r_na = na_value
    # read the shapefile
    if isinstance(shape_inpt, str):
        gdfbnd = gp.read_file(shape_inpt)
    else:
        try:
            gdfbnd = shape_inpt
        except Exception:
            # if no shape is given no clipping will be done:
            print('No readable clipping file defined')
            exit()
    #reproject depending on what have entered as information
    if r_proj.find('epsg')!=-1:
        try: 
            epsg_code=int(r_proj[5:])
            gdfbnd.to_crs(epsg=epsg_code,inplace=True)
            print('Boundary GIS vectorfile was reprojected to Raster projection')
        except Exception as e:
                print('error converting epsg code')
                print(e)
    else:
        #if not we simply use projection str
        gdfbnd=gdfbnd.to_crs(crs=r_proj)
        print('Boundary GIS vectorfile was reprojected to Raster projection')
    #gdfbnd.to_file('test_example.shp')
    # Find the clipping window
    cellsize = abs(min((abs(r_transform[1]), abs(r_transform[5]))))
    BoundsClipDs = np.add(gdfbnd.geometry.total_bounds,
                          np.array([-1, -1, 1, 1]) * buffrcllsz * cellsize)
    # find the positions of the new boundary
    colums = ((
        BoundsClipDs[[0, 2]] - r_transform[0]) / r_transform[1]).astype(int)
    rows = ((
        BoundsClipDs[[1, 3]] - r_transform[3]) / r_transform[5]).astype(int)
    # get the new dataarray (first row than columns)
    r_clip_data = r_data[rows[1]:rows[0], colums[0]:colums[1]]
    #write a new transform
    r_clip_transform = (r_transform[0] + colums[0] * r_transform[1],
                        r_transform[1], 0,
                        r_transform[3] + rows[1] * r_transform[5], 0,
                        r_transform[5])

    return r_clip_data, r_clip_transform, colums, rows


# A tool to create footprint cells (polygons of each radolan dataset
def create_footprint_cells(transform=None, data_size=None, proj_crs=None):
    """
    This tool creates footprint cells(rectangles) from given rastermetadata
    Output is a geodataframe
    data_size= No of rows/columns,e.g using the shape function
    """
    ulx, xres, xskew, uly, yskew, yres = transform

    # ULCorner=src.transform * (0, 0)
    CellBoundsLR = np.linspace(ulx, ulx + (data_size[1] * xres),
                               data_size[1] + 1)
    CellBoundsUB = np.linspace(uly, uly + (data_size[0] * yres),
                               data_size[0] + 1)

    # create boundaries of each cell
    cellbounds = np.zeros((data_size[1] * data_size[0], 6))
    cellbounds[:, (0, 3)] = np.array(
        list(product(CellBoundsLR[:-1], CellBoundsUB[:-1])))
    cellbounds[:, (2, 1)] = np.array(
        list(product(CellBoundsLR[1:], CellBoundsUB[1:])))
    cellbounds[:, (4, 5)] = np.array(
        list(
            product(
                range(0,
                      len(CellBoundsLR) - 1), range(0,
                                                    len(CellBoundsUB) - 1))))

    # Create a pandas dataframe
    df = pd.DataFrame({
        'left': cellbounds[:, 0],
        'bottom': cellbounds[:, 1],
        'right': cellbounds[:, 2],
        'top': cellbounds[:, 3],
        'Index_column': cellbounds[:, 4].astype(int),
        'Index_row': cellbounds[:, 5].astype(int)
    })
    b = [
        box(l, b, r, t)
        for l, b, r, t in zip(cellbounds[:, 0], cellbounds[:, 1],
                              cellbounds[:, 2], cellbounds[:, 3])
    ]
    footprint_cells = gp.GeoDataFrame(df, geometry=b)
    footprint_cells.crs = proj_crs
    return footprint_cells


def radoname_to_date(filename, date_type):
    """
    Filenames of the hourly datasets to datestr
    """

    datestr = str()
    for number in re.findall(r'\d+', filename):
        datestr = datestr + number
    if date_type == 'hours':
        date = datetime.strptime(datestr, '%Y%m%d%H')
    if date_type == 'days':
        date = datetime.strptime(datestr, '%Y%m%d')
    if date_type == 'minutes':
        date = datetime.strptime(datestr, '%Y%m%d%H%M')
    return date


# The Actual downloader for Radolan to archive
def rado_io(start_date='20171231',
            end_date='20180101',
            shapefile='.\Examples\einzugsgebiet.shp'):
    """
    The Actual downloader for Radolan to stream and clipping
    """
    # create the dates
    start_date = datetime.strptime(start_date, '%Y%m%d')
    end_date = datetime.strptime(end_date, '%Y%m%d')

    dts = [
        dt.strftime('%Y%m%d')
        for dt in daterange(start_date, end_date, relativedelta(days=1))
    ]
    dts_historical = [
        dt.strftime('%Y%m%d')
        for dt in daterange(start_date, end_date, relativedelta(months=1))
    ]
    years = list(range(start_date.year, end_date.year + 1))

    #define the radolan_projection
    rado_proj = get_radoproj()

    # Connect to the Server
    server = 'opendata.dwd.de'
    connected = False

    while not connected:
        try:
            ftp = FTP(server)
            ftp.login()
            connected = True

        except:
            time.sleep(5)
            print('Reconnect to Server')
            pass

    # Set Initial run to true
    initDf = True
    # avoid downloading 2018 multiple times
    if (datetime.now().year in years) and (datetime.now().year - 1 in years):
        years = years[:-1]

    for year in years:
        #check whether data is recent
        if (year == datetime.now().year or year == (datetime.now().year - 1)):
            ftp.cwd(
                '/climate_environment/CDC/grids_germany/hourly/radolan/recent/asc/'
            )
            files = ftp.nlst()
            for dt, file in product(dts, files):
                if dt in file:
                    print('Retrieving {}...'.format(file))
                    archive = BytesIO()
                    ftp.retrbinary("RETR " + file, archive.write)
                    archive.seek(0)
                    archive_daily = tarfile.open(fileobj=archive)
                    #extract file to bytestream
                    for member in archive_daily.getmembers():
                        radolan_io = archive_daily.extractfile(member.name)
                        with MemoryFile(radolan_io) as memfile:
                            with memfile.open() as rado_ds:
                                rado_data = rado_ds.read()[0]
                                #depending whether we get first dataset or not we do
                                #different calculations
                                if initDf:
                                    NaN_Value = rado_ds.nodata
                                    afn_transform = rado_ds.transform
                                    rado_transform = (afn_transform[2],
                                                      afn_transform[0], 0,
                                                      afn_transform[5], 0,
                                                      afn_transform[4])
                                    # do the complicated buffer clipping
                                    rado_clip_data, rado_clip_transform, cols, rows = buffered_raster_clipping(
                                        rado_data,
                                        shape_inpt=shapefile,
                                        raster_transfrm=rado_transform,
                                        raster_proj=rado_proj)

                                    #generate the footprint cells
                                    radocells = create_footprint_cells(
                                        transform=rado_clip_transform,
                                        data_size=rado_clip_data.shape,
                                        proj_crs=rado_proj)
                                    #initialize the merged dataset
                                    rado_stacked_data = rado_clip_data
                                    #the dates
                                    rado_dates = [
                                        radoname_to_date(
                                            member.name, 'minutes')
                                    ]
                                    initDf = False
                                # if we initialised already, computation is easy
                                else:
                                    rado_clip_data = rado_data[rows[1]:rows[0],
                                                               cols[0]:cols[1]]
                                    rado_stacked_data = np.dstack(
                                        (rado_stacked_data, rado_clip_data))
                                    rado_dates.append(
                                        radoname_to_date(
                                            member.name, 'minutes'))
                    print('Processing {}...finished'.format(file))
        # the historical case
        else:
            ftp.cwd(
                '/climate_environment/CDC/grids_germany/hourly/radolan/historical/asc/{}/'.
                format(year))
            files = ftp.nlst()
            for dt, file in product(dts_historical, files):

                if dt[:-2] in file:
                    print('Retrieving {}...'.format(file))
                    archive = BytesIO()
                    ftp.retrbinary("RETR " + file, archive.write)
                    archive.seek(0)
                    archive_monthly = tarfile.open(fileobj=archive)
                    #double_zipped so we need to get daily archives
                    for members_daily in archive_monthly.getnames():
                        # check whether the day is in our time span
                        members_date = radoname_to_date(members_daily, 'days')
                        if members_date >= start_date and members_date <= end_date:
                            tar_daily = archive_monthly.extractfile(
                                members_daily)
                            tar_daily.seek(0)
                            archive_daily = tarfile.open(fileobj=tar_daily)
                            #extract file to bytestream
                            for member in archive_daily.getmembers():
                                radolan_io = archive_daily.extractfile(
                                    member.name)
                                with MemoryFile(radolan_io) as memfile:
                                    with memfile.open() as rado_ds:
                                        rado_data = rado_ds.read()[0]
                                        #depending whether we get first dataset or not we do
                                        #different calculations
                                        if initDf:
                                            NaN_Value = rado_ds.nodata
                                            afn_transform = rado_ds.transform
                                            rado_transform = (afn_transform[2],
                                                              afn_transform[0],
                                                              0,
                                                              afn_transform[5],
                                                              0,
                                                              afn_transform[4])
                                            # do the complicated buffer clipping
                                            rado_clip_data, rado_clip_transform, cols, rows = buffered_raster_clipping(
                                                rado_data,
                                                shape_inpt=shapefile,
                                                raster_transfrm=rado_transform,
                                                raster_proj=rado_proj)

                                            #generate the footprint cells
                                            radocells = create_footprint_cells(
                                                transform=rado_clip_transform,
                                                data_size=rado_clip_data.shape,
                                                proj_crs=rado_proj)
                                            #initialize the merged dataset
                                            rado_stacked_data = rado_clip_data
                                            #the dates
                                            rado_dates = [
                                                radoname_to_date(
                                                    member.name, 'minutes')
                                            ]
                                            initDf = False
                                        # if we initialised already, computation is easy
                                        else:
                                            rado_clip_data = rado_data[rows[
                                                1]:rows[0], cols[0]:cols[1]]
                                            rado_stacked_data = np.dstack(
                                                (rado_stacked_data,
                                                 rado_clip_data))
                                            rado_dates.append(
                                                radoname_to_date(
                                                    member.name, 'minutes'))
                    print('Processing {}...finished'.format(file))
    ftp.quit()
    return rado_stacked_data, rado_dates, radocells


# Quit
def map_radonum_on_cellgrd(rado_stacked_data,
                           rado_dates,
                           radocells,
                           numerator=10,
                           Output=True):
    """
    A tool which maps all retrieved precipitation grids on the poylgoncells
    """
    #sort cells
    try:
        radocells = radocells.sort_values(['Index_row', 'Index_column'])
    except Exception:
        print('No index row/column exist, generate a better radocellpolygrid')
        return None
    precip_cols = list()
    #add new columns for the precipitation
    for i in range(0, len(rado_dates)):
        precip_col_nms = precip_cols.append(
            rado_dates[i].strftime("%y%m%d%H%M"))
        radocells[precip_cols[i]] = np.true_divide(
            rado_stacked_data[:, :, i].reshape(-1, 1)[:, 0], numerator)
    if Output:
        # try to create the directory
        try:
            os.mkdir('Data')
        except OSError:
            pass
        radocells.to_file('.\Data\Radoprecipitationgrid.shp')
    return radocells, precip_col_nms


#
def map_cellgrd_on_polyg(radocells, shape_inpt=None, outpt_proj='epsg:25833'):
    """
    Maps the radocellgrid on a shapefile (e.g.basin)
    radocells should be a geodataframe and shapefile can be either a 
    geodataframe or a external shapefile
    Data will be reproject on radocells crs
    """
    # read the shapefile
    if isinstance(shape_inpt, str):
        gdfbnd = gp.read_file(shape_inpt)
    else:
        try:
            gdfbnd = shape_inpt
        except Exception:
            # if no shape is given no clipping will be done:
            print('No readable clipping file defined')
            exit()
    # reproject
    gdfbnd = gdfbnd.to_crs(crs=outpt_proj)
    radocells = radocells.to_crs(crs=outpt_proj)
    print('polygons reprojected to', outpt_proj)
    #calculate the area of the radocells
    radocells['gridcellarea'] = radocells['geometry'].area
    # replace gridcodes with new ID, IDs are in order of polygons
    gdfbnd['basinID'] = range(1, gdfbnd.shape[0] + 1)
    # clip radolan data with input shapefile -->most time consuming timestep
    gdfclip = gp.overlay(
        radocells,
        gdfbnd,
        how='intersection',
        make_valid=True,
        use_sindex=None)

    return gdfclip, gdfbnd


def compute_polyg_precip(gdfclip,
                         gdfbnd,
                         precip_colmns='AllDigits',
                         Output=True,
                         outpt_proj='epsg:25833',outpt_nm='radohydro')):
    """
    This Function basically sums of all cells which belong to the same basin ID (polygon
    A Weighted average approach is used
    Output is both a csv time series and and geodataframe
    precip_colmns can be either defined by giving a list of column names or by
    reading all columns which have only digits in column name (mode='AllDigits')
    """
    # first we identify the header of the precip_clms
    if precip_colmns == 'AllDigits':
        precipitationCols = [
            column for column in gdfclip.columns if column.isdigit()
        ]
    else:
        precipitationCols = precip_colmns
    #reproject both datasets
    gdfbnd = gdfbnd.to_crs(crs=outpt_proj)
    gdfclip = gdfclip.to_crs(crs=outpt_proj)
    print('polygons reprojected to', outpt_proj)
    #We compute the cellweights of each clipped cell
    cellweights = (
        gdfclip['geometry'].area / gdfclip['gridcellarea']).to_numpy()
    #get basin ID vector and also the unique values preserving the row
    basinIDs = gdfclip['basinID'].to_numpy()
    basinIDsunique = pd.unique(gdfclip['basinID'])
    #Get numpy array from the dataframe for all precipitation data as well as for the Basin ID and weights, both together does not work :-(
    precipitationcellsextracted = gdfclip.loc[:, precipitationCols[
        0]:precipitationCols[-1]].to_numpy() * gdfclip['gridcellarea'][0]
    #Multiply by cellweights and add basinID
    precipitationcellweighted = precipitationcellsextracted * cellweights[:,
                                                                          None]  #http://scipy-lectures.org/advanced/advanced_numpy/#broadcasting
    #add the weights
    precipitationcellweighted = np.hstack((precipitationcellweighted,
                                           cellweights.reshape(-1, 1)))
    #add basinIDs as a new column
    precipitationcellweighted = np.hstack((precipitationcellweighted,
                                           basinIDs.reshape(-1, 1)))
    #Finally we sum each row to get value for each basin and each time step https://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.reduceat.html
    precipitationbasin = np.add.reduceat(
        precipitationcellweighted,
        np.r_[0, np.where(np.diff(precipitationcellweighted[:, -1]))[0] + 1],
        axis=0)
    #Compute the average weight
    precipitationbasin[:,
                       -2] = precipitationbasin[:,
                                                -2] * basinIDsunique / precipitationbasin[:,
                                                                                          -1]
    #replace sum of unique basin id by actual basin id
    precipitationbasin[:, -1] = basinIDsunique
    #retrieve output_area
    polygonarea = np.vstack((gdfbnd['basinID'].to_numpy(),
                             gdfbnd.area.to_numpy())).T
    precipitationbasin = np.hstack(
        (precipitationbasin,
         polygonarea[(basinIDsunique.flatten() - 1)][:, 1].reshape(-1, 1)))
    #compute to mm/h
    precipitationbasin[:, :len(precipitationCols)] = (
        precipitationbasin[:, :len(precipitationCols)].T /
        precipitationbasin[:, -1]).T
    # output if desired
    if Output == True:
        # try to create the directory
        try:
            os.mkdir('Data')
        except OSError:
            pass

        #write out the numpy array
        for basin in precipitationbasin:
            with open(
                    '.\Data\\'+outpt_nm+'_{!s}.csv'.format(basin[-2]),
                    'w',
                    newline='') as fout:
                fout.write('basin ID: {:d}\n'.format(int(basin[-2])))
                fout.write('average_cellweight: {0:.3f}\n'.format(basin[-3]))
                fout.write('basin_area: {0:.3f}\n'.format(basin[-1]))
                fout.write('Time[yymmddhh],rainfall[mm/h]\n')
                np.savetxt(
                    fout,
                    np.column_stack((precipitationCols,
                                     np.around(
                                         basin[:len(precipitationCols)],
                                         decimals=3))),
                    delimiter=",",
                    fmt="%s")

        if len(precipitationCols) < 500:
            #Write out the updated Polygonshape
            #add new information to old dataframe and write it out as shp
            precipitationbasinsorted = precipitationbasin[
                precipitationbasin[:, -2].argsort()]
            for i in range(0, len(precipitationCols)):
                gdfbnd[precipitationCols[i]] = precipitationbasinsorted[:, i]
            gdfbnd['BasinIDNew'] = precipitationbasinsorted[:, -2]
            gdfbnd['average_cellweight'] = precipitationbasinsorted[:, -3]
            gdfbnd['basin_area'] = precipitationbasinsorted[:, -1]
            gdfbnd.to_file(
                os.path.join(os.getcwd(), 'Data', 'precip_polygon.shp'))

        print('\nOutput saved to /Data')

    return precipitationbasin


def rasterizegeo(shp_inpt='.\example\einzugsgebiet.shp',pixel_sz=(100,100),atrbt_nm='gridcode',Output=True,Output_nm='test.tif'):
    """
    Created on Tue Oct 22 12:15:13 2019
    A small tool which converts shapefiles and a selectable feature to rasterfiles
    Heavily exploits rasterize.features.rasterize functionality
    @author: nixdorf
    """
    
    #%% import external libraries
    import rasterio
    import affine # transformation issues
    import geopandas as gpd
    from rasterio import features

    #%% start processinghttps://sigon.gitlab.io/post/2019-05-02-rasterize/
    #open geometry to geopandas
    if isinstance(shp_inpt,str):
        gdf_inpt=gpd.GeoDataFrame.from_file(shp_inpt)
    else:
        gdf_inpt=shp_inpt.copy()
    
    rst_transform=affine.Affine(pixel_sz[0],0,gdf_inpt.total_bounds[0],0,-pixel_sz[1],gdf_inpt.total_bounds[3])
    rst_shape=(int((gdf_inpt.total_bounds[2]-gdf_inpt.total_bounds[0])/pixel_sz[0]),int((gdf_inpt.total_bounds[3]-gdf_inpt.total_bounds[1])/pixel_sz[1]))
    rst_data = features.rasterize(shapes = gdf_inpt[['geometry', atrbt_nm]].values.tolist(),out_shape=(rst_shape[1],rst_shape[0]), transform=rst_transform,fill=-9999)
    
    
    
    if Output:    
        with rasterio.open(
            Output_nm, 'w',
            driver='GTiff',
            transform = rst_transform,
            dtype=rasterio.float64,
            crs=gdf_inpt.crs['init'],
            count=1,
            nodata=-9999,
            width=rst_shape[0],
            height=rst_shape[1]) as dst:
                dst.write(rst_data.astype(float), indexes=1)
                


#define a main function which calls all other subfunctions
def radohydro(start_date='20171230',
              end_date='20180102',
              shape_inpt='.\Examples\einzugsgebiet.shp',
              outpt_proj='epsg:25833',
              Output=True):
    """
    Couples the four main function for the entire workflow
    """

    rado_stacked_data, rado_dates, radocellgrid = rado_io(
        start_date=start_date, end_date=end_date, shapefile=shape_inpt)
    radocell_precip, precip_col_nms = map_radonum_on_cellgrd(
        rado_stacked_data,
        rado_dates,
        radocellgrid,
        numerator=10,
        Output=False)
    radocell_clip, gdf_shape = map_cellgrd_on_polyg(
        radocell_precip, shape_inpt=shape_inpt, outpt_proj=outpt_proj)
    compute_polyg_precip(
        radocell_clip,
        gdf_shape,
        precip_colmns='AllDigits',
        Output=Output,
        outpt_proj=outpt_proj)
    return None


#%% Test our functions
if __name__ == "__main__":
    radohydro()
