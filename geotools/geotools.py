# -*- coding: utf-8 -*-
"""
a module wich writes the pythonsubroutines to process retrieved datasets
for the smart_monitoring_toolbox


"""
__author__ = "Erik Nixdorf"
__propertyof__ = "Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. "
__email__ = "erik.nixdorf@ufz.de, marco.hannemann@ufz.de"
__version__ = "0.1"
#%% load required packages
import numpy as np
import pandas as pd
import os
from osgeo import gdal
import geopandas as gp
from itertools import product
from shapely.geometry import box
import rasterio
import affine  # transformation issues
import geopandas as gpd
from rasterio import features
#%%

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
        r_na = na_value  #
    if isinstance(shape_inpt, str):
        gdfbnd = gp.read_file(shape_inpt)
    else:
        try:
            gdfbnd = shape_inpt
        except Exception:
            # if no shape is given no clipping will be done:
            print('No readable clipping file defined')
            exit()
    #Unify both CRS systems
    crs_dest=rasterio.crs.CRS.from_user_input(r_proj)
    crs_shp=rasterio.crs.CRS.from_user_input(gdfbnd.crs)
    if crs_dest.to_string() != crs_shp.to_string():
        try:
            # if not we simply use projection str
            gdfbnd = gdfbnd.to_crs(crs=crs_dest)
            print('Boundary GIS vectorfile was reprojected to Raster projection')
        except Exception as e:
            print(e)
            print('projection is not provided as crs object')
    else:
        print('projection system is similar')

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
    try:
        footprint_cells.crs = proj_crs
    except:
        print('Projection is not provided as crs object')
    return footprint_cells

def map_arraystack_on_cellgrd(array_stack,
                           stack_dates,
                           cellgrd,
                           numerator=10,
                           Output=False,
                           outpt_proj='epsg:25833',number_frmt=np.float32):
    """
    A tool which maps all retrieved precipitation grids on the poylgoncells
    """
    #sort cells
    try:
        cellgrd = cellgrd.sort_values(['Index_row', 'Index_column'])
    except Exception:
        print('No index row/column exist, generate a better radocellpolygrid')
        return None
    data_cols = list()
    #add new columns for the precipitation
    for i in range(0, len(stack_dates)):
        data_col_nms = data_cols.append(
            stack_dates[i].strftime("%y%m%d%H%M"))
        cellgrd[data_cols[i]] = np.true_divide(
            array_stack[:, :, i].reshape(-1, 1)[:, 0], numerator).astype(
                number_frmt)

    if Output:
        # try to create the directory
        try:
            os.mkdir('Data')
        except OSError:
            pass
        cellgrd = cellgrd.to_crs({'init': outpt_proj})
        cellgrd.to_file('.\Data\datacellgrid.shp')
    return cellgrd, data_col_nms


#
def map_cellgrd_on_polyg(cellgrd, shape_inpt=None, outpt_proj='epsg:25833'):
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
    # reproject both to output projection
    gdfbnd = gdfbnd.to_crs({'init': outpt_proj})
    cellgrd = cellgrd.to_crs({'init': outpt_proj})
    print('polygons reprojected to', outpt_proj)
    #calculate the area of the radocells
    cellgrd['gridcellarea'] = cellgrd['geometry'].area
    # replace gridcodes with new ID, IDs are in order of polygons
    gdfbnd['basinID'] = range(1, gdfbnd.shape[0] + 1)
    # clip radolan data with input shapefile -->most time consuming timestep
    gdfclip = gp.overlay(
        cellgrd,
        gdfbnd,
        how='intersection',
        make_valid=True,
        use_sindex=None)

    return gdfclip, gdfbnd


def compute_polyg_values(gdfclip,
                         gdfbnd,header='rainfall[mm/h]',
                         datacol_type='AllDigits',
                         Output=True,
                         outpt_proj='epsg:25833',
                         outpt_nm='radohydro'):
    """
    This Function basically sums of all cells which belong to the same basin ID (polygon
    A Weighted average approach is used
    Output is both a csv time series and and geodataframe
    precip_colmns can be either defined by giving a list of column names or by
    reading all columns which have only digits in column name (mode='AllDigits')
    """
    # first we identify the header of the precip_clms
    if datacol_type == 'AllDigits':
        datacols = [
            column for column in gdfclip.columns if column.isdigit()
        ]
    else:
        datacols = datacol_type
    #reproject both datasets
    gdfbnd = gdfbnd.to_crs({'init': outpt_proj})
    gdfclip = gdfclip.to_crs({'init': outpt_proj})
    print('polygons reprojected to', outpt_proj)
    #We compute the cellweights of each clipped cell
    cellweights = (
        gdfclip['geometry'].area / gdfclip['gridcellarea']).to_numpy()
    #get basin ID vector and also the unique values preserving the row
    basinIDs = gdfclip['basinID'].to_numpy()
    basinIDsunique = pd.unique(gdfclip['basinID'])
    #Get numpy array from the dataframe for all data as well as for the Basin ID and weights, both together does not work :-(
    datacellsextracted = gdfclip.loc[:, datacols[
        0]:datacols[-1]].to_numpy() * gdfclip['gridcellarea'][0]
    #Multiply by cellweights and add basinID
    datacellweighted = datacellsextracted * cellweights[:, None]  #http://scipy-lectures.org/advanced/advanced_numpy/#broadcasting
    #add the weights
    datacellweighted = np.hstack((datacellweighted,
                                           cellweights.reshape(-1, 1)))
    #add basinIDs as a new column
    datacellweighted = np.hstack((datacellweighted,
                                           basinIDs.reshape(-1, 1)))
    #Finally we sum each row to get value for each polygon and each time step https://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.reduceat.html
    polyg_values = np.add.reduceat(
        datacellweighted,
        np.r_[0, np.where(np.diff(datacellweighted[:, -1]))[0] + 1],
        axis=0)
    #Compute the average weight
    polyg_values[:, -2] = polyg_values[:, -2] * basinIDsunique / polyg_values[:,-1]
    #replace sum of unique basin id by actual basin id
    polyg_values[:, -1] = basinIDsunique
    #retrieve output_area
    polygonarea = np.vstack((gdfbnd['basinID'].to_numpy(),
                             gdfbnd.area.to_numpy())).T
    polyg_values = np.hstack(
        (polyg_values,
         polygonarea[(basinIDsunique.flatten() - 1)][:, 1].reshape(-1, 1)))
    #compute to intense parameter again by divding by area
    polyg_values[:, :len(datacols)] = (
        polyg_values[:, :len(datacols)].T /
        polyg_values[:, -1]).T
    # output if desired
    if Output == True:
        # try to create the directory
        try:
            os.mkdir('Data')
        except OSError:
            pass

        #write out the numpy array
        for polyg_value in polyg_values:
            with open(
                    '.\Data\\' + outpt_nm + '_{!s}.csv'.format(polyg_value[-2]),
                    'w',
                    newline='') as fout:
                fout.write('basin ID: {:d}\n'.format(int(polyg_value[-2])))
                fout.write('average_cellweight: {0:.3f}\n'.format(polyg_value[-3]))
                fout.write('basin_area: {0:.3f}\n'.format(polyg_value[-1]))
                fout.write('Time[yymmddhh],'+header)
                np.savetxt(
                    fout,
                    np.column_stack((datacols,
                                     np.around(
                                         polyg_value[:len(datacols)],
                                         decimals=3))),
                    delimiter=",",
                    fmt="%s")

        if len(datacols) < 500:
            #Write out the updated Polygonshape
            #add new information to old dataframe and write it out as shp
            polyg_valuessorted = polyg_values[
                polyg_values[:, -2].argsort()]
            for i in range(0, len(datacols)):
                gdfbnd[datacols[i]] = polyg_valuessorted[:, i]
            gdfbnd['BasinIDNew'] = polyg_valuessorted[:, -2]
            gdfbnd['average_cellweight'] = polyg_valuessorted[:, -3]
            gdfbnd['basin_area'] = polyg_valuessorted[:, -1]
            gdfbnd.to_file(
                os.path.join(os.getcwd(), 'Data', 'polygon_values_'+header+'.shp'))

        print('\nOutput saved to /Data')

    return polyg_values


def rasterizegeo(shp_inpt='.\example\einzugsgebiet.shp',
                 pixel_sz=(100, 100),
                 atrbt_nm='gridcode',
                 Output=True,
                 Output_nm='test.tif'):
    """
    Created on Tue Oct 22 12:15:13 2019
    A small tool which converts shapefiles and a selectable feature to rasterfiles
    Heavily exploits rasterize.features.rasterize functionality
    @author: nixdorf
    """


    #%% start processinghttps://sigon.gitlab.io/post/2019-05-02-rasterize/
    #open geometry to geopandas
    if isinstance(shp_inpt, str):
        gdf_inpt = gpd.GeoDataFrame.from_file(shp_inpt)
    else:
        gdf_inpt = shp_inpt.copy()

    rst_transform = affine.Affine(pixel_sz[0], 0, gdf_inpt.total_bounds[0], 0,
                                  -pixel_sz[1], gdf_inpt.total_bounds[3])
    rst_shape = (int(
        (gdf_inpt.total_bounds[2] - gdf_inpt.total_bounds[0]) / pixel_sz[0]),
                 int((gdf_inpt.total_bounds[3] - gdf_inpt.total_bounds[1]) /
                     pixel_sz[1]))
    rst_data = features.rasterize(
        shapes=gdf_inpt[['geometry', atrbt_nm]].values.tolist(),
        out_shape=(rst_shape[1], rst_shape[0]),
        transform=rst_transform,
        fill=-9999)

    if Output:
        with rasterio.open(
                Output_nm,
                'w',
                driver='GTiff',
                transform=rst_transform,
                dtype=rasterio.float64,
                crs=gdf_inpt.crs['init'],
                count=1,
                nodata=-9999,
                width=rst_shape[0],
                height=rst_shape[1]) as dst:
            dst.write(rst_data.astype(float), indexes=1)
