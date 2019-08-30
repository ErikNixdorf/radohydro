#!/usr/bin/env python3.5

"""This script intersects a given shapefile of one or multiple basins with the Radolan ASCII data by DWD and returns a
   .csv with the rainfall intensity for each basin in the given period as output. Optionally the clipped shape file
   and the geoTIFFs can be requested
   Major Changes in V0.3->Performance by vectorization
"""

__author__ = "Erik Nixdorf, Marco Hannemann"
__propertyof__="Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. "
__email__ = "erik.nixdorf@ufz.de, marco.hannemann@ufz.de"
__version__ = "0.3"
import time
import os
import shutil
import glob
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

# generator to create list of dates
def daterange(start, end, delta):
    current = start
    while current <= end:
        yield current
        current += delta

# get and unpack radolan data to /RadolanData
def processradolan(start_date, end_date, shapefile, output, buffercell=2):
#start_date='2018123122'
#end_date='2019010101'
#shapefile='einzugsgebiet.shp'
#output=True
#buffercell=2

    SystemPath = os.getcwd()
    os.chdir(SystemPath)
    
    start_date = datetime.strptime(start_date,'%Y%m%d%H')
    end_date = datetime.strptime(end_date,'%Y%m%d%H')
    
    dts = [dt.strftime('%Y%m%d') for dt in daterange(start_date.replace(hour=0), end_date.replace(hour=0), relativedelta(days=1))]
    dts_historical = [dt.strftime('%Y%m%d') for dt in daterange(start_date.replace(hour=0), end_date.replace(hour=0), relativedelta(months=1))]
    years = list(range(start_date.year, end_date.year+1))

    # create local directories
    try:
        os.mkdir('zip')
        os.chdir('zip')
        os.mkdir('recent')
        os.mkdir('historical')
        os.chdir('../')
    except OSError:
        pass
    try:    
        os.mkdir('RadolanData')
    except OSError:
        pass
    try:
        os.mkdir('Data')
    except OSError:
        pass

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
    
    ## download archives
    os.chdir('zip')
    # avoid downloading 2018 multiple times
    if (datetime.now().year in years) and (datetime.now().year - 1 in years):
        years = years[:-1]
    
    for year in years:
        if (year == datetime.now().year or year == (datetime.now().year - 1)):
            ftp.cwd('/climate_environment/CDC/grids_germany/hourly/radolan/recent/asc/')
            os.chdir('recent')
            files=ftp.nlst()
            for dt, file in product(dts, files):
                if dt in file:
                    with open(file, 'wb') as f:
                        print('Downloading ' + file)
                        while True  :
                            try:          
                                ftp.retrbinary("RETR " + file, f.write)
                            except:            
                                print('retry connection to server')
                                continue
            
                            break
            os.chdir('../')

        else:
            ftp.cwd('/climate_environment/CDC/grids_germany/hourly/radolan/historical/asc/{}/'.format(year))
            os.chdir('historical')
            files=ftp.nlst()
            for dt, file in product(dts_historical, files):
                if dt[:-2] in file:
                    with open(file, 'wb') as f:
                        print('Downloading ' + file)
                        while True  :
                            try:          
                                ftp.retrbinary("RETR " + file, f.write)
                            except:            
                                print('retry connection to server')
                                continue
            
                            break
            os.chdir('../')

    ftp.quit()
    
    # extract archives
    os.chdir(SystemPath + '/zip')
    
    if os.listdir('historical'):
        os.chdir('historical')
        files = os.listdir()
        print(files)
        for file in files:
            tar = tarfile.open(file)
            tar.extractall()
            tar.close()
            os.remove(file)
            print('Unpacking {}...'.format(file))
        for dt in dts:
            files = [file for file in os.listdir() if dt in file]
            os.chdir(SystemPath + '/RadolanData')
            for file in files:
                tar = tarfile.open(os.path.join(SystemPath,'zip','historical',file))
                tar.extractall()
                tar.close()
                os.remove(os.path.join(SystemPath,'zip','historical',file))
                print('Unpacking {}...'.format(file))
                os.chdir(os.path.join(SystemPath,'zip','historical'))
    os.chdir(SystemPath)
    os.chdir('zip')
    
    if os.listdir('recent'):
        files = os.listdir('recent')
        os.chdir(SystemPath + '/RadolanData')
        for file in files:
            tar = tarfile.open(os.path.join(SystemPath,'zip','recent',file))
            tar.extractall()
            tar.close()
            os.remove(os.path.join(SystemPath,'zip','recent',file))
            print('Unpacking {}...'.format(file))
    
    os.chdir(SystemPath)
    shutil.rmtree('zip')
    
    # get cellsize from Radolan Raster
    os.chdir('RadolanData')
    files = sorted(os.listdir())
    
    # read cellsize from file
    with open(files[0], 'r') as f:
        for line in f:
            if 'cellsize' in line:
                cellsize = int(re.search(r'\d+', line).group())
                break
    
    #delete data outside of date range
    for file in files[:int(start_date.hour)]:
        os.remove(file)
    
    for file in files[-(24-int(end_date.hour)):]:
        os.remove(file)
    #print('Download complete')
    os.chdir('../')
    ## define radolan projection
    # Native projection ("RADOLAN-Projektion"):
    R = 6370040
    # Erdradius (Kugel)
    R = "+a="+str(R)+" +b="+str(R)
    # zu 'R' zusammenfassen -> Kugel
    nat_proj = "+proj=stere +lon_0=10.0 +lat_0=90.0 +lat_ts=60.0 "+R+" +units=m"
    
    ## start action
    # Check the Boundaryfile <--str or not
    gdfbnd = gp.read_file(shapefile)
    src_epsg = gdfbnd.crs.get('init')
    
    # this gives us the flexibility for the clipping
    BoundsClipDs = np.add(gdfbnd.geometry.total_bounds, np.array([-1,-1,1,1])*buffercell*cellsize)
    
    os.chdir('RadolanData')
    files = sorted(glob.glob("*.asc"))
    initDf = False
    
    for file in files:
        # os.system(PyPath+"\\gdalwarp -overwrite -s_srs EPSG:31467 -t_srs EPSG:31467 -of GTiff " +ResultsPath+"\\file " +ResultsPath+file[:-5]+"GK4.tiff")
    
        print('Processing {}...'.format(file))
        if output == True:
            RadolanTif = gdal.Warp(SystemPath + '\\Data\\' + file[:-4] + '.tiff', file, srcSRS=nat_proj, dstSRS=src_epsg, format="Gtiff", srcNodata=-1, outputBounds=BoundsClipDs, outputBoundsSRS=src_epsg, dstNodata=-9999)
        else:
            RadolanTif = gdal.Warp('', file, srcSRS=nat_proj, dstSRS=src_epsg, format="vrt", srcNodata=-1, outputBounds=BoundsClipDs, outputBoundsSRS=src_epsg, dstNodata=-9999)
        ulx, xres, xskew, uly, yskew, yres  = RadolanTif.GetGeoTransform()
    
        if initDf == False:
    
            # ULCorner=src.transform * (0, 0)
            CellBoundsLR = np.linspace(ulx, ulx + (RadolanTif.RasterXSize * xres), RadolanTif.RasterXSize+1)
            CellBoundsUB = np.linspace(uly, uly + (RadolanTif.RasterYSize * yres), RadolanTif.RasterYSize+1)
    
            # create boundaries of each cell
            cellbounds=np.zeros((RadolanTif.RasterXSize*RadolanTif.RasterYSize,6))
            cellbounds[:,(0,3)]=np.array(list(product(CellBoundsLR[:-1],CellBoundsUB[:-1])))
            cellbounds[:,(2,1)]=np.array(list(product(CellBoundsLR[1:],CellBoundsUB[1:])))
            cellbounds[:,(4,5)]=np.array(list(product(range(0,len(CellBoundsLR)-1),range(0,len(CellBoundsUB)-1))))
    
            # Create a pandas dataframe
            df = pd.DataFrame({'left': cellbounds[:, 0], 'bottom': cellbounds[:,1], 'right':cellbounds[:,2], 'top':cellbounds[:,3], 'Index_column':cellbounds[:,4].astype(int), 'Index_row':cellbounds[:,5].astype(int)})
            b = [box(l, b, r, t) for l, b, r, t in zip(cellbounds[:,0], cellbounds[:,1], cellbounds[:,2], cellbounds[:,3])]
            gdf = gp.GeoDataFrame(df, geometry=b)
            df['gridcellarea'] = gdf['geometry'].area
            gdf.crs = gdfbnd.crs
            #get the area of the gridcells
            gridcellarea=df['gridcellarea'][0]
    
        initDf = True
    
        # write precipitation data into gdf dataframe
        precipitation = RadolanTif.ReadAsArray()
        gdf = gdf.sort_values(['Index_row','Index_column'])
        gdf[file[5:11] + file[12:16]] = np.true_divide(precipitation.reshape(-1,1)[:,0], 10)
    
    # replace gridcodes with new ID, IDs are in order of polygons
    gdfbnd['basinID'] = range(1, gdfbnd.shape[0]+1)
    #Get the area of the boundary polygon
    BndPolygonArea=np.vstack((gdfbnd['basinID'].to_numpy(),gdfbnd.area.to_numpy())).T
    # clip radolan data with input shapefile -->most time consuming timestep
    gdfclip = gp.overlay(gdf, gdfbnd, how='intersection', make_valid=True, use_sindex=None)
    
    # get column names for columns that contain precipitation data
    precipitationCols = [column for column in gdfclip.columns if column.isdigit()]
    
    
    #We compute the cellweights of each clipped cell
    cellweights = (gdfclip['geometry'].area/gdfclip['gridcellarea']).to_numpy()
    #get basin ID vector and also the unique values preserving the row
    basinIDs=gdfclip['basinID'].to_numpy()
    basinIDsunique=pd.unique(gdfclip['basinID'])
    #Get numpy array from the dataframe for all precipitation data as well as for the Basin ID and weights, both together does not work :-(
    precipitationcellsextracted=gdfclip.loc[:,precipitationCols[0]:precipitationCols[-1]].to_numpy()*gridcellarea
    #Multiply by cellweights and add basinID
    precipitationcellweighted=precipitationcellsextracted*cellweights[:, None] #http://scipy-lectures.org/advanced/advanced_numpy/#broadcasting
    #add the weights
    precipitationcellweighted=np.hstack((precipitationcellweighted,cellweights.reshape(-1,1)))
    #add basinIDs as a new column 
    precipitationcellweighted=np.hstack((precipitationcellweighted,basinIDs.reshape(-1,1)))
    #Finally we sum each row to get value for each basin and each time step https://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.reduceat.html
    precipitationbasin = np.add.reduceat(precipitationcellweighted, np.r_[0, np.where(np.diff(precipitationcellweighted[:, -1]))[0] + 1], axis=0)
    #Compute the average weight
    precipitationbasin[:,-2]=precipitationbasin[:,-2]*basinIDsunique/precipitationbasin[:,-1]
    #replace sum of unique basin id by actual basin id
    precipitationbasin[:,-1]=basinIDsunique
    #Sort the BndPolygonarea for final printing and add to final numpy array
    precipitationbasin=np.hstack((precipitationbasin,BndPolygonArea[(basinIDsunique.flatten()-1)][:,1].reshape(-1,1)))
    #compute to mm/h
    precipitationbasin[:,:len(precipitationCols)]=(precipitationbasin[:,:len(precipitationCols)].T/precipitationbasin[:,-1]).T

    os.chdir('../Data')
    if output == True:
        #write out the numpy array 
        for basin in precipitationbasin:
            with open('radolanData_{!s}.csv'.format(basin[-2]), 'w', newline='') as fout:
                fout.write('basin ID: {:d}\n'.format(int(basin[-2])))
                fout.write('average_cellweight: {0:.3f}\n'.format(basin[-3]))
                fout.write('basin_area: {0:.3f}\n'.format(basin[-1]))
                fout.write('Time[yymmddhh],rainfall[mm/h]\n')
                np.savetxt(fout,np.column_stack((precipitationCols, np.around(basin[:len(precipitationCols)],decimals=3))), delimiter=",", fmt="%s")

        if len(precipitationCols) < 500:
            #Write out the updated Polygonshape
            #add new information to old dataframe and write it out as shp
            precipitationbasinsorted=precipitationbasin[precipitationbasin[:,-2].argsort()]
            for i in range(0,len(precipitationCols)):
                gdfbnd[precipitationCols[i]]=precipitationbasinsorted[:,i]
            gdfbnd['BasinIDNew']=precipitationbasinsorted[:,-2]
            gdfbnd['average_cellweight']=precipitationbasinsorted[:,-3]
            gdfbnd['basin_area']=precipitationbasinsorted[:,-1]
            gdfbnd.to_file(os.path.join(os.getcwd(),'Data',shapefile[:-4]+'_Updated.shp'))
    shutil.rmtree(os.path.join(SystemPath,'RadolanData'))
    print('\nOutput saved to /Data')
    
    #return the numpy array, last column =area, sc last =basinID, thd last=average weight, 1:fourthlast=precipitation info
    
    return precipitationbasin

if __name__ == "__main__":
    processradolan('2017123020', '2018010220', os.path.join(os.getcwd(),'Examples','einzugsgebiet.shp'), output=True)
