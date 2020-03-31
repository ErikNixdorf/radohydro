#!/usr/bin/env python3.7
"""This script intersects a given shapefile of one or multiple basins with the Radolan ASCII data by DWD and returns a
   .csv with the rainfall intensity for each basin in the given period as output. Optionally the clipped shape file
   and the geoTIFFs can be requested
   Major Changes in V0.3:
       -Performance by vectorization
   Major Changes in V04:
       -Performance (factor 3) by having fully streambased solution
       -Script splitted into functions to increase flexibility
	Minor in V041: Fixed reprojection bugs, allows no shape inpt, improved memory usage
   Major Changes in V05:allows download from regnie data as well
"""

__author__ = "Erik Nixdorf, Marco Hannemann"
__propertyof__ = "Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ. "
__email__ = "erik.nixdorf@ufz.de, marco.hannemann@ufz.de"
__version__ = "0.5"
import time
from ftplib import FTP
import tarfile
from datetime import datetime
from dateutil.relativedelta import relativedelta
from itertools import product
from pyproj import CRS
import re
import numpy as np
from io import BytesIO
from rasterio.io import MemoryFile
import rasterio
import sys
import gzip
import time
# import required geometrytools from geotools
import geotools.geotools as gs


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
    Defines the native radoproj as a projection CRS    
    """
    # Native projection ("RADOLAN-Projektion"):
    R = earthradius
    # Erdradius (Kugel)
    R = "+a=" + str(R) + " +b=" + str(R)
    # zu 'R' zusammenfassen -> Kugel
    nat_proj = rasterio.crs.CRS.from_proj4(
        "+proj=stere +lon_0=10.0 +lat_0=90.0 +lat_ts=60.0 " + R + " +units=m")
    return nat_proj


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


def connect_ftp(server = 'opendata.dwd.de',connected=False):
    while not connected:
        try:
            ftp = FTP(server)
            ftp.login()
            connected = True
    
        except:
            time.sleep(5)
            print('Reconnect to Server')
            pass
    return ftp


def regnierow_to_array(line,no_per_line=611):
    """
    A function which processes linestring from regnie to np array
    Fixed a stupid bug in regnie that if daily precip is larger 1000 cm, no space
    is between values and hence can't be split easisly by numpy
    assumptions: 1)never will be a daily rainfall >1999 cm
    2)nearby cells of rainfall >1000 will at least have 100 cm'    
    """
    #subfunction fo find all occurences of substring in string
    def findOccurrences(s, ch):
        return [i for i, letter in enumerate(s) if letter == ch]

    #%% separate all numbers to one space
    line=str(line).replace('-',' -').strip()
    line=line.replace('-',' -').strip()
    line=line.replace('   ',' ').strip()
    line=line.replace('  ',' ').strip()
    line=line.replace(' ',' ').strip()
    value_list=line.split()
    #%% process to numpy array
    if len(value_list)==no_per_line:
        line_array=np.array(value_list)
    else:
        print('local rain >100mm detected, repair linestring')
        #we have to correct digits by looping over the list
        value_list_updated=list()
        for value in value_list:
            if len(value)>4:
                #if value can be divided by four most easy
                if np.mod(len(value),4)==0:
                    value_list_updated.extend([value[i:i+4] for i in range(0, len(value), 4)])
                else:
                    #we if the first value is either minus or 1-3, the three digit
                    #rainfall is in the last three digits
                    if value[0] in ['-','1','3','4']:
                        value_list_updated.append(value[-3:])
                        value=value[:-3]
                        value_list_updated.extend([value[i:i+4] for i in range(0, len(value), 4)])
                    else:
                        #in the other case the first value has three digits
                        value_list_updated.append(value[:3])
                        value=value[3:]
                        value_list_updated.extend([value[i:i+4] for i in range(0, len(value), 4)])
            else:
                value_list_updated.append(value) 
        line_array=np.array(value_list_updated)
        if len(value_list_updated)!=no_per_line:
            print('still problem, improve correction function')

    return line_array             


# The Actual downloader for Radolan to archive
def rado_io(start_date='20171231',
            end_date='20180101',
            shapefile='.\Examples\einzugsgebiet.shp', number_frmt=np.int16):
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
    server='opendata.dwd.de'
    ftp=connect_ftp(server = server,connected = False)  

    # Set Initial run to true
    initDf = True
    # avoid downloading 2018 multiple times
    if (datetime.now().year in years) and (datetime.now().year - 1 in years):
        years = years[:-1]

    for year in years:
        #check whether data is recent
        if year == datetime.now().year:
            ftp.cwd(
                '/climate_environment/CDC/grids_germany/hourly/radolan/recent/asc/'
            )
            files = ftp.nlst()
            for dt, file in product(dts, files):
                if dt in file:
                    print('Retrieving {}...'.format(file))
                    retrieved = False
                    archive = BytesIO()
                    # try to retrieve file
                    while not retrieved:
                        try:
                            ftp.retrbinary("RETR " + file, archive.write)
                            retrieved = True
                        except:
                            print('reconnect to ftp')
                            ftp = FTP(server)
                            ftp.login()
                            ftp.cwd(
                                '/climate_environment/CDC/grids_germany/hourly/radolan/recent/asc/'
                            )

                    archive.seek(0)
                    archive_daily = tarfile.open(fileobj=archive)
                    #extract file to bytestream
                    for member in archive_daily.getmembers():
                        radolan_io = archive_daily.extractfile(member.name)
                        with MemoryFile(radolan_io) as memfile:
                            with memfile.open() as rado_ds:
                                rado_data = rado_ds.read()[0].astype(number_frmt)
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
                                    # if a shapefile exist
                                    if shapefile is not None:
                                        rado_clip_data, rado_clip_transform, cols, rows = gs.buffered_raster_clipping(
                                            rado_data,
                                            shape_inpt=shapefile,
                                            raster_transfrm=rado_transform,
                                            raster_proj=rado_proj)
                                    else:
                                        rado_clip_data = rado_data
                                        rado_clip_transform = rado_transform
                                        rows = [rado_data.shape[0], 0]
                                        cols = [0, rado_data.shape[1]]
                                    #generate the footprint cells
                                    radocells = gs.create_footprint_cells(
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
                                    try:
                                        rado_stacked_data = np.dstack(
                                            (rado_stacked_data,
                                             rado_clip_data))
                                    except Exception as e:
                                        print(e)
                                        sys.exit(
                                            'Memory Error :-(, buy more RAM')
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
                    retrieved = False
                    archive = BytesIO()
                    # try to retrieve file
                    while not retrieved:
                        try:
                            ftp.retrbinary("RETR " + file, archive.write)
                            retrieved = True
                        except:
                            print('reconnect to ftp')
                            ftp = FTP(server)
                            ftp.login()
                            ftp.cwd(
                                '/climate_environment/CDC/grids_germany/hourly/radolan/historical/asc/{}/'.
                                format(year))
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
                            print('extract daily file', members_daily)
                            for member in archive_daily.getmembers():
                                radolan_io = archive_daily.extractfile(
                                    member.name)
                                with MemoryFile(radolan_io) as memfile:
                                    with memfile.open() as rado_ds:
                                        rado_data = rado_ds.read()[0].astype(
                                            number_frmt)
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
                                            # if a shapefile exist
                                            if shapefile is not None:
                                                rado_clip_data, rado_clip_transform, cols, rows = gs.buffered_raster_clipping(
                                                    rado_data,
                                                    shape_inpt=shapefile,
                                                    raster_transfrm=
                                                    rado_transform,
                                                    raster_proj=rado_proj)
                                            else:
                                                rado_clip_data = rado_data
                                                rado_clip_transform = rado_transform
                                                rows = [rado_data.shape[0], 0]
                                                cols = [0, rado_data.shape[1]]

                                            #generate the footprint cells
                                            radocells = gs.create_footprint_cells(
                                                transform=rado_clip_transform,
                                                data_size=rado_clip_data.shape,
                                                proj_crs=rado_proj,
                                                rado_divisor=1000)
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
                                            try:
                                                rado_stacked_data = np.dstack(
                                                    (rado_stacked_data,
                                                     rado_clip_data))
                                            except Exception as e:
                                                print(e)
                                                sys.exit(
                                                    'Memory Error :-(, buy more RAM'
                                                )
                                            rado_dates.append(
                                                radoname_to_date(
                                                    member.name, 'minutes'))
                    rado_dates = sorted(rado_dates)
                    print('Processing {}...finished'.format(file))
    try:
        ftp.quit()
    except Exception as e:
        print(e)
    # repar the radocell crs
    #define the radolan_projection
    rado_proj_string='+proj=stere +lat_0=90 +lat_ts=90 +lon_0=10 +k=0.93301270189 + x_0=0 +y_0=0 +a=6370040 +b=6370040 +to_meter=1000 +no_defs'
    radocells.crs = CRS(rado_proj_string)
    return rado_stacked_data, rado_dates, radocells


def regnie_io(start_date='20171231',
            end_date='20180101',
            shapefile='.\Examples\einzugsgebiet.shp', number_frmt=np.float16):
    """
    process data from daily regnie dataset
    """
    # create the dates
    start_date = datetime.strptime(start_date, '%Y%m%d')
    end_date = datetime.strptime(end_date, '%Y%m%d')
    
    years = list(range(start_date.year, end_date.year + 1))
    
    regnie_proj=CRS('epsg:4326')
    #regnie_transform=(55.083333,1/60,0,5.833333,0,1/120)
    regnie_transform=(5.833333-1/120,+1/60,0,55.083333+1/240,0,-1/120)
    regnie_nan = -999
    # Connect to the Server
    server='opendata.dwd.de'
    ftp=connect_ftp(server = server,connected = False)
    
    # Set Initial run to true
    initDf = True
    
    for year in years:
            try:
                ftp.cwd('/climate_environment/CDC/grids_germany/daily/regnie')
            except:
                print('reconnect to ftp')
                server='opendata.dwd.de'
                ftp=connect_ftp(server = server,connected = False) 
                ftp.cwd('/climate_environment/CDC/grids_germany/daily/regnie')
            files = ftp.nlst()
            matching_file = [s for s in files if 'ra'+str(year)+'m' in s][0]
            #retrieve file
            print('Retrieving {}...'.format(matching_file))
            retrieved = False
            archive = BytesIO()
            # try to retrieve file
            while not retrieved:
                try:
                    ftp.retrbinary("RETR " + matching_file, archive.write)
                    retrieved = True
                except:
                    print('reconnect to ftp')
                    ftp = FTP(server)
                    ftp.login()
                    ftp.cwd(
                        '/climate_environment/CDC/grids_germany/hourly/regnielan/historical/asc/{}/'.
                        format(year))
            archive.seek(0)
            archive_yearly = tarfile.open(fileobj=archive)
            #double_zipped so we need to get daily archives
            for members_daily in archive_yearly.getnames():
                # check whether the day is in our time span
                members_date = datetime.strptime(members_daily[2:8],'%y%m%d')
                if members_date >= start_date and members_date <= end_date:                
                    tar_daily = archive_yearly.extractfile(
                        members_daily)
                    tar_daily.seek(0)
                    #get each row as a list
                    regnie_io = gzip.open(tar_daily).read().decode("utf-8").splitlines()
                    #last line is not necessary
                    regnie_io=regnie_io[:-1]
                    #convert to numpy array
                    regnie_data=np.array([]).reshape(0,611).astype(number_frmt)
                    for line in regnie_io:
                        line_array=regnierow_to_array(line,no_per_line=611).astype(number_frmt)
                        regnie_data = np.vstack([regnie_data,line_array])
                    regnie_data=np.where(regnie_data==regnie_nan, np.nan, regnie_data)
                    #write the information for transform
                    if initDf:
                        # do the complicated buffer clipping
                        # if a shapefile exist
                        if shapefile is not None:
                            regnie_clip_data, regnie_clip_transform, cols, rows = gs.buffered_raster_clipping(
                                regnie_data,
                                shape_inpt=shapefile,
                                raster_transfrm=
                                regnie_transform,
                                raster_proj=regnie_proj)
                        else:
                            regnie_clip_data = regnie_data
                            regnie_clip_transform = regnie_transform
                            rows = [regnie_data.shape[0], 0]
                            cols = [0, regnie_data.shape[1]]
    
                        #generate the footprint cells
                        regniecells = gs.create_footprint_cells(
                            transform=regnie_clip_transform,
                            data_size=regnie_clip_data.shape,
                            proj_crs=regnie_proj,
                                rado_divisor=1)
                        #initialize the merged dataset
                        regnie_stacked_data = regnie_clip_data
                        #the dates
                        regnie_dates = [members_date]
                        initDf = False
                        # if we initialised already, computation is easy
                    else:
                        regnie_clip_data = regnie_data[rows[
                            1]:rows[0], cols[0]:cols[1]]
                        try:
                            regnie_stacked_data = np.dstack(
                                (regnie_stacked_data,
                                 regnie_clip_data))
                        except Exception as e:
                            print(e)
                            sys.exit(
                                'Memory Error :-(, buy more RAM'
                            )
                        regnie_dates.append(members_date)
                    print('processing for day',members_date.date(), 'done')
            regnie_dates = sorted(regnie_dates)
            print('Processing {}...finished'.format(matching_file))
    try:
        ftp.quit()
    except Exception as e:
        print(e)
    #return values
    return regnie_stacked_data, regnie_dates, regniecells
    
    
#define a main function which calls all other subfunctions
def radohydro(start_date='20190201',
              end_date='20200301',
              shape_inpt='.\\Examples\\Mueglitz_Basin.shp',datasource='regnie',
              shape_integration=True,
              outpt_proj='epsg:25833',
              Output=True,
              outpt_nm='precip_Mueglitz',number_frmt=np.int16):
    """
    Couples the four main function for the entire workflow
    #'.\Examples\Mueglitz_Basin.shp'
    """
    # retrieve data
    if datasource=='radolan':
        rado_stacked_data, rado_dates, radocellgrid = rado_io(
            start_date=start_date, end_date=end_date, shapefile=shape_inpt,number_frmt=number_frmt)
    if datasource=='regnie':
        number_frmt=np.float64
        rado_stacked_data, rado_dates, radocellgrid = regnie_io(
            start_date=start_date, end_date=end_date, shapefile=shape_inpt,number_frmt=number_frmt)
    #map data on vectorgrid
    radocell_precip, precip_col_nms = gs.map_arraystack_on_cellgrd(
        rado_stacked_data,
        rado_dates,
        radocellgrid,
        numerator=10,
        Output=Output,
        outpt_proj=outpt_proj,number_frmt=number_frmt)
    # delete stacked numpy array
    del rado_stacked_data
    # Integrate values to boundary shape, if needed
    if shape_integration:
        radocell_clip, gdf_shape = gs.map_cellgrd_on_polyg(
            radocell_precip, shape_inpt=shape_inpt, outpt_proj=outpt_proj)
        
        gs.compute_polyg_values(
            radocell_clip,
            gdf_shape,header='rainfall[mm]',
            datacol_type='AllDigits',
            Output=Output,
            outpt_proj=outpt_proj,
            outpt_nm=outpt_nm)
    else:
        print('No integration of data on geometry requested')
    return None


#%% Test our functions
if __name__ == "__main__":
    radohydro()
