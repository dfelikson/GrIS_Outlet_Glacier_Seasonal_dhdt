#!/usr/bin/env python

from osgeo import gdal
from osgeo import gdalconst

import numpy
from bilinear_interpolate import *

def interp1985DEM(x,y): #,method='bilinear'):
    dem_filename = './Cheat_matfiles/aerodem_1985_utm22_1_epsg3413.lzw.tif'
    
    raster = gdal.Open(dem_filename, gdalconst.GA_ReadOnly)
    rasterBand = raster.GetRasterBand(1)
    rasterBandArray = rasterBand.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize)
    
    # No data values
    rasterBandArray[rasterBandArray == -9999] = np.nan
    rasterBandArray[rasterBandArray < -3e+38] = np.nan
    rasterBandArray[rasterBandArray == 0] = np.nan
    rasterBandArray[np.abs(rasterBandArray - 0) < np.finfo(float).eps] = numpy.nan

    gt = raster.GetGeoTransform()
    
    z = np.nan * np.ones(x.shape)
    for i in range(len(x)):
        xi = x[i]
        yi = y[i]
        imagex = (xi - gt[0]) / gt[1]
        imagey = (yi - gt[3]) / gt[5]

        if (imagex >= 0) and (imagey >= 0) and (imagex < raster.RasterXSize) and (imagey < raster.RasterYSize):
            zi = bilinear_interpolate(rasterBandArray, imagex, imagey, nodataValue=np.nan)
        
        z[i] = zi
    
    # Not sure why there's any data being returned for these nodes,
    # but this is a hardcoded way to cleanup bad data
    z[x > -140000] = np.nan
    
    return z


def surface_shift(md, surface1, surface2, method='linear_fit'):
    # Find the surface elevation differences
    pos = md.mask.ice_levelset < 0
    surface_diffs = surface2[pos]-surface1[pos]
    
    # Combine the two surfaces
    surf_shifted = np.nan * np.ones( md.mesh.numberofvertices )
    pos = ~np.isnan(surface2)
    surf_shifted[pos] = surface2[pos]
    
    if method == 'mean_of_diffs':
        surf_shift = np.nanmean(surface2[pos]-surface1[pos])
        
        # Combine the two surfaces
        pos = np.isnan(surface2)
        surf_shifted[pos] = surface1[pos] - surf_shift

    elif method == 'linear_fit':
        pos = (~np.isnan(surface2 - surface1)) & (surface2 > 700) & (surface2 < 1000)
        p = np.polyfit(surface1[pos], surface2[pos], 1)
                            
        # Combine the two surfaces
        pos = np.isnan(surface2)
        surf_shifted[pos] = np.polyval(p, surface1[pos])

    return surf_shifted
