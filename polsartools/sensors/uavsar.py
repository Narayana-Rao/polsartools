import numpy as np
import glob
from osgeo import gdal,osr
import os
import sys
import simplekml
from polsartools.utils.utils import time_it
def write_bin_uav_old(file,wdata,lat,lon,dx,dy):
    
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([lon,dx,0,lat,0,dy])
    outdata.SetProjection("EPSG:4326")
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!

def write_bin_uav(file, wdata, lat, lon, dx, dy, sensor_type="UAVSAR"):
    [rows, cols] = wdata.shape
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, cols, rows, 1, gdal.GDT_Float32)
    if dy<0:
        dy = -dy
    outdata.SetGeoTransform([lon, dx, 0, lat, 0, dy])
    outdata.SetProjection("EPSG:4326")
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache()
    
    # Write the header file
    header_filename = file.replace('.bin', '.hdr') 
    with open(header_filename, 'w') as header_file:
        header_file.write("ENVI\n")
        header_file.write(f"description = {{{file}}}\n")
        header_file.write(f"samples = {cols}\n")
        header_file.write(f"lines = {rows}\n")
        header_file.write("bands = 1\n")
        header_file.write("header offset = 0\n")
        header_file.write("file type = ENVI Standard\n")
        header_file.write("data type = 4\n")  # Float32
        header_file.write("interleave = bsq\n")
        # header_file.write("wavelength units = METERS\n")
        header_file.write("byte order = 0\n")
        header_file.write(f"sensor type = {sensor_type}\n")
        map_info = f"map info = {{Geographic Lat/Lon, 1, 1, {lon}, {lat}, {dx}, {dy}, WGS-84}}\n"
        header_file.write(map_info)
        # Write the projection information
        # header_file.write("projection = GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS_84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.017453292519943295],AUTHORITY[\"EPSG\",\"4326\"]]\n")
        header_file.write(f"band names = {{{os.path.basename(file)}}}\n")

def create_kml_polygon(corner_coords, output_filename):
    
    kml = simplekml.Kml()
    pol = kml.newpolygon(name="Polygon")
    polygon_coords = corner_coords + [corner_coords[0]]
    pol.outerboundaryis.coords = polygon_coords
    pol.style.polystyle.color = simplekml.Color.changealphaint(150, simplekml.Color.red)  
    kml.save(output_filename)

def create_extent(annFile):
    inFolder = os.path.dirname(annFile)
    annFile = open(annFile, 'r')
    for line in annFile:
        if "Approximate Upper Left Latitude" in line:
            uly = float(line.split('=')[1].split(';')[0])
        if "Approximate Upper Left Longitude" in line:
            ulx = float(line.split('=')[1].split(';')[0])
        if "Approximate Lower Right Latitude" in line:
            lry = float(line.split('=')[1].split(';')[0])
        if "Approximate Lower Right Longitude" in line:
            lrx = float(line.split('=')[1].split(';')[0])
        if "Approximate Lower Left Latitude" in line:
            lly = float(line.split('=')[1].split(';')[0])
        if "Approximate Lower Left Longitude" in line:
            llx = float(line.split('=')[1].split(';')[0]) 
        if "Approximate Upper Right Longitude" in line:
            urx = float(line.split('=')[1].split(';')[0])
        if "Approximate Upper Right Latitude" in line:
            ury = float(line.split('=')[1].split(';')[0]) 

    corner_coordinates = [
        (ulx, uly),  # (longitude, latitude) coordinate 1
        (urx, ury),  # coordinate 2
        (lrx, lry),  # coordinate 3
        (llx, lly),  # coordinate 4
    ]
    # print(corner_coordinates)

    output_file = os.path.join(inFolder,"scene_extent.kml")
    create_kml_polygon(corner_coordinates, output_file)
    
@time_it    
def uavsar_grd(annFile):
    inFolder = os.path.dirname(annFile)
    create_extent(annFile)
    annFile = open(annFile, 'r')
    for line in annFile:
        if "grd_mag.set_rows" in line:
            rows = int(line.split('=')[1].split(';')[0])
        if "grd_mag.set_cols" in line:
            cols = int(line.split('=')[1].split(';')[0])
        if "grd_mag.row_addr" in line:
            lat = float(line.split('=')[1].split(';')[0])
        if "grd_mag.col_addr" in line:
            lon = float(line.split('=')[1].split(';')[0])
        if "grd_mag.row_mult" in line and "(deg/pixel)" in line:
            dy = float(line.split('=')[1].split(';')[0])
        if "grd_mag.col_mult" in line and "(deg/pixel)" in line:
            dx = float(line.split('=')[1].split(';')[0])

    outFolder = inFolder+'/C3'
    if not os.path.isdir(outFolder):
        os.mkdir(outFolder)

    hhhh = np.fromfile(glob.glob(inFolder+'/*HHHH*.grd')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C11.bin',hhhh,lat,lon,dx,dy)
    del hhhh
    vvvv = np.fromfile(glob.glob(inFolder+'/*VVVV*.grd')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C33.bin',vvvv,lat,lon,dx,dy)
    del vvvv
    hvhv = np.fromfile(glob.glob(inFolder+'/*HVHV*.grd')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C22.bin',2*hvhv,lat,lon,dx,dy)
    del hvhv
    hhhv = np.fromfile(glob.glob(inFolder+'/*HHHV*.grd')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C12_real.bin',np.real(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C12_imag.bin',np.imag(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    del hhhv
    hhvv = np.fromfile(glob.glob(inFolder+'/*HHVV*.grd')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C13_real.bin',np.real(hhvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C13_imag.bin',np.imag(hhvv),lat,lon,dx,dy)
    del hhvv
    hvvv = np.fromfile(glob.glob(inFolder+'/*HVVV*.grd')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C23_real.bin',np.real(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C23_imag.bin',np.imag(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    del hvvv

    file = open(outFolder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()  
    print("Extracted C3 files to %s"%outFolder)
@time_it  
def uavsar_mlc(annFile):
    create_extent(annFile)
    inFolder = os.path.dirname(annFile)
    create_extent(annFile)
    annFile = open(annFile, 'r')
    for line in annFile:
        if "mlc_mag.set_rows" in line:
            rows = int(line.split('=')[1].split(';')[0])
        if "mlc_mag.set_cols" in line:
            cols = int(line.split('=')[1].split(';')[0])
        if "grd_mag.row_addr" in line:
            lat = float(line.split('=')[1].split(';')[0])
        if "grd_mag.col_addr" in line:
            lon = float(line.split('=')[1].split(';')[0])
        if "grd_mag.row_mult" in line and "(deg/pixel)" in line:
            dy = float(line.split('=')[1].split(';')[0])
        if "grd_mag.col_mult" in line and "(deg/pixel)" in line:
            dx = float(line.split('=')[1].split(';')[0])
        # if "set_plon" in line  and "(deg)" in line:
        #     lon = float(line.split('=')[1].split(';')[0])
        # if "set_plat" in line  and "(deg)" in line:
        #     lat = float(line.split('=')[1].split(';')[0])
        # if "mlc_mag.row_mult" in line and "(m/pixel) " in line:
        #     dy = float(line.split('=')[1].split(';')[0])*0.00001
        # if "mlc_mag.col_mult" in line and "(m/pixel) " in line:
        #     dx = float(line.split('=')[1].split(';')[0])*0.00001

    outFolder = inFolder+'/C3'
    if not os.path.isdir(outFolder):
        os.mkdir(outFolder)

        
    hhhh = np.fromfile(glob.glob(inFolder+'/*HHHH*.mlc')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C11.bin',hhhh,lat,lon,dx,dy)

    del hhhh
    vvvv = np.fromfile(glob.glob(inFolder+'/*VVVV*.mlc')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C33.bin',vvvv,lat,lon,dx,dy)
    del vvvv
    hvhv = np.fromfile(glob.glob(inFolder+'/*HVHV*.mlc')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C22.bin',2*hvhv,lat,lon,dx,dy)
    del hvhv
    hhhv = np.fromfile(glob.glob(inFolder+'/*HHHV*.mlc')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C12_real.bin',np.real(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C12_imag.bin',np.imag(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    del hhhv
    hhvv = np.fromfile(glob.glob(inFolder+'/*HHVV*.mlc')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C13_real.bin',np.real(hhvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C13_imag.bin',np.imag(hhvv),lat,lon,dx,dy)
    del hhvv
    hvvv = np.fromfile(glob.glob(inFolder+'/*HVVV*.mlc')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C23_real.bin',np.real(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C23_imag.bin',np.imag(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    del hvvv

    file = open(outFolder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()   
    print("Extracted C3 files to %s"%outFolder)
