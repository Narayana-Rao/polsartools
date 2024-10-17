import numpy as np
import glob
from osgeo import gdal,osr
import os
import sys
import simplekml
def write_bin_uav(file,wdata,lat,lon,dx,dy):
    
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([lon,dx,0,lat,0,dy])
    outdata.SetProjection("EPSG:4326")
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!


def create_kml_polygon(corner_coords, output_filename):
    # Create a KML object
    kml = simplekml.Kml()

    # Create a polygon from the corner coordinates
    pol = kml.newpolygon(name="Polygon")
    
    # Close the polygon by repeating the first point at the end
    polygon_coords = corner_coords + [corner_coords[0]]
    
    # Assign the polygon coordinates (in (longitude, latitude) format)
    pol.outerboundaryis.coords = polygon_coords

    # Set a basic style for the polygon (optional)
    pol.style.polystyle.color = simplekml.Color.changealphaint(150, simplekml.Color.red)  # Set transparency and color
    pol.style.polystyle.outline = 1  # Outline on
    
    # Save the KML file
    kml.save(output_filename)
    # print(f"KML file saved as {output_filename}")
def create_extent(inFolder):
    annFile = open(glob.glob(inFolder+'/*.ann')[0], 'r')
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
    print(corner_coordinates)

    output_file = os.path.join(inFolder,"scene_extent.kml")
    create_kml_polygon(corner_coordinates, output_file)
    
    
def uavsar_grd_c3(inFolder):
    create_extent(inFolder)
    annFile = open(glob.glob(inFolder+'/*.ann')[0], 'r')
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

def uavsar_mlc_c3(inFolder):
    create_extent(inFolder)
    annFile = open(glob.glob(inFolder+'/*.ann')[0], 'r')
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
