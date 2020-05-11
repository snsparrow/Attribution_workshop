#!/usr/bin/env python
#############################################################################
# Program : create_mask_txt_shp.py
# Author  : Sihan Li
# Date    : 06/04/2019
# Purpose : Create a netCDF mask for a shapefile
# Updates : Updated 11/05/2020 Sarah Sparrow to include command line argument parsing
#############################################################################
import numpy as np
import argparse
from matplotlib.path import Path
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import unicodedata

from osgeo import ogr

###################################################################################

def load_grid(fname,latname='latitude01',lonname='longitude01'):
	# load region grid and returns list of points (lon,lat)
	with Dataset(fname,'r') as f:
		# Load 1d arrays of lon and lat
		lat=f.variables[latname][:]
		lon=f.variables[lonname][:]
		
		if len(lat.shape)==2:
			# 2D lat and lon:
			lonxx=lon
			latyy=lat
		else:
			# Create 2D arrays of lon and lat
			lonxx,latyy=np.meshgrid(lon,lat)
	return lonxx,latyy

##################################################################################

# Function to create mask given polygons object and points array
def create_mask(polygons,points,nlat,nlon):
	# Convert polygons to mask (true if inside the polygon region)
	# add the masks for multiple polygons together
	for i,polygon in enumerate(polygons):
		# Determine if  points inside polygon
		tmp_mask = polygon.contains_points(points)
		# Reshape mask to dimensions of the grid
		tmp_mask=np.reshape(tmp_mask,[nlat,nlon])
		try:
			mask=tmp_mask | mask
		except:
			mask=tmp_mask

	return ~mask # Invert the mask so true is outside the region

#################################################################################

def load_shapefile(shapefile,fieldname,field_list=None):
	
	print('Loading Shapefile ',shapefile)
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(shapefile, 0)
	layer = dataSource.GetLayer()
	counties={}
	boundaries=[]
	types=[]

	layerDefinition = layer.GetLayerDefn()


	for i in range(layerDefinition.GetFieldCount()):
    		print(layerDefinition.GetFieldDefn(i).GetName())
	print("\n")
	for feature in layer:
		try:
			region=feature.GetField(fieldname)
		except ValueError:
			print(feature.items())
			raise Exception('Error, field "'+fieldname+'" does not exist in the shapefile')
		print(region)
		if field_list is not None and region not in field_list:
			# Skip this region
			continue


		geometry=feature.GetGeometryRef()
		boundary=eval(feature.geometry().Boundary().ExportToJson())
		#geometry = json['geometry']

		if boundary['type']=='LineString':
			polygons=[Path(np.array(boundary['coordinates']))]
		elif boundary['type']=='MultiLineString':
			polygons=[]
			for p in boundary['coordinates']:
				polygons.append(Path(np.array(p)))
		else:
			print('Error: unknown geometry')
			continue

		if region is not None and boundary is not None:
			counties[region]=polygons
	
	return counties

################################################################################
		
def create_netcdf(template,data,outname,template_var='pr'):
	# create outfile object
	outfile=Dataset(outname,'w')

	# Create dimensions copied from template file
	temp=template.variables[template_var]
	for dim in temp.dimensions:
		if dim[:3]=='lat' or dim[:3] =='lon':
			leng=len(template.dimensions[dim])
		
			outfile.createDimension(dim[:3],leng)
			outfile.createVariable(dim[:3],'f',(dim[:3],))
			outfile.variables[dim[:3]][:]=template.variables[dim][:]
			#print template.variables[dim].__dict__
			for att in template.variables[dim].ncattrs():
				outfile.variables[dim[:3]].__setattr__(att,template.variables[dim].__getattribute__(att))

	# Create data variable (named region_mask)
	outfile.createVariable('region_mask','f',['lat','lon'])
	outfile.variables['region_mask'][:]=(data-1)*-1
	
	#outfile.flush()
	outfile.close()

#
###############################################################################

# Create a number of masks, from the shapefile, for a specific grid
# Area outside the polygons is True/1, area inside the polygons is False/0
#
# 
# f_grid: filename of netcdf file contatining grid information
# latname, lonname, template_var: variable names for latitude, longitude and a template variable in f_grid netcdf file
# shapefile: path of shapefile containing polygons
# fieldname: attribute name in shapefile used to identify each field. 
# field_list: (optional)- specify a list (subset) of fields to create masks for, otherwise masks will be created covering all fields
# plot, netcdf_out: (optional) booleans- whether or not to create output plot and/or output netcdf file
#
# Returns: dictionary of region_name:mask_array pairs
#
def create_masks(f_grid, shapefile, fieldname, field_list=None, latname='lat',lonname='lon',template_var='pr', plot=False, netcdf_out = False):

	# first create folder (if needed)
	if netcdf_out and not os.path.exists('masks_netcdf'):
		os.mkdir('masks_netcdf')

	# Load Shape file
	regions=load_shapefile(shapefile,fieldname,field_list=field_list)
	
	# Load lat lon grid (for mask)
	lonxx,latyy=load_grid(f_grid,latname=latname,lonname=lonname)
	nlat,nlon=lonxx.shape
	# Update lon to be from -180 to 180 
	# NOTE: (this is only if the shapefile uses lat coordinates from -180-180 )
	# Comment out otherwise
	lonxx[lonxx>180]=lonxx[lonxx>180]-360
	# Turn lat and lon into a list of coordinates
	points = np.vstack((lonxx.flatten(),latyy.flatten())).T

	# Either loop over all regions, or list of regions specified by 'field_list'
	if field_list == None:
		field_list = iter(regions.keys())

	# Dictionary of masks
	masks={}
	
	# Do the loop
	print('Looping over regions and creating gridded masks')
	for region in field_list:
		region_ascii = unicodedata.normalize('NFKD',str(region)).encode('ascii','ignore')

		print(fieldname,'=',region)
		polygons = regions[region]
		# Create mask out of polygon, matching points from grid
		mask = create_mask(polygons,points,nlat,nlon)
		
		# Add to dictionary
		masks[region_ascii] = mask

		if netcdf_out:
			create_netcdf(Dataset(f_grid,'r'),mask,'masks_netcdf/mask_'+region+'.nc', template_var=template_var)

		if plot:
			plot_mask(lonxx,latyy, mask,region,plot)


	return masks

def plot_mask(lons,lats, mask,region,plot):
	# first create folder (if needed)
	if plot and not os.path.exists('plots'):
		os.mkdir('plots')
	
	plt.clf()	
	ax = plt.axes(projection=ccrs.PlateCarree())
	ax.coastlines()	
	ax.add_feature(cfeature.BORDERS)

	cm=ax.contourf(lons, lats, mask) 
		
	plt.title('Mask: '+region)
	plt.savefig('plots/mask_'+region+'.png')


def main():
    #Read in the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--f_grid", help="The netcdf file (including path) that acts as a template file for the desired grid of the mask.")
    parser.add_argument("--shapefile", help="The shapefile (including path) to convert to a netcdf template.")
    parser.add_argument("--shapefield", help="The shapefile field to use define layers")
    parser.add_argument("--template_var",  help="The template variable")
    parser.add_argument("--latname", default="lat", help="The latitude name in the template file, default = 'lat'")
    parser.add_argument("--lonname", default="lon", help="The longitude name in the template file, default = 'lon'")
    parser.add_argument("--field_list", default=None, nargs='+', help="Optionally specify a list (subset) of fields to create masks for, otherwise masks will be created covering all fields")
    args = parser.parse_args()

    print(args.field_list)
    # Create mask for each layer in shapefile
    create_masks(args.f_grid,args.shapefile,args.shapefield,field_list=args.field_list,latname=args.latname,lonname=args.lonname,template_var=args.template_var,plot=True, netcdf_out=True)

#Washerboard function that allows main() to run on running this file
if __name__=='__main__':
    main()
