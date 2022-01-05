#
# Prep_DOH365.py : prepare DOH TM for construction
#
import pandas as pd
import geopandas as gpd
import xlrd
import pyproj
import numpy as np
import fiona
fiona.supported_drivers['KML'] = 'rw'
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def PPM( k ):
    return (k-1.)*1000.*1000

P4_TM  = '+proj=tmerc +lat_0={lat_0} +lon_0={lon_0} +k_0={k_0}'\
        ' +x_0={x_0} +y_0={y_0} +a={a} +b={b} +units=m +no_defs'
#################################################################
FILE=r'DOH/GPS 107.xlsx'
df = pd.read_excel( FILE, engine='openpyxl', skiprows=4, header=None ) 
TM_k0 = df[9].iloc[9]
TM_FE = df[9].iloc[3]
TM_FN = df[9].iloc[7]
df.drop( labels=[9, 10,11], axis=1, inplace=True )
df.columns=['Stn','WGS_Lat','WGS_Lng','hae', 'UTM_N','UTM_E', 'MSL', 
               'TM_N', 'TM_E']
#################################################################
W84 = pyproj.CRS(4326)  # WGS84
UTM = pyproj.CRS(32647)  # UTM zone47
TM_PARAM = { 'lat_0': 0.0, 'lon_0': 99.0,   # LOCAL central meridian
             'k_0'  : TM_k0, 'x_0'  : TM_FE, 'y_0'  : TM_FN,
             'a'    : pyproj.Geod(ellps='WGS84').a, 
             'b'    : pyproj.Geod(ellps='WGS84').b, }
P4_TM = P4_TM.format( **TM_PARAM )
######################################
#P4_TM = 'proj=lcc +lon_0=0.0 +lat_1=19.55 +lat_0=19.55 '\
#         '+k_0=1.0000847 +x_0=-9800000 +y_0=-2875000 '\
#         '+a=6378137 +b=6356752.3142 +units=m +no_defs '
######################################
TMC = pyproj.CRS.from_proj4( P4_TM )
TR_UTM_TMC = pyproj.Transformer.from_crs( UTM, TMC )
TR_UTM_W84 = pyproj.Transformer.from_crs( UTM, W84 )

tm = TR_UTM_TMC.transform( df.UTM_E, df.UTM_N )
df['P4TM_E'], df['P4TM_N'] = tm[0], tm[1]

ll = TR_UTM_W84.transform( df.UTM_E, df.UTM_N )
df['P4lng'], df['P4lat'] = ll[1], ll[0]

#################################################################
df['diffE'] = df['TM_E']-df['P4TM_E']
df['diffN'] = df['TM_N']-df['P4TM_N']
print( df[['MSL','diffE','diffN']].describe() )

#################################################################
p4tmc = pyproj.Proj( P4_TM.format( **TM_PARAM ) )
df['PSF'] = p4tmc.get_factors( 
                df.P4lng, df.P4lat, radians=False ).meridional_scale
df['PSF_ppm'] = df['PSF'].apply( PPM )
for col in ('MSL','TM_E','TM_N'):
    df[col] = df[col].map( '{:,.3f}'.format )
print( df[['WGS_Lat', 'WGS_Lng', 'MSL', 'TM_E','TM_N', 'PSF_ppm' ]])

df = gpd.GeoDataFrame( df, crs='epsg:4326',
        geometry=gpd.points_from_xy( df.P4lng, df.P4lat ) )
df['MSL'] = df['MSL'].astype( float )
df.to_file('df_DOH107.gpkg', layer='GPS0107', driver='GPKG' )

df.to_file('df_DOH107.kml', driver='KML' )

df_csv = df[['Stn','UTM_N','UTM_E','MSL']].copy()
df_csv['ZONE'] = 47
df_csv = df_csv[['Stn','ZONE', 'UTM_N','UTM_E','MSL']]
for col in ('UTM_N','UTM_E','MSL'):
    df_csv[col] = df_csv[col].map( '{:.3f}'.format )
df_csv.to_csv( 'DOH107_GNSS_UTM47.csv', index_label='Seq' )
#import pdb; pdb.set_trace()



