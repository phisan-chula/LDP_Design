#
#
# ProjDesign : ProjDesign.py reads the design parameter in YAML format.
#              Map projection paramaters in section "LDP" has from
#              defacto Proj4 standard format.  Results will be produced in 
#              summary of Point Scale Factor (PSF) Height Scale
#              Facotr (HSF) and Combined Scale Factor (CSF). 
#              Visual color-shaded PNG will be produced namely
#     - TOPO.png : MSL of the analysed ROI and HSF
#     - UTM.png  : PSF and CSF for UTM zone47/48 projection
#     - TMC.png  : PSF and CSF for Transverse Mercator 
#     - LCC.png  : PSF and CSF for Lambert Conical Conformal projection
#     - OMC.png  : PSF and CSF for Oblique Transverse Mercator projection
#
# P.Santitamnont (2021) phisan.chula@gmail.com
#
#
import os,sys,gc,yaml
from pyproj import Proj, CRS, Transformer, Geod
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from math import cos,sin,radians,pow
import pygeodesy.geoids as geoids
import pygeodesy.geohash as gh
import pygeodesy.dms as dms
import pygeodesy as pgd
import rasterio as rio
from shapely.geometry import Polygon,LineString
from collections import OrderedDict, namedtuple
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
import warnings

WGS= CRS( 4326 )
UTMz47 = Proj( 'EPSG:32647' )
UTMz48 = Proj( 'EPSG:32647' )
ECEF = Proj(proj='geocent', ellps='WGS84', datum='WGS84')

######################################################################
class REPORT:
    def __init__( self, QUIET=False ):
        self.LINES = list()
        self.QUIET = QUIET
        pass
    def PRN( self, *args ):
        txt = ''
        for arg in args:
            txt = txt + str( arg )
        self.LINES.append( txt )
        if not self.QUIET:
            print( txt )
    def to_file(self, OUT_FILE ):
        LINES = "\n".join( self.LINES )
        with open(OUT_FILE,'w') as f:
            f.writelines( LINES ) 
RP = REPORT( )

##############################################################
def PPM( k ):
    return (k-1.)*1000.*1000

##############################################################
def RoundToMin( dd ): 
    DD,MM = divmod( int(dd*60+0.5)  ,60)  # round to 'minute'
    ddmm = DD+MM/60.
    ddmm_str = dms.toDMS(ddmm)[:-5] 
    return ddmm_str, ddmm 
##############################################################
def UTMZone( lng,lat ):
    zone = int((int(lng) + 180) / 6) + 1
    if lat > 0:
        return f'EPSG:326{zone}'
    else:
        return f'EPSG:327{zone}'

####################################################################
class MapDesign:
    def __init__(self, MAPFILE, RESULT_DIR):
        self.RESULT_DIR = RESULT_DIR
        with open( MAPFILE ) as f:
            self.YAML = yaml.load( f, Loader=yaml.FullLoader )
        self.FILE_ROI = RESULT_DIR.parents[0].joinpath( self.YAML['KML'] )
        if 'BUFF' in self.YAML.keys():
            BUFF_M=self.YAML['BUFF']
        if 'GRID' in self.YAML.keys():
            GRID_N=self.YAML['GRID']
        self.df_sampl = self.GenSamplPnt( GRID_N=GRID_N,BUFF_M=BUFF_M )
        self.df_sampl = self.GetDEM( self.df_sampl )

        self.MapInitial()

        self.df_sampl['Rg'] = self.df_sampl['lat'].apply(
                  pgd.Ellipsoids.WGS84.rocGauss )
        self.df_sampl['HSF'] = PPM( 
                self.df_sampl['Rg']/(self.df_sampl['Rg']+self.df_sampl['HAE']) )

    ##############################################################
    def MapInitial(self):
        MeanLat = self.df_sampl.geometry.y.mean()
        self.MeanLat = RoundToMin( MeanLat )
        dlat_min = self.df_sampl.geometry.y.min()
        dlat_13 = (self.df_sampl.geometry.y.max()-dlat_min)/3
        self.SP1 = RoundToMin( dlat_min + dlat_13 )
        self.SP2 = RoundToMin( dlat_min + 2*dlat_13 )
        self.MeanLng = RoundToMin( self.df_sampl.geometry.x.mean() )
        self.MeanMSL = self.df_sampl.MSL.mean()
        self.MeanHAE = self.df_sampl.HAE.mean()
        self.MeanRg = pgd.Ellipsoids.WGS84.rocGauss(MeanLat)
        k_0 = 1.+ self.MeanHAE/self.MeanRg
        self.k_0 = round( k_0, 7 )
        #import pdb;pdb.set_trace()
        return 

    ##############################################################
    def GenSamplPnt(self, GRID_N=3, BUFF_M=300):
        ''' TENTH   1=10 ; 2=100 ; 3=1000
        '''
        RP.PRN( f'Reading KML {self.FILE_ROI} ...' )
        df = gpd.read_file( self.FILE_ROI , driver='KML' )
        if 'GRID' in self.YAML.keys():
            GRID_N = self.YAML['GRID']
        GEOM1st = df.iloc[0].geometry
        if GEOM1st.geom_type == 'LineString':
            #import pdb;pdb.set_trace()
            RP.PRN( 'Input Linestring with length {:,.0f} m.'.\
                    format( GEOM1st.length*111_000) )
            warnings.simplefilter('ignore')
            df_buf = df.buffer( BUFF_M/111_000, cap_style=2, join_style=2 )
            df_buf = gpd.GeoDataFrame( crs=4326, 
                                geometry=gpd.GeoSeries(df_buf) )
        elif GEOM1st.geom_type == 'Polygon':
            df_buf = df.copy()
        else:
            RP.PRN('*** geometry type not defined ***')
            raise '***ERROR***'

        self.UTM_EPSG = UTMZone( GEOM1st.centroid.x, GEOM1st.centroid.x )
        df_samp = self.RectangleSampling( df_buf, GRID_N )
        df_samp = gpd.sjoin( df_samp, df_buf, how='inner', op='intersects')
        df_samp.reset_index( inplace=True )
        def Geoh(row):
            return gh.encode( row.geometry.y, row.geometry.x, precision=7 )
        df_samp['Geohash'] = df_samp.apply( Geoh , axis='columns' )
        df_samp = df_samp[['Geohash','geometry']]
        return df_samp

    ##############################################################
    def RectangleSampling(self, df_poly,GRID_N ):
        #import pdb;pdb.set_trace()
        df_poly = df_poly.to_crs( crs=self.UTM_EPSG )
        bnd = df_poly.iloc[0].geometry.bounds
        pwn = pow(10,GRID_N)
        E_ = np.arange( round( bnd[0],-GRID_N), 
                        round( bnd[2],-GRID_N )+pwn, pwn )
        N_ = np.arange( round( bnd[1],-GRID_N), 
                        round( bnd[3],-GRID_N )+pwn, pwn )
        E_,N_ = np.meshgrid( E_,N_, indexing='xy' )
        df_samp = pd.DataFrame( {'UTM_E': E_.ravel(), 'UTM_N': N_.ravel()} )
        df_samp = gpd.GeoDataFrame( df_samp, crs=self.UTM_EPSG,
                geometry=gpd.points_from_xy( df_samp.UTM_E, df_samp.UTM_N ) )
        return df_samp.to_crs( 4326 )
        
    ##############################################################
    def GetDEM(self, df_samp ):
        coords = list( zip( df_samp.geometry.x, df_samp.geometry.y ) )
        with rio.open( 'DEM/SRTM15p_v21_Thai.tif' ,'r' ) as src:
            dem = src.read(1)
            try:
                df_samp['MSL'] = [msl[0] for msl in src.sample( coords )]
            except:
                RP.PRN('*** ERROR Profile Data ***')
        #import pdb;pdb.set_trace()
        if 'CLAMP_MSL' in self.YAML.keys():
            ClmMSL   = self.YAML['CLAMP_MSL']
            RP.PRN( '***WARNING*** CLAMP_MSL : {} m. (DEM biased)'.\
                    format(ClmMSL) ) 
            AvgMSL  = self.df_sampl.MSL.mean()
            df_samp['MSL'] = df_samp['MSL']-AvgMSL+ClmMSL
        df_samp['lng'] = df_samp.geometry.x
        df_samp['lat'] = df_samp.geometry.y
        # use internal pyGeodesy/GeoidGPM()
        Model = geoids.GeoidPGM('DEM/egm2008-5.pgm')
        df_samp['UNDUL'] = Model.height( 
                 df_samp.geometry.y , df_samp.geometry.x )
        df_samp['HAE'] = df_samp['UNDUL'] + df_samp['MSL']
        df_samp['HAE'] = df_samp['HAE'].astype( 'int')
        return df_samp

    ##############################################################
    def CalcProjScaleFactor( self, crs: CRS, SYMBOL=None ):
        df = self.df_sampl
        if SYMBOL is None:
            PRJ = crs.to_dict()['proj']    
        else:
            PRJ = SYMBOL
        TRANSF = Transformer.from_crs( 4326 , crs, always_xy=True ) 
        lnglat = list( zip(df.lng,df.lat) )
        df[ [ f'{PRJ}_E', f'{PRJ}_N'] ] = np.array( list( 
            TRANSF.itransform( lnglat, switch=False, radians=False ) ) )
        #import pdb;pdb.set_trace()
        facts = Proj(crs).get_factors( 
                list(df.lng), list(df.lat), radians=False )
        df[f'{PRJ}_sf'] = facts.meridional_scale
        df[f'{PRJ}_PSF'] = df[f'{PRJ}_sf'].apply( PPM )
        df[f'{PRJ}_CSF'] = df[f'{PRJ}_PSF'] + df[f'HSF'] 
        return
