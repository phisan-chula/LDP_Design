#
#
# uiLDP_PostProc : Map Project Design by reading the design parameter.
#              Map projection paramaters are defacto Proj4 
#              starndard format.  Results will be produced in 
#              summary of Point Scale Factor (PSF) Height Scale
#              Facotr (HSF) and Combined Scale Factor (CSF). 
#              A CSV file will be read converted to a LDP specified 
# P.Santitamnont (2021) phisan.chula@gmail.com
#
#
import os,sys,gc
from pyproj import Proj, CRS, Transformer, Geod
import pandas as pd
import geopandas as gpd
import numpy as np
import yaml
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
from matplotlib import ticker
from math import cos,sin,radians
import pygeodesy.geoids as geoids
import pygeodesy as pgd
import pyproj
import argparse
from collections import OrderedDict, namedtuple
from ProjDesign import *

####################################################################
def MakeHSF(df_samp ):
    # use internal pyGeodesy/GeoidGPM()
    Model = geoids.GeoidPGM('DEM/egm2008-5.pgm')
    df_samp['UNDUL'] = Model.height(
             df_samp.geometry.y , df_samp.geometry.x )
    df_samp['HAE'] = df_samp['UNDUL'] + df_samp['MSL']
    df_samp['HAE'] = df_samp['HAE'].astype( 'int')
    def SF(row):
        RG  = pgd.Ellipsoids.WGS84.rocGauss(row.lat)
        HSF = PPM( RG/(RG+row.HAE) )
        return pd.Series( [RG,HSF] )
    df_samp[['RG','HSF']] = df_samp.apply( SF, axis=1 )
    #import pdb;pdb.set_trace()
    return df_samp

####################################################################
def PlotPoint( df,SYM,LDP,RESULT_DIR ):
    E = f'{SYM}_E' ;  N = f'{SYM}_N'
    fig,ax = plt.subplots( figsize=(10,10))
    df.plot.scatter( x=E ,y=N,marker='^',s=128,alpha=0.5,ax=ax )
    ax.set_xlabel( f'{E} (m)' )
    ax.set_ylabel( f'{N} (m)' )
    def SqLim( ):
        dx= df[E].max()-df[E].min(); dy= df[N].max()-df[N].min()
        sq = 1.4*max(dx,dy)/2  # over 20%
        x_m = df[E].mean(); y_m = df[N].mean()
        return ( x_m-sq, x_m+sq, y_m-sq, y_m+sq)
    limx1,limx2,limy1,limy2 = SqLim()
    #import pdb;pdb.set_trace()
    ax.set_xlim( limx1,limx2 )
    ax.set_ylim( limy1,limy2 )
    tick = StrMethodFormatter('{x:,.0f}')
    ax.xaxis.set_major_formatter( tick )
    ax.yaxis.set_major_formatter( tick )
    for i,row in df.iterrows():
        ann = ax.annotate( row.Stn,
                  xy=(row[E], row[N]), xycoords='data',
                  xytext=(32, -32),  textcoords='offset points',
                  size=20, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                       connectionstyle="arc3,rad=-0.2", fc="w"),)
    plt.grid(True)
    plt.setp(ax.get_xticklabels(), **{"rotation" :90})
    def_str = LDP.definition ; pos = len(def_str)//2
    plt.title( def_str[:pos]+'\n'+def_str[pos:] )
    fig.tight_layout()
    plt.savefig( RESULT_DIR.joinpath( f'PNT_{SYM}.png' ) ) 
    #plt.show()

####################################################################
pd.set_option("display.max_rows", None)
####################################################################
if __name__=="__main__":
    RP.PRN('********************* uiPostProc ************************')
    RP.PRN('*** P.Santitamnont (phisan.chula@gmail.com) Oct.2021 ***')
    RP.PRN('*********************************************************')
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="YAML configuration FILE",
                        type=pathlib.Path )
    ARGS = parser.parse_args()
    #import pdb; pdb.set_trace()
    YAMLFILE = ARGS.file

    ###############################################################
    RP.PRN( f'Reading YAML {YAMLFILE} ...')
    RESULT_DIR = YAMLFILE.parents[0].joinpath('RESULTS')
    if pathlib.Path( RESULT_DIR ).mkdir( parents=True, exist_ok=True ):
        RP.PRN(f'Creating {RESULT_DIR} ...')
    RP.PRN(f'Writing all results in directory/folder {RESULT_DIR} ...')

    ################################################################
    with open( YAMLFILE ) as f:
        YAML = yaml.load( f, Loader=yaml.FullLoader )
    if 'CSV_UTM' not in YAML.keys():
        print('*** Warning no CSV_UTM section in config file...')
        print('***    No transformation from UTM to any LDPs...')
        sys.exit( -1 )
    CSV = RESULT_DIR.parents[0].joinpath( YAML['CSV_UTM'] )
    df = pd.read_csv( CSV )
    assert( len(df.ZONE.value_counts()==1) ) # all 46/47
    assert( df.ZONE.iloc[0] in (47,48) )
    UTM_EPSG = 'EPSG:326'+ str(df.ZONE.iloc[0])
    df = gpd.GeoDataFrame( df , crs=UTM_EPSG,
            geometry=gpd.points_from_xy( x=df.UTM_E, y=df.UTM_N ) )
    df = df.to_crs('EPSG:4326')
    df['lng'] = df.geometry.x; df['lat'] = df.geometry.y
    df = MakeHSF( df )

    #################################################################
    UTM = pyproj.CRS(UTM_EPSG)  # UTM zone47/47
    for SYM in YAML['LDP'].keys():
        ldp = pyproj.CRS.from_proj4( YAML['LDP'][SYM] )
        print( 115*'='); print( YAML['LDP'][SYM] )
        TR_UTM_ldp = pyproj.Transformer.from_crs( UTM, ldp, always_xy=True )
        ldp_EN = TR_UTM_ldp.transform( df.UTM_E, df.UTM_N )
        LDP_E = f'{SYM}_E';  LDP_N = f'{SYM}_N'
        df[LDP_E] = ldp_EN[0] ; df[LDP_N] = ldp_EN[1]
        #############################################################
        ldp = pyproj.Proj(  YAML['LDP'][SYM] )
        df['sf'] = ldp.get_factors( df.lng, df.lat, 
                            radians=False ).meridional_scale
        df['PSF'] = df['sf'].apply( PPM )
        df['CSF'] = df['PSF'] + df['HSF']
        #############################################################
        PlotPoint( df, SYM, ldp,  RESULT_DIR )

        #############################################################
        df_tab = df[['Stn','ZONE','UTM_E','UTM_N','MSL',LDP_E,LDP_N,
                     'HSF','PSF','CSF']].copy()
        for col in df_tab.columns:
            if col in [ LDP_E ,LDP_N , 'MSL', 'UTM_E', 'UTM_N']:
                df_tab[col] = df_tab[col].map('{:,.3f}'.format)
            elif col[-2:]=='SF' and len(col)==3:
                df_tab[col] = df_tab[col].map('{:.0f}'.format)
        df_tab.index += 1
        from tabulate import tabulate
        TABLE = tabulate( df_tab, df_tab.columns, tablefmt='pretty')
        print( TABLE )
        gc.collect()
        #import pdb; pdb.set_trace()

#####################################################################

