#
#
# uiProjDesign : Map Project Design by reading the design parameter.
#              Map projection paramaters are defacto Proj4 
#              starndard format.  Results will be produced in 
#              summary of Point Scale Factor (PSF) Height Scale
#              Facotr (HSF) and Combined Scale Factor (CSF). 
#              Visual color-shaded PNG will be produced namely
#     - TOPO.png : MSL of the analysed ROI and HSF
#     - UTM.png  : PSF and CSF for UTM zone47/48 projection
#     - TMC.png  : PSF and CSF for Transverse Mercator 
#     - LCC.png  : PSF and CSF for Labert Conical Conformal projection
#     - OMC.png  : PSF and CSF for Oblique Transverse Mercator projection
#
# P.Santitamnont (2021) phisan.chula@gmail.com
#
#
import os,sys,gc
import pyproj
from pyproj import Proj, CRS, Transformer, Geod
import pandas as pd
import geopandas as gpd
import numpy as np
import yaml
from pathlib import Path
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from math import cos,sin,radians
import pygeodesy.geoids as geoids
from collections import OrderedDict, namedtuple
from tabulate import tabulate
import argparse
from ProjDesign import *
from ProjPltPPM import *

pd.set_option("display.max_columns", None)
pd.options.display.float_format ='{:,.0f}'.format

#####################################################################
def ReviewTopography( ldp ):
    RP.PRN('*** [1] Topography Reviewing ... ***')  
    RP.PRN( ldp.df_sampl[['MSL','UNDUL','HAE','HSF']].describe() )
    RP.PRN('============================================================')
    RP.PRN(f'Mean topo MSL over ROI : {ldp.MeanMSL:.0f}  m.' )
    RP.PRN(f'Mean topo HAE over ROI : {ldp.MeanHAE:.0f} m.'  )
    RP.PRN(f'k_0 PP coincided with topo : {ldp.k_0} (exact 7-digit)' )
    RP.PRN('Mean Latitude of ROI    : ', ldp.MeanLat )
    #RP.PRN('Suggested first standard parallel      : ', ldp.SP1 )
    #RP.PRN('Suggested secondary standard parallel  : ', ldp.SP2 )
    RP.PRN('Mean Longitude of ROI   : ', ldp.MeanLng )
    RP.PRN('============================================================')
    if 'LDP' not in ldp.YAML.keys(): # if no LDP sect, create one for UTM
        ldp.YAML['LDP'] = dict()
    pltppm = ProjectPlotPPM( 'TOPO', None, ldp.df_sampl, RESULT_DIR) # 
    #import pdb; pdb.set_trace()
    ldp_ = OrderedDict( (['UTM', CRS(ldp.UTM_EPSG)], ) )
    ldp_.update( OrderedDict(ldp.YAML['LDP'] ) )
    RP.PRN('*** [2] SUGGESTION DESIGN ... ***')  
    LDP_SUGGEST=\
f'''LDP :
   TMC : +proj=tmerc +lat_0={ldp.MeanLat[1]} +lon_0={ldp.MeanLng[1]} 
         +k_0={ldp.k_0} +x_0=+0 +y_0=+0 +ellps=GRS80 +units=m +no_defs
   LCC : proj=lcc +lon_0={ldp.MeanLng[1]} +lat_0={ldp.MeanLat[1]}
         +lat_1={ldp.MeanLat[1]}  +k_0={ldp.k_0} 
         +x_0=+0 +y_0=+0 +ellps=GRS80 +units=m +no_defs
   OMC : proj=omerc +lat_0={ldp.MeanLat[1]} +lonc={ldp.MeanLng[1]} 
         +x_0=+0 +y_0=+0 +alpha=+45 +k_0={ldp.k_0} 
         +ellps=GRS80 +units=m +no_defs '''
    RP.PRN( LDP_SUGGEST )
    return ldp_

####################################################################
if __name__=="__main__":
    #VER = pyproj.show_versions()
    RP.PRN('********************** uiDesign ************************')
    RP.PRN('*** P.Santitamnont (phisan.chula@gmail.com) Oct.2021 ***')
    RP.PRN('********************************************************')
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="YAML configuration FILE",
                        type=pathlib.Path )
    parser.add_argument("-p", "--process", help="processing the map design",
                        action="store_true")
    ARGS = parser.parse_args()
    YAMLFILE = ARGS.file
    ###############################################################
    RP.PRN( f'Reading YAML {YAMLFILE} ...')
    RESULT_DIR = YAMLFILE.parents[0].joinpath('RESULTS')
    if Path( RESULT_DIR ).mkdir( parents=True, exist_ok=True ):
        RP.PRN(f'Creating {RESULT_DIR} ...')
    RP.PRN(f'Writing all results in directory/folder {RESULT_DIR} ...')
    ldp  = MapDesign( YAMLFILE, RESULT_DIR )
    ldp_ = ReviewTopography( ldp )
    if not ARGS.process:
        sys.exit(-1)
    ###############################################################
    RP.PRN('*** [3] PROCESSING LDP and PSF/CSF Analysis ... ***')  
    for SYM, p4txt in ldp_.items():
        if isinstance( p4txt, CRS):
            crs = p4txt
        else:
            #import pdb; pdb.set_trace()
            crs = CRS.from_proj4( p4txt)
        title = SYM+' : ' + crs.to_proj4()
        RP.PRN( 70*':'); RP.PRN( title ); RP.PRN( 70*':' )
        ldp.CalcProjScaleFactor( crs, SYMBOL=SYM )
        COLS = [ f'{SYM}_E',f'{SYM}_N', f'{SYM}_PSF',f'{SYM}_CSF' ]
        for col in COLS:
            ldp.df_sampl.style.format( { col: '{:.0f}' })
        RP.PRN( ldp.df_sampl[ COLS ].describe() )
        if 0:
            fig,ax = plt.subplots()
            ldp.df_sampl[f'{SYM}_CSF'].hist(ax=ax)
            plt.title( f'{SYM}_CSF' )
            plt.savefig( RESULT_DIR.joinpath( f'hst_CSF_{SYM}.png') )
        pltppm = ProjectPlotPPM( SYM, crs, ldp.df_sampl, RESULT_DIR ) 
        RP.PRN( f'Plotting result in {pltppm.PLOT_FIG} ...' )
        gc.collect()

    ###############################################################
    RP.PRN('*** [4] Write report and sampling GIS file ... ***')  
    RP.to_file( RESULT_DIR.joinpath( 'Report.txt' ) )
    sampl = ldp.df_sampl.copy()
    for col in sampl.columns:
        if col[-4:]=='_CSF' or col[-4:]=='_PSF' or col=='HSF' :
            sampl[col] = sampl[col].round(1)
    sampl.to_file( RESULT_DIR.joinpath('ldp_sampling.gpkg') , 
                            driver="GPKG", layer='sampl')
