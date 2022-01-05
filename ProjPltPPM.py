#
#
# ProjectPlotPPM : plot PPM or MSL color-shaded map for map projection
#                  design
#  13 Apr 2021 : Phisan Santitamnont , Chulalongkorn University
#                Phisan.Chula@gmail.com
#
import os,sys
from pyproj import Proj, CRS, Transformer, Geod
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Polygon,LineString
from math import cos,sin,radians
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

################################################################
class ProjectPlotPPM:
    ''' ProjectPlotPPM will make dual-plot (ax1 and ax2)
    '''
    def __init__(self, SYMB, crs, DF_LDP, RESULT_DIR):
        #import pdb; pdb.set_trace()
        self.crs = crs
        self.DF_LDP = DF_LDP
        if SYMB=='TOPO':
            COL1 = 'MSL'         ; COL2 = 'HSF'
            TITLE1 = 'Mean Sea Level : MSL (meter) '
            TITLE2 = 'Height Scale Factor from HAE : HSF in ppm'
        else:  # UTM etc...
            COL1 = f'{SYMB}_PSF' ; COL2 = f'{SYMB}_CSF'
            TITLE1 = 'Projection Point Scale Factor : PSF in ppm'
            TITLE2 = 'Combined Scale Factor : CSF in ppm'

        fig,(ax1,ax2) = plt.subplots( 1,2, figsize=(24,12 ) )
        sc1 =ax1.scatter( DF_LDP.lng,DF_LDP.lat,c=DF_LDP[COL1], cmap='terrain' )
        col_sch = 'binary' if SYMB=='UTM' else 'Spectral'
        sc2 = ax2.scatter( DF_LDP.lng,DF_LDP.lat,c=DF_LDP[COL2],cmap=col_sch )
        for ax,sc,TITLE in [ (ax1,sc1,TITLE1) ,(ax2,sc2,TITLE2) ] :
            self._PlotMapAxis( ax )
            ax.set_title( TITLE)
            clb = fig.colorbar( sc, ax=ax )
            clb_title = 'MSL(meter)' if 'MSL' in TITLE else 'ppm'
            clb.ax.set_title( clb_title )
        if crs is None:
            title = 'Topography SRTM15+ (GSD=450m) and Geoid Model EGM-2008'
        else:
            self._PlotProjectionAxis( crs, ax1, ax2 )
            title = SYMB + ' : ' + crs.to_proj4()
        plt.suptitle( title, size=16, y=0.990 ) 
        plt.tight_layout()
        #import pdb; pdb.set_trace()
        if SYMB=='TOPO':
            self.PLOT_FIG = f'{RESULT_DIR}/HSF_{SYMB}.png' 
        else:
            self.PLOT_FIG = f'{RESULT_DIR}/CSF_{SYMB}.png' 
        plt.savefig( self.PLOT_FIG )
        #plt.show()

    ##############################################################
    def _PlotMapAxis( self, ax):
        df = self.DF_LDP
        ax.plot(df.lng.mean(),df.lat.mean(), marker="+" , 
                 color='k', mew=4, ms=24, alpha=0.5 )
        #################################################
        axis_size = max(df.lng.max()-df.lng.min(),df.lat.max()-df.lat.min()) 
        half = axis_size/2.
        ax.set_xlim( df.lng.mean()-half, df.lng.mean()+half ) 
        ax.set_ylim( df.lat.mean()-half, df.lat.mean()+half ) 
        #################################################
        ax.ticklabel_format( useOffset=False)
        ax.tick_params(axis='x',rotation=45)
        ax.set( xlabel='Longitude', ylabel='Latitude' )
        ax.set_aspect('equal', 'box')
        ax.grid( True )

    ##############################################################
    def _PlotProjectionAxis(self, crs, ax1, ax2  ):
        #import pdb; pdb.set_trace()
        crs_prop = crs.to_dict()
        kws = dict( c='r', ls='dashdot' , lw=6, alpha=0.3 )
        if crs_prop['proj'] == 'utm':
            pass
        elif crs_prop['proj'] == 'tmerc':
            for ax in (ax1,ax2):
                ax.axvline( crs_prop['lon_0'], **kws )
        elif (crs_prop['proj']=='lcc') and not ('lat_2' in crs_prop.keys()):
            for ax in (ax1,ax2):
                ax.axhline( crs_prop['lat_0'], **kws )
        elif (crs_prop['proj']=='lcc') and ( 'lat_2' in crs_prop.keys() ):
            for ax in (ax1,ax2):
                ax.axhline( crs_prop['lat_2'], **kws )   
                ax.axhline( crs_prop['lat_0'], **kws )   
                ax.axhline( crs_prop['lat_1'], **kws )   
        elif crs_prop['proj']=='omerc':
            lat_0 = crs_prop['lat_0'] 
            lonc  = crs_prop['lonc'] 
            alpha = crs_prop['alpha'] 
            xmin, ymin, xmax, ymax =  self.DF_LDP.total_bounds
            poly = Polygon([ [xmin,ymin],[xmin,ymax],[xmax,ymax],[xmax,ymin] ])
            r = 10*max( xmax-xmin, ymax-ymin )
            rcos,rsin = r*cos( radians(alpha) ),r*sin( radians(alpha) )
            line = LineString([(lonc-rsin,lat_0-rcos),(lonc+rsin,lat_0+rcos)])
            pnts = line.intersection( poly )
            ((x1,y1),(x2,y2))  = list(pnts.coords)
            for ax in (ax1,ax2):
                ax.plot( [x1,x2],[y1,y2], **kws )
        else:
            #import pdb;pdb.set_trace()
            print(f'***ERROR*** unsupported {PROJ.definition}' )
            raise

##############################################################
if __name__=="__main__":
    ppm = ProjectPlotPPM()





