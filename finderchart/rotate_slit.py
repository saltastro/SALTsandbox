import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import astropy.coordinates as astcoo

from finder_chart import *

###########################################################################q######################################################################################################################################      
def get_coords( obj_name ):
    ''' Attempts a SIMBAD Sesame query with the given object name. '''
    try: 
        print( '\nQuerying SIMBAD database for {}...'.format(repr(obj_name)) )
        coo = astcoo.name_resolve.get_icrs_coordinates( obj_name )
        ra = coo.ra.to_string( unit='h', precision=2, sep=' ', pad=1 )
        dec = coo.dec.to_string( precision=2, sep=' ', alwayssign=1, pad=1 )
        
        print( 'The following ICRS J2000.0 coordinates were retrieved:\nRA = {}, DEC = {}\n'.format(ra, dec) )
        return coo, ra, dec
        
    except Exception as err:     #astcoo.name_resolve.NameResolveError
        print( 'ERROR in retrieving coordinates...' )
        print( err )
        return None, None, None

 
#######################################################################################################################
class RotateSlit( object ):

    #====================================================================================================
    def __init__(self, plot, ra, dec, pa=0):
        self.selection = None
        self.fig = plot._figure
        self.ax = ax = fig.axes[0]

        #slit width and height
        w, h = 10./60, 20./3600          

        #determine rotation between pixel and world (sky) coordinates
        px = np.array( plot.world2pixel([ra,ra+w],[dec, dec]) )
        dpx = px[:,0] - px[:,1]
        self.pa0 = pa0 = np.degrees( np.arctan2(dpx[1], dpx[0]) )       #origin for position angle

        self.background = fig.canvas.copy_from_bbox(fig.bbox)

        plot.show_rectangles(ra, dec, w, h, edgecolor='r', alpha=0.7, lw = 2, picker=10)

        self.fig2data = self.ax.transData.inverted()  #transformation from display to pixel coordinates
        self.pxc = xc, yc = plot.world2pixel( ra, dec )        #slit center in pixel coords   

        self.slit = slit = plot.get_layer('rectangle_set_1')
        trans = Affine2D().rotate_deg_around(xc, yc, pa+pa0) + ax.transData
        slit.set_transform( trans )

        plot.add_label(1.05, -0.05, "PA=%3.1f" %pa, relative=True, style='italic', weight='bold')
        self.pa_label = plot.get_layer('label_7')

        #fig.canvas.draw()

    #====================================================================================================
    def on_pick(self, event):
        self.selection = event.artist


    #====================================================================================================
    def on_motion(self, event):

        ax = self.ax
        if event.inaxes != ax:
            return

        if self.selection:
            xc, yc = pxc = self.pxc
            slit = self.slit
            xy = x,y = self.fig2data.transform((event.x, event.y)) - pxc
            th = np.degrees( np.arctan2(y,x) )
            pa = th + self.pa0

            fig.canvas.restore_region(self.background)
            trans = Affine2D().rotate_deg_around(xc, yc, pa) + ax.transData
            slit.set_transform( trans )
            self.pa_label.set_text( "PA=%3.1f" %th )
            #print( self.pa_label.get_text() )

            self.slit.set_alpha(0.3)
            ax.draw_artist(self.slit)
            ax.draw_artist(self.pa_label)
            fig.canvas.blit(self.fig.bbox)

            #self.fig.canvas.draw()

    #====================================================================================================
    def on_release(self, event):

        if self.selection:
            fig.canvas.restore_region(self.background)
            self.slit.set_alpha(0.7)
            ax.draw_artist(self.slit)
            ax.draw_artist(self.pa_label)
            fig.canvas.blit(self.fig.bbox)

        self.selection = None

    #====================================================================================================
    def connect(self):
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        #self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
       

#parse ommand line arguments
import argparse, os
parser = argparse.ArgumentParser(description='')
parser.add_argument('-pi', '--principle-investigator', nargs='+', dest='pi', default=[''],
                        help='The principle investigator for the proposal')
parser.add_argument('-c', '--proposal-code', dest='propcode',  
                        help='The proposal code')
parser.add_argument('-t', '--target', nargs='+', dest='target', 
                        help='Any named target')
parser.add_argument('-o', '--output-dir', dest='outdir', default='',
                        help='Optional output directory for saving the figure.')
args = parser.parse_args()

print( args )

obj_name = ' '.join( args.target )
pi = ' '.join( args.pi )
propcode = args.propcode



if os.path.isdir( args.outdir ):
    from matplotlib import rc
    rc( 'savefig', directory=args.outdir )

#Setup access through UCT proxy 
#proxy = 'www.ast.uct.ac.za:33128'
#proxy = urllib2.ProxyHandler({'http': proxy})
#opener = urllib2.build_opener(proxy)
#urllib2.install_opener(opener)

#get coordinates from object name
coo, _, _ = get_coords( obj_name )
ra, dec = coo.ra.deg, coo.dec.deg

#get image
allowed_servers = ('poss2ukstu_blue', 'poss2ukstu_red', 'poss2ukstu_ir',
                    'poss1_blue', 'poss1_red',
                    'all'                    )
for imserver in allowed_servers:
    try:
        print( 'Retrieving FITS image from DSS server: %s' %imserver )
        hdu = get_dss(imserver, ra, dec)
        break
    except Exception as err:
        print( 'DSS image retrieval failed with:\n%s\n' %err )


title = '{} ({}; {})'.format( obj_name, propcode, pi )
plot = init_plot(hdu, imserver, title, ra, dec)

fig = plot._figure
ax = fig.axes[0]

r = RotateSlit(plot, ra, dec, 0)
r.connect()

plt.show()
