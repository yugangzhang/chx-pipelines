from __future__ import absolute_import, division, print_function

from databroker import DataBroker as db, get_images, get_table, get_events
from filestore.api import register_handler, deregister_handler
from filestore.retrieve import _h_registry, _HANDLER_CACHE

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import skxray.core.roi as roi
from datetime import datetime



import h5py
from filestore.retrieve import HandlerBase
from eiger_io.pims_reader import EigerImages

EIGER_MD_DICT = {
    'y_pixel_size': 'entry/instrument/detector/y_pixel_size',
    'x_pixel_size': 'entry/instrument/detector/x_pixel_size',
    'detector_distance': 'entry/instrument/detector/detector_distance',
    'incident_wavelength': 'entry/instrument/beam/incident_wavelength',
    'frame_time': 'entry/instrument/detector/frame_time',
    'beam_center_x': 'entry/instrument/detector/beam_center_x',
    'beam_center_y': 'entry/instrument/detector/beam_center_y',
    'count_time': 'entry/instrument/detector/count_time',
    'pixel_mask': 'entry/instrument/detector/detectorSpecific/pixel_mask',
}

class FixedEigerImages(EigerImages):
    def __init__(self, path, metadata):
        super().__init__(path)
        self._metadata = metadata
    
    @property
    def md(self):
        return self._metadata
    
    @property
    def dtype(self):
        return self.pixel_type
    
    @property
    def shape(self):
        return self.frame_shape

class LazyEigerHandler(HandlerBase):
    specs = {'AD_EIGER'} | HandlerBase.specs
    def __init__(self, fpath, frame_per_point, mapping=None):
        # create pims handler
        self.vals_dict = EIGER_MD_DICT.copy()
        if mapping is not None:
            self.vals_dict.update(mapping)
        self._base_path = fpath
        self.fpp = frame_per_point

    def __call__(self, seq_id):
        import h5py
        master_path = '{}_{}_master.h5'.format(self._base_path, seq_id)
        md = {}
        print('hdf5 path = %s' % master_path)
        with h5py.File(master_path, 'r') as f:
            md = {k: f[v].value for k, v in self.vals_dict.items()}
        # the pixel mask from the eiger contains:
        # 1  -- gap
        # 2  -- dead
        # 4  -- under-responsive
        # 8  -- over-responsive
        # 16 -- noisy
        pixel_mask = md['pixel_mask']
        pixel_mask[pixel_mask>0] = 1
        pixel_mask[pixel_mask==0] = 2
        pixel_mask[pixel_mask==1] = 0
        pixel_mask[pixel_mask==2] = 1
        md['framerate'] = 1./md['frame_time']
        # TODO Return a multi-dimensional PIMS seq
        return FixedEigerImages(master_path, md)

deregister_handler('AD_EIGER')
_HANDLER_CACHE.clear()
register_handler('AD_EIGER', LazyEigerHandler)


class Reverse_Coordinate(object):
    def __init__(self, indexable, mask):
        self.indexable = indexable
        self.mask = mask
        self.shape = indexable.shape
        self.length= len(indexable)
    def __getitem__(self, key ):      
        if self.mask is not None:
            img =self.indexable[key] * self.mask  
        else:
            img = self.indexable[key]
            
        if len(img.shape) ==3:
            img_=img[:,::-1,:]
        if len(img.shape)==2:
            img_=img[::-1,:] 
        return img_
 

class Masker(object):
    def __init__(self, indexable, mask):
        self.indexable = indexable
        self.mask = mask
    def __getitem__(self, key):        
        img =self.indexable[key] * self.mask        
        return img

  


def view_image(imgsr,i):
    #from ipywidgets import interact
    fig, ax = plt.subplots()
    ax.imshow(imgsr[i], interpolation='nearest', cmap='viridis',
                  origin='lower', norm= LogNorm(vmin=0.001, vmax=1e1 ) )
    ax.set_title("Browse the Image Stack")
    plt.show()
    
    
    
import time 
def view_image_movie(imgsr,sleeps=1, ims=0, ime = 1):    
    fig, ax = plt.subplots()  
    for i in range( ims, ime  ):
        ax.imshow(imgsr[i],  interpolation='nearest', cmap='viridis',
                  origin='lower', norm= LogNorm( vmin=0.001, vmax=1e1 ) )
        ax.set_title("images_%s"%i)
        time.sleep( sleeps )        
        plt.draw()
        if i!=ime-1:
            ax.cla()
        
        
        
def average_img( imgs, Ns=None,Ne = None ):
    ''' Do imgs average,
        Optiions:
        imgs: the image seriers
        Ns: the start image
        Ne: the last image
        e.g.,
        ave = average_img(imgs)'''
    import numpy as np 
    ave = np.zeros_like(imgs[0],dtype =float)
    #if Ns is None:Ns=0
    #if Ne is None:Ne=len(imgs)
    #if Ne>len(imgs):Ne=len(imgs)
    for i in range(Ns,Ne):
        ave += imgs[i]
    ave /= (Ne-Ns)
    return ave
        
    
    
    
import xray_vision
import xray_vision.mpl_plotting as mpl_plot  
from xray_vision.mpl_plotting import speckle
from xray_vision.mask.manual_mask import ManualMask

import skxray.core.roi as roi

import skxray.core.correlation as corr
import skxray.core.utils as utils


#GiSAXS
##########################################


def make_gisaxs_grid( qr_w= 10, qz_w = 12, dim_r =100,dim_z=120):
    y, x = np.indices( [dim_z,dim_r] )
    Nr = int(dim_r/qp_w)
    Nz = int(dim_z/qz_w)
    noqs = Nr*Nz
    
    ind = 1
    for i in range(0,Nr):
        for j in range(0,Nz):        
            y[ qr_w*i: qr_w*(i+1), qz_w*j:qz_w*(j+1)]=  ind
            ind += 1 
    return y 


def get_incident_angles( inc_x0, inc_y0, refl_x0, refl_y0, pixelsize=[75,75], Lsd=5.0):
    ''' giving: incident beam center: bcenx,bceny
                reflected beam on detector: rcenx, rceny
                sample to detector distance: Lsd, in meters
                pixelsize: 75 um for Eiger4M detector
        get incident_angle (alphai), the title angle (phi)
    '''
    px,py = pixelsize
    phi = np.arctan2( (refl_x0 - inc_x0)*px *10**(-6), (refl_y0 - inc_y0)*py *10**(-6) )    
    alphai = np.arctan2( (refl_y0 -inc_y0)*py *10**(-6),  Lsd ) /2.     
    #thetai = np.arctan2(  (rcenx - bcenx)*px *10**(-6), Lsd   ) /2.  #??   
    
    return alphai,phi 
    
       
def get_reflected_angles(inc_x0, inc_y0, refl_x0, refl_y0, thetai=0.0,
                         pixelsize=[75,75], Lsd=5.0,dimx = 2070.,dimy=2167.):
    
    ''' giving: incident beam center: bcenx,bceny
                reflected beam on detector: rcenx, rceny
                sample to detector distance: Lsd, in meters                
                pixelsize: 75 um for Eiger4M detector
                detector image size: dimx = 2070,dimy=2167 for Eiger4M detector
        get  reflected angle alphaf (outplane)
             reflected angle thetaf (inplane )
    '''    
    
    alphai, phi =  get_incident_angles( inc_x0, inc_y0, refl_x0, refl_y0, pixelsize, Lsd)
    print ('The incident_angle (alphai) is: %s'%(alphai* 180/np.pi))
    px,py = pixelsize
    y, x = np.indices( [dimy,dimx] )    
    #alphaf = np.arctan2( (y-inc_y0)*py*10**(-6), Lsd )/2 - alphai 
    alphaf = np.arctan2( (y-inc_y0)*py*10**(-6), Lsd )  - alphai 
    thetaf = np.arctan2( (x-inc_x0)*px*10**(-6), Lsd )/2 - thetai   
    
    return alphaf,thetaf, alphai, phi    
    
    

def convert_gisaxs_pixel_to_q( inc_x0, inc_y0, refl_x0, refl_y0, 
                               pixelsize=[75,75], Lsd=5.0,dimx = 2070.,dimy=2167.,
                              thetai=0.0, lamda=1.0 ):
    
    ''' giving: incident beam center: bcenx,bceny
                reflected beam on detector: rcenx, rceny
                sample to detector distance: Lsd, in meters                
                pixelsize: 75 um for Eiger4M detector
                detector image size: dimx = 2070,dimy=2167 for Eiger4M detector                
                wavelength: angstron               
                
        get: q_parallel (qp), q_direction_z (qz)
                
    '''         
    
    
    alphaf,thetaf,alphai, phi = get_reflected_angles( inc_x0, inc_y0, refl_x0, refl_y0, thetai, pixelsize, Lsd,dimx,dimy)
       
    pref = 2*np.pi/lamda
    
    qx = np.cos( alphaf)*np.cos( 2*thetaf) - np.cos( alphai )*np.cos( 2*thetai)  
    qy_ = np.cos( alphaf)*np.sin( 2*thetaf) - np.cos( alphai )*np.sin ( 2*thetai)    
    qz_ = np.sin(alphaf) + np.sin(alphai)   
    
    qy = qz_* np.sin( phi) + qy_*np.cos(phi) 
    qz = qz_* np.cos( phi) - qy_*np.sin(phi)   
    
    qr = np.sqrt( qx**2 + qy**2 ) 
    
    
    return qx*pref  , qy*pref  , qr*pref  , qz*pref  
        
    
    
    
def get_qedge( qstart,qend,qwidth,noqs,  ):
    ''' DOCUMENT make_qlist( )
    give qstart,qend,qwidth,noqs
    return a qedge by giving the noqs, qstart,qend,qwidth.
           a qcenter, which is center of each qedge 
    KEYWORD:  None    ''' 
    import numpy as np 
    qcenter = np.linspace(qstart,qend,noqs)
    #print ('the qcenter is:  %s'%qcenter )
    qedge=np.zeros(2*noqs) 
    qedge[::2]= (  qcenter- (qwidth/2)  ) #+1  #render  even value
    qedge[1::2]= ( qcenter+ qwidth/2) #render odd value
    return qedge, qcenter    
    
    
def get_qmap_label( qmap, qedge ):
    import numpy as np
    '''give a qmap and qedge to bin the qmap into a label array'''
    edges = np.atleast_2d(np.asarray(qedge)).ravel()
    label_array = np.digitize(qmap.ravel(), edges, right=False)
    label_array = np.int_(label_array)
    label_array = (np.where(label_array % 2 != 0, label_array, 0) + 1) // 2
    label_array = label_array.reshape( qmap.shape )
    return label_array
        
    
    
def get_qzrmap(label_array_qz, label_array_qr, qz_center, qr_center   ):
    qzmax = label_array_qz.max()
    label_array_qr_ = np.zeros( label_array_qr.shape  )
    ind = np.where(label_array_qr!=0)
    label_array_qr_[ind ] =  label_array_qr[ ind ] + 1E4  #add some large number to qr
    label_array_qzr = label_array_qz * label_array_qr_  
    
    #convert label_array_qzr to [1,2,3,...]
    uqzr = np.unique( label_array_qzr )[1:]
    
    uqz = np.unique( label_array_qz )[1:]
    uqr = np.unique( label_array_qr )[1:]
    #print (uqzr)
    label_array_qzr_ = np.zeros_like( label_array_qzr )
    newl = np.arange( 1, len(uqzr)+1)
    
    qzc =list(qz_center) * len( uqr )
    qrc= [  [qr_center[i]]*len( uqz ) for i in range(len( uqr ))  ]
    
    for i, label in enumerate(uqzr):
        #print (i, label)
        label_array_qzr_.ravel()[ np.where(  label_array_qzr.ravel() == label)[0] ] = newl[i]    
    
    
    return np.int_(label_array_qzr_), np.array( qzc ), np.concatenate(np.array(qrc ))
    
        
    
def get_qr_intensity_series( qr, data,vert_rect, show_roi=True ):
    V_K_label_array = roi.rectangles(vert_rect, data.shape)  #(y,x, hight, wdith)
    qr_ = qr  *V_K_label_array
    data_ = data*V_K_label_array 
    if False:
        fig, ax = plt.subplots()
        im = plt.imshow(data_,origin='lower',norm= LogNorm( vmin=.1, vmax=1e0 ) )
        fig.colorbar(im)
        plt.show()
    
    data_ave = np.average( data_, axis=0)
    qr_ave = np.average( qr_, axis=0)
    
    if show_roi:
        fig, ax = plt.subplots()
        im = plt.imshow(data_,origin='lower',norm= LogNorm( vmin=.1, vmax=1e0 ) )
        fig.colorbar(im)
        plt.show()
    
   
def get_qr_intensity( qr, data,vert_rect, show_roi=True ):
    
    V_K_label_array = roi.rectangles(vert_rect, data.shape)  #(y,x, hight, wdith)
    
    if show_roi:
        data_ = data*V_K_label_array 
        fig, ax = plt.subplots()
        im = plt.imshow(data_,origin='lower',norm= LogNorm( vmin=.1, vmax=1e0 ) )
        fig.colorbar(im)
        plt.show()    
        
    fig, ax = plt.subplots()
    for i, vr in enumerate( vert_rect):
        print (i, vr)
        V_K_label_array_i = roi.rectangles((vr,), data.shape)  #(y,x, hight, wdith)     
        roi_pixel_num = np.sum( V_K_label_array_i, axis=0)
        qr_ = qr  *V_K_label_array_i
        data_ = data*V_K_label_array_i  
    
        qr_ave = np.sum( qr_, axis=0)/roi_pixel_num
        data_ave = np.sum( data_, axis=0)/roi_pixel_num  
 
        ax.plot( qr_ave, data_ave,  label= 'interest_roi_%i'%i)
    ax.set_xlabel( r'$q_r$', fontsize=15)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend() 
    
   
#GiSAXS End
###############################

def show_label_array_on_image(ax, image, label_array, cmap=None,norm=None,
                              imshow_cmap='gray', **kwargs):  #norm=LogNorm(), 
    """
    This will plot the required ROI's(labeled array) on the image
    Additional kwargs are passed through to `ax.imshow`.
    If `vmin` is in kwargs, it is clipped to minimum of 0.5.
    Parameters
    ----------
    ax : Axes
        The `Axes` object to add the artist too
    image : array
        The image array
    label_array : array
        Expected to be an unsigned integer array.  0 is background,
        positive integers label region of interest
    cmap : str or colormap, optional
        Color map to use for plotting the label_array, defaults to 'None'
    imshow_cmap : str or colormap, optional
        Color map to use for plotting the image, defaults to 'gray'
    norm : str, optional
        Normalize scale data, defaults to 'Lognorm()'
    Returns
    -------
    im : AxesImage
        The artist added to the axes
    im_label : AxesImage
        The artist added to the axes
    """
    ax.set_aspect('equal')
    im = ax.imshow(image, cmap=imshow_cmap, interpolation='none',norm=LogNorm(norm), 
                   **kwargs)  #norm=norm,
    im_label = mpl_plot.show_label_array(ax, label_array, cmap=cmap, norm=norm,
                                **kwargs)  # norm=norm,
    
    
    return im, im_label    
    
    
 




import logging
import time

#import numpy as np

#from . import utils as core
#from . import roi





from lmfit import minimize, Model, Parameters

logger = logging.getLogger(__name__)


def multi_tau_auto_corr(num_levels, num_bufs, labels, images):
    ##comments, please add start_image, end_image, the default as None
    
    
    from skxray.core import roi
    from skxray.core  import utils as core
    
    """
    This function computes one-time correlations.
    It uses a scheme to achieve long-time correlations inexpensively
    by downsampling the data, iteratively combining successive frames.
    The longest lag time computed is num_levels * num_bufs.
    Parameters
    ----------
    num_levels : int
        how many generations of downsampling to perform, i.e.,
        the depth of the binomial tree of averaged frames
    num_bufs : int, must be even
        maximum lag step to compute in each generation of
        downsampling
    labels : array
        labeled array of the same shape as the image stack;
        each ROI is represented by a distinct label (i.e., integer)
    images : iterable of 2D arrays
        dimensions are: (rr, cc)
    Returns
    -------
    g2 : array
        matrix of normalized intensity-intensity autocorrelation
        shape (num_levels, number of labels(ROI))
    lag_steps : array
        delay or lag steps for the multiple tau analysis
        shape num_levels
    Notes
    -----
    The normalized intensity-intensity time-autocorrelation function
    is defined as
    :math ::
        g_2(q, t') = \frac{<I(q, t)I(q, t + t')> }{<I(q, t)>^2}
    ; t' > 0
    Here, I(q, t) refers to the scattering strength at the momentum
    transfer vector q in reciprocal space at time t, and the brackets
    <...> refer to averages over time t. The quantity t' denotes the
    delay time
    This implementation is based on code in the language Yorick
    by Mark Sutton, based on published work. [1]_
    References
    ----------
    .. [1] D. Lumma, L. B. Lurio, S. G. J. Mochrie and M. Sutton,
        "Area detector based photon correlation in the regime of
        short data batches: Data reduction for dynamic x-ray
        scattering," Rev. Sci. Instrum., vol 70, p 3274-3289, 2000.
    """
    # In order to calculate correlations for `num_bufs`, images must be
    # kept for up to the maximum lag step. These are stored in the array
    # buffer. This algorithm only keeps number of buffers and delays but
    # several levels of delays number of levels are kept in buf. Each
    # level has twice the delay times of the next lower one. To save
    # needless copying, of cyclic storage of images in buf is used.

    if num_bufs % 2 != 0:
        raise ValueError("number of channels(number of buffers) in "
                         "multiple-taus (must be even)")

    if hasattr(images, 'frame_shape'):
        # Give a user-friendly error if we can detect the shape from pims.
        if labels.shape != images.frame_shape:
            raise ValueError("Shape of the image stack should be equal to"
                             " shape of the labels array")

    # get the pixels in each label
    label_mask, pixel_list = roi.extract_label_indices(labels)

    num_rois = np.max(label_mask)

    # number of pixels per ROI
    num_pixels = np.bincount(label_mask, minlength=(num_rois+1))
    num_pixels = num_pixels[1:]

    if np.any(num_pixels == 0):
        raise ValueError("Number of pixels of the required roi's"
                         " cannot be zero, "
                         "num_pixels = {0}".format(num_pixels))

    # G holds the un normalized auto-correlation result. We
    # accumulate computations into G as the algorithm proceeds.
    G = np.zeros(((num_levels + 1)*num_bufs/2, num_rois),
                 dtype=np.float64)

    # matrix of past intensity normalizations
    past_intensity_norm = np.zeros(((num_levels + 1)*num_bufs/2, num_rois),
                                   dtype=np.float64)

    # matrix of future intensity normalizations
    future_intensity_norm = np.zeros(((num_levels + 1)*num_bufs/2, num_rois),
                                     dtype=np.float64)

    # Ring buffer, a buffer with periodic boundary conditions.
    # Images must be keep for up to maximum delay in buf.
    buf = np.zeros((num_levels, num_bufs, np.sum(num_pixels)),
                   dtype=np.float64)

    # to track processing each level
    track_level = np.zeros(num_levels)

    # to increment buffer
    cur = np.ones(num_levels, dtype=np.int64)

    # to track how many images processed in each level
    img_per_level = np.zeros(num_levels, dtype=np.int64)

    start_time = time.time()  # used to log the computation time (optionally)

    for n, img in enumerate(images):

        cur[0] = (1 + cur[0]) % num_bufs  # increment buffer

        # Put the image into the ring buffer.
        buf[0, cur[0] - 1] = (np.ravel(img))[pixel_list]

        # Compute the correlations between the first level
        # (undownsampled) frames. This modifies G,
        # past_intensity_norm, future_intensity_norm,
        # and img_per_level in place!
        _process(buf, G, past_intensity_norm,
                 future_intensity_norm, label_mask,
                 num_bufs, num_pixels, img_per_level,
                 level=0, buf_no=cur[0] - 1)

        # check whether the number of levels is one, otherwise
        # continue processing the next level
        processing = num_levels > 1

        # Compute the correlations for all higher levels.
        level = 1
        while processing:
            if not track_level[level]:
                track_level[level] = 1
                processing = False
            else:
                prev = 1 + (cur[level - 1] - 2) % num_bufs
                cur[level] = 1 + cur[level] % num_bufs

                buf[level, cur[level] - 1] = (buf[level - 1, prev - 1] +
                                              buf[level - 1,
                                                  cur[level - 1] - 1])/2

                # make the track_level zero once that level is processed
                track_level[level] = 0

                # call the _process function for each multi-tau level
                # for multi-tau levels greater than one
                # Again, this is modifying things in place. See comment
                # on previous call above.
                _process(buf, G, past_intensity_norm,
                         future_intensity_norm, label_mask,
                         num_bufs, num_pixels, img_per_level,
                         level=level, buf_no=cur[level]-1,)
                level += 1

                # Checking whether there is next level for processing
                processing = level < num_levels

    # ending time for the process
    end_time = time.time()

    logger.info("Processing time for {0} images took {1} seconds."
                "".format(n, (end_time - start_time)))

    # the normalization factor
    if len(np.where(past_intensity_norm == 0)[0]) != 0:
        g_max = np.where(past_intensity_norm == 0)[0][0]
    else:
        g_max = past_intensity_norm.shape[0]

    # g2 is normalized G
    g2 = (G[:g_max] / (past_intensity_norm[:g_max] *
                       future_intensity_norm[:g_max]))

    # Convert from num_levels, num_bufs to lag frames.
    tot_channels, lag_steps = core.multi_tau_lags(num_levels, num_bufs)
    lag_steps = lag_steps[:g_max]

    return g2, lag_steps


def _process(buf, G, past_intensity_norm, future_intensity_norm,
             label_mask, num_bufs, num_pixels, img_per_level, level, buf_no):
    """
    Internal helper function. This modifies inputs in place.
    This helper function calculates G, past_intensity_norm and
    future_intensity_norm at each level, symmetric normalization is used.
    Parameters
    ----------
    buf : array
        image data array to use for correlation
    G : array
        matrix of auto-correlation function without
        normalizations
    past_intensity_norm : array
        matrix of past intensity normalizations
    future_intensity_norm : array
        matrix of future intensity normalizations
    label_mask : array
        labels of the required region of interests(roi's)
    num_bufs : int, even
        number of buffers(channels)
    num_pixels : array
        number of pixels in certain roi's
        roi's, dimensions are : [number of roi's]X1
    img_per_level : array
        to track how many images processed in each level
    level : int
        the current multi-tau level
    buf_no : int
        the current buffer number
    Notes
    -----
    :math ::
        G   = <I(\tau)I(\tau + delay)>
    :math ::
        past_intensity_norm = <I(\tau)>
    :math ::
        future_intensity_norm = <I(\tau + delay)>
    """
    img_per_level[level] += 1

    # in multi-tau correlation other than first level all other levels
    #  have to do the half of the correlation
    if level == 0:
        i_min = 0
    else:
        i_min = num_bufs//2

    for i in range(i_min, min(img_per_level[level], num_bufs)):
        t_index = level*num_bufs/2 + i

        delay_no = (buf_no - i) % num_bufs

        past_img = buf[level, delay_no]
        future_img = buf[level, buf_no]

        #  get the matrix of auto-correlation function without normalizations
        tmp_binned = (np.bincount(label_mask,
                                  weights=past_img*future_img)[1:])
        G[t_index] += ((tmp_binned / num_pixels - G[t_index]) /
                       (img_per_level[level] - i))

        # get the matrix of past intensity normalizations
        pi_binned = (np.bincount(label_mask,
                                 weights=past_img)[1:])
        past_intensity_norm[t_index] += ((pi_binned/num_pixels
                                         - past_intensity_norm[t_index]) /
                                         (img_per_level[level] - i))

        # get the matrix of future intensity normalizations
        fi_binned = (np.bincount(label_mask,
                                 weights=future_img)[1:])
        future_intensity_norm[t_index] += ((fi_binned/num_pixels
                                           - future_intensity_norm[t_index]) /
                                           (img_per_level[level] - i))

    return None  # modifies arguments in place!

def interp_zeros(  data ): 
    from scipy.interpolate import interp1d
    gf = data.ravel() 
    indice, = gf.nonzero() 
    start, stop = indice[0], indice[-1]+1 
    dx,dy = data.shape 
    x=np.arange( dx*dy ) 
    f = interp1d(x[indice], gf[indice]) 
    gf[start:stop] = f(x[start:stop]) 
    return gf.reshape([dx,dy]) 
 

 