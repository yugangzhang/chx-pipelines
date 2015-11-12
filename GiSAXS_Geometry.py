def make_gisaxs_grid( qp_w= 10, qz_w = 12, dim_p =100,dim_z=120):
    y, x = np.indices( [dim_z,dim_p] )
    Np = int(dim_p/qp_w)
    Nz = int(dim_z/qz_w)
    noqs = Np*Nz
    ind = 1
    for i in range(0,Np):
        for j in range(0,Nz):        
            y[ qp_w*i: qp_w*(i+1), qz_w*j:qz_w*(j+1)]=  ind
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
    print ('incident_angle (alphai) is: %s'%(alphai* 180/np.pi))
    px,py = pixelsize
    y, x = np.indices( [dimy,dimx] )    
    alphaf = np.arctan2( (y-inc_y0)*py*10**(-6), Lsd ) - alphai 
    thetaf = np.arctan2( (x-inc_x0)*px*10**(-6), Lsd ) - thetai   
    
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
    qp = np.sqrt( qx**2 + qy**2 ) 
    
    
    return qx*pref  , qy*pref  , qp*pref  , qz*pref 

