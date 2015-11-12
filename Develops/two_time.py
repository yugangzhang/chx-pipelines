import numpy as np
import sys
import time
import skxray.core.roi as roi



def autocor_two_time( num_buf,  ring_mask, imgs, num_lev=None, start_img=None, end_img=None    ):
    

    #print (dly)
    if start_img is None:start_img=0
    if end_img is None:
        try:
            end_img= len(imgs)
        except:
            end_img= imgs.length
            
    #print (start_img, end_img)    
    noframes = end_img - start_img +1
    
    if num_lev is None:num_lev = int(np.log( noframes/(num_buf-1))/np.log(2) +1) +1
    print ( 'The lev number is %s'%num_lev)
    
    dly, dict_dly = delays( num_lev, num_buf, time=1 )
    #print (dly.max())
    
    qind, pixelist = roi.extract_label_indices(   ring_mask  )
    noqs = np.max(qind)    
    nopr = np.bincount(qind, minlength=(noqs+1))[1:]
    nopixels = nopr.sum() 
    
    start_time = time.time()
    
    buf=np.zeros([num_lev,num_buf,nopixels])  #// matrix of buffers, for store img
    
    
    cts=np.zeros(num_lev)
    cur=np.ones(num_lev) * num_buf
    countl = np.array( np.zeros(  num_lev ),dtype='int')  
    
    g12 =  np.zeros( [ noframes, noframes, noqs] )      
    
    num= np.array( np.zeros(  num_lev ),dtype='int')          
    time_ind ={key: [] for key in range(num_lev)}   
    
    ttx=0        
    for n in range( start_img, end_img ):   ##do the work here
        
        cur[0]=1+cur[0]%num_buf  # increment buffer  
        img = imgs[n] 
        
        #print ( 'The insert image is %s' %(n) )
    
        buf[0, cur[0]-1 ]=  (np.ravel(img))[pixelist]
        img=[] #//save space 
        countl[0] = 1+ countl[0]
        current_img_time = n - start_img +1
    
        process_two_time(lev=0, bufno=cur[0]-1,n=current_img_time,
                        g12=g12, buf=buf, num=num, num_buf=num_buf, noqs=noqs, qind=qind, nopr=nopr, dly=dly)     
        time_ind[0].append(  current_img_time   )
        processing=1
        lev=1
        while processing:
            if cts[lev]:
                prev=  1+ (cur[lev-1]-1-1+num_buf)%num_buf
                cur[lev]=  1+ cur[lev]%num_buf
                countl[lev] = 1+ countl[lev]                                
                buf[lev,cur[lev]-1] = ( buf[lev-1,prev-1] + buf[lev-1,cur[lev-1]-1] ) /2.
                cts[lev]=0                
                t1_idx=   (countl[lev]-1) *2
                current_img_time = ((time_ind[lev-1])[t1_idx ] +  (time_ind[lev-1])[t1_idx +1 ] )/2. 
                time_ind[lev].append(  current_img_time      )  
                process_two_time(lev=lev, bufno=cur[lev]-1,n=current_img_time,
                        g12=g12, buf=buf, num=num, num_buf=num_buf, noqs=noqs, qind=qind, nopr=nopr, dly=dly)  
                lev+=1
                #//Since this level finished, test if there is a next level for processing
                if lev<num_lev:processing = 1
                else:processing = 0                                
            else:
                cts[lev]=1      #// set flag to process next time
                processing=0    #// can stop until more images are accumulated              
 
        
        if  n %(noframes/10) ==0:
            sys.stdout.write("#")
            sys.stdout.flush()                
    
    
    for q in range(noqs):            
        x0 =  g12[:,:,q]
        g12[:,:,q] = np.tril(x0) +  np.tril(x0).T - np.diag( np.diag(x0) )            
    elapsed_time = time.time() - start_time
    print ('Total time: %.2f min' %(elapsed_time/60.))
    
    
    return g12, elapsed_time/60.



    
    
    
            
def process_two_time(lev, bufno,n ,    
                     g12, buf, num, num_buf,noqs,qind,nopr, dly ):
    num[lev]+=1  
    if lev==0:imin=0
    else:imin= int(num_buf/2 )
    for i in range(imin, min(num[lev],num_buf) ):
        ptr=lev*int(num_buf/2)+i    
        delayno=(bufno-i)%num_buf #//cyclic buffers            
        IP=buf[lev,delayno]
        IF=buf[lev,bufno]
        I_t12 =  (np.histogram(qind, bins=noqs, weights= IF*IP))[0]
        I_t1  =  (np.histogram(qind, bins=noqs, weights= IP))[0]
        I_t2  =  (np.histogram(qind, bins=noqs, weights= IF))[0]
        tind1 = (n-1)
        tind2=(n -dly[ptr] -1)
        
        if not isinstance( n, int ):                
            nshift = 2**(lev-1)                
            for i in range( -nshift+1, nshift +1 ):
                #print tind1+i
                g12[ int(tind1 + i), int(tind2 + i) ] =I_t12/( I_t1 * I_t2) * nopr
        else:
                #print tind1
            g12[ tind1, tind2 ]  =   I_t12/( I_t1 * I_t2) * nopr       
        
        
        

    
def delays( num_lev=3, num_buf=4, time=1 ): 
    ''' DOCUMENT delays(time=)
        return array of delays.
        KEYWORD:  time: scale delays by time ( should be time between frames)
     '''
    if num_buf%2!=0:print ("nobuf must be even!!!"    )
    dly=np.zeros( (num_lev+1)*int(num_buf/2) +1  )        
    dict_dly ={}
    for i in range( 1,num_lev+1):
        if i==1:imin= 1
        else:imin= int(num_buf/2)+1
        ptr=(i-1)*int(num_buf/2)+ np.arange(imin,num_buf+1)
        dly[ptr]= np.arange( imin, num_buf+1) *2**(i-1)            
        dict_dly[i] = dly[ptr-1]            
        dly*=time
        #print (i, ptr, imin)
    return dly, dict_dly
            
     
    

class Get_Pixel_Array(object):
    def __init__(self, indexable, pixelist):
        self.indexable = indexable
        self.pixelist = pixelist
        self.shape = indexable.shape
        try:
            self.length= len(indexable)
        except:
            self.length= indexable.length           
            
    def get_data(self ): 
        data_array = np.zeros([ self.length,len(self.pixelist)])
        for key in range(self.length):
            data_array[key] = np.ravel( self.indexable[key])[self.pixelist]  
        return data_array
    
    

    
 
    
    
    
def autocor_arrays_two_time( seg1, pixelist,qind, seg2=None,                         
                            get_half=False,get_whole=False,up_half=True,
            print_=True,    ):
        #for same seg1 and seg2, use get_half=True, get_whole=True

        if seg2 is None:
            seg2=seg1 #half_flag=True #half_flag for only calculate half

        start_time = time.time()    
        m,n = seg1.shape
        #print m,n
        noqs = len( np.unique(qind) )
        nopr = np.bincount(qind, minlength=(noqs+1))[1:] #qind start from 1   
        seg1f = np.ravel( seg1 )                
        qinds = np.zeros( [ m, len(pixelist)],dtype=int)
        for i in range( m ):qinds[i]=qind + ( max(qind) ) *i     #qind start from 1   
        G_t1 = np.bincount( qinds.ravel(), weights = seg1f )[1:]      #qind start from 1          
        G_t1= G_t1.reshape( [m, noqs] )      
        
        g12s =  np.zeros( [m,m, noqs] )
        cal=-1
        
        for i in range( m ):
            if get_half:
                termin=i+1;
            else:
                termin=m
                #print ('not get half')
            if up_half:
                #print ('cal up_half')
                seg12i = ( seg1[:termin] * seg2[i]).ravel()
                qindsi=qinds[:termin]
            else: #for down_half
                seg12i = ( seg1[termin-1 : ] * seg2[i]).ravel()
                qindsi=qinds[: m- (termin-1) ]
            G_t12i = np.bincount( qindsi.ravel(), weights = seg12i )[1:]
            #print G_t12i.shape, qindsi.flatten().shape, seg12i.shape
            #print seg12i.shape, G_t12i.shape            
            mG = (G_t12i.shape)[0]
            #print seg1[ termin-1  : ].shape            
            G_t12i= G_t12i.reshape( [ int(mG/noqs), noqs] )            
            seg2i =  np.bincount( qind.ravel(), weights = seg2[i] )[1:]

            #print seg12i.shape,qindsi.shape, mG,G_t12i.shape,seg2i.shape
            
            #print G_t12i.shape            
            #print G_t12i.shape, seg2i.shape, G_t1.shape, termin, g12s.shape           
            #print termin-1
            if up_half:
                g12s[i][:termin] = G_t12i/( G_t1[:termin] * seg2i ) * nopr
            else:
                g12s[i][termin-1:] = G_t12i/( G_t1[termin-1:] * seg2i ) * nopr
            
            if print_:
                if  int(i /(m/10.)) >cal:                
                    sys.stdout.write("#")
                    sys.stdout.flush()
                    cal = int(i /(m/10.))
                    #elapsed_time = time.time() - start_time
                    #print ('%s:  cal time: %.2f min' %( int(i /(m/10.)),elapsed_time/60.)  )
                    
##        if up_half==False:
##            g12s_=g12s.copy()
##            g12s=g12s*0
##            for q in range(noqs):            
##                x0 =  g12s_[:,:,q]                
##                g12s[:,:,q] =   tril(x0).T 
                
            
        if get_whole:
            for q in range(noqs):  
                #print (q)
                x0 =  g12s[:,:,q]
                if up_half:
                    g12s[:,:,q] = np.tril(x0) +  np.tril(x0).T - np.diag(np.diag(x0))
                else:
                    g12s[:,:,q] = np.triu(x0) +  np.triu(x0).T - np.diag(np.diag(x0))
        if print_:
            elapsed_time = time.time() - start_time
            print ('Total time: %.2f min' %(elapsed_time/60.)  )
            
        return g12s



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
 

def autocor_large_arrays_two_time( data, pixelist, qind, divide=4,  
                 only_one_nodia=True,  get_whole=True,print_=True  ): 
    #currently, data is the pims data
    
    #noframes/divde should be int
    try:
        noframes = len(data)
    except:
        noframes = data.length    
    if noframes/divide - int( noframes/divide ) != 0:
        print ('noframes/divde should be int!!! Please give another divide number!!!')          
    start_time = time.time() 
    noqs = len( np.unique(qind) )
    g12L =  np.zeros( [noframes,noframes, noqs] )
    step_fram = int(noframes/divide)    
    data_div = np.zeros( [step_fram, len(pixelist), divide])    
    for block in range( divide):
        imgs_ = data[block*step_fram : (block+1)* step_fram]
        #print (imgs_.shape)
        imgsr = Reverse_Coordinate(imgs_, mask=None)
        data_div[:,:,block]= Get_Pixel_Array( imgsr, pixelist).get_data()
    
    #print (data_div.shape)
        
    cal=-1
    m=0
    st = divide**2/20.
    
    for block1 in range(divide):
        data1 = data_div[:,:,block1]
        for block2 in range(block1, divide):
            #print (block1, block2)
            if print_:
                m+=1
                if  int( m/st ) >cal:                
                    sys.stdout.write("#")
                    sys.stdout.flush()
                    cal = int( m/st  ) 
                        
            if block1==block2: #this is the diagonal part
                #print (block1, block2)
                fm1,fm2 = [block1*step_fram , (block1+1)* step_fram]
                #print fm1,fm2
                g12L[fm1:fm2,fm1:fm2,:]= autocor_arrays_two_time(                        
                        seg1 = data1, pixelist=pixelist,qind=qind, seg2=None, 
                        get_half=True,get_whole=False,up_half=True,
                        print_= False)  

            else:  #this is the no-diagonal part              
                if not only_one_nodia:  #cal all nodiagon
                    data2 = data_div[:,:,block2]
                    fm1,fm2 = [block1*step_fram , (block1+1)* step_fram]
                    fm3,fm4 = [block2*step_fram , (block2+1)* step_fram]
                    g12L[fm3:fm4,fm1:fm2,:] =autocor_arrays_two_time(    
                        seg1=data1, pixelist=pixelist,qind=qind,seg2=data2,
                        get_half=False, get_whole=False,  up_half=True,                      
                       print_= False)
                else: # cal only one nodiagon
                    block2x=block1+1
                    data2 = data_div[:,:,block2x]
                    fm1,fm2 = [block1*step_fram , (block1+1)* step_fram]
                    fm3,fm4 = [block2x*step_fram , (block2x+1)* step_fram]
                    g12L[fm3:fm4,fm1:fm2,:] =autocor_arrays_two_time(
                        seg1=data1,pixelist=pixelist,qind=qind,seg2=data2,
                        get_half=True,get_whole=False,up_half=False,
                       print_= False)
                    
    if get_whole:
        for q in range(noqs):  
            #print (q)
            x0 =  g12L[:,:,q]
            g12L[:,:,q] = np.tril(x0) +  np.tril(x0).T - np.diag(np.diag(x0))                   
    if print_:
        elapsed_time = time.time() - start_time
        print ('Total time: %.2f min' %(elapsed_time/60.)  )

    return g12L












def test():
    pixelist_q1 = pixelist[ np.where( qind ==1)[0] ]
    seg_q1 =   Get_Pixel_Array( imgsr, pixelist_q1).get_data()
    nopr_q1 = len(pixelist_q1)
    sum1 = (np.average( seg_q1, axis=1)).reshape( 1, seg_q1.shape[0]   )  
    sum2 = sum1.T
    m= np.dot(   seg_q1, seg_q1.T)  /sum1  / sum2  / nopr_q1


















