def graph_failed_objects(fore_mask,ra,dec,tilename,nested):
    '''
    Function to generate plots given a foreground star mask, RAs, Decs, and Flags for a tile. 
    Also specify if the HEALPix map is nested or not.
    Assumes that the flags have already been applied and objects being passed are ones that failed.
    
    Variables:
    fore_mask - path to your foreground star mask
    ra - array of RAs for failed objects
    dec - array of Decs for failed objects
    tilename - DECam Tilename
    nested - is the foreground mask nested or not
    '''
    import os
    import healpy as hp
    import healsparse as hsp
    import matplotlib.pyplot as plt
    import numpy as np
    import skyproj
    
    STOR_DIR = '/data/des81.a/data/kherron/projects/LSS_Data/STORE.fits'

    mask = hp.read_map(fore_mask)
    nside = hp.npix2nside(len(mask))
    print(nside)
    #Mask out everything but our tile and change sentinel values
    ras, decs = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)),lonlat=True)
    sel = (ras < min(ra)) | (ras > max(ra)) | (decs < min(dec)) | (decs > max(dec))
    mask[sel]=hp.UNSEEN
    mask[mask == 0] = hp.UNSEEN
    
    #Write changes and generate healsparse map
    hp.write_map(STOR_DIR,mask,overwrite=True,nest=False)
    hspmap = hsp.HealSparseMap.read(STOR_DIR,
                                   nside_coverage=nside)
    
    #Get our foreground array
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    sp = skyproj.Skyproj(ax=ax)
    im,ras,decs,val = sp.draw_hspmap(hspmap,cmap='jet',alpha=0.8)
    #plt.show()
    plt.clf()
    
    #Get our 2D Hist array
    _ = plt.figure(figsize=(8,8))
    sel = (ras > min(ra)) & (ras < max(ra)) & (decs > min(dec)) & (decs < max(dec))
    pix = []
    for i in range(len(ra)):
        pix = np.append(pix,hp.ang2pix(nside,ra[i],dec[i],lonlat=True))
        
    mask = hp.read_map(fore_mask,nest=nested)
    
    indices = []
    for i in range(len(pix)):
        indices = np.append(indices, mask[int(pix[i])]==0)
    indices=indices.astype(bool)
    fig3 = plt.hist2d(ra[indices],dec[indices],bins=np.shape(val)[0],norm=LogNorm(),cmap='jet')
    plt.gca().invert_xaxis()
    #plt.show()
    plt.clf()
    
    #Finally plot it all together
    fig = plt.figure(figsize=(10,10))
    
    plt.title(tilename)
    
    sel = fig3[0] != 0
    ma = np.ma.MaskedArray(data=fig3[0],mask =sel )    
    plt.imshow(np.flip(np.rot90(ma),axis=1),cmap='jet')
    plt.imshow(np.flip(np.flip(val,axis=1),axis=0),alpha=0.7,cmap='jet',norm=LogNorm())
    
    plt.xticks(ticks=[])
    plt.yticks(ticks=[])
    plt.show()
