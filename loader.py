import sys
import os
import numpy as np
from misc import make_path, default_maxedge, default_minpers, headfolder
import warnings
warnings.filterwarnings("ignore", category=UserWarning)



'''
Script to load in various datasets whose topology we want to examine,
as well as topological data (e.g. loops) we have computed on these.

NB: both raw and filtered data should, when loaded, have shape (D, time),
    where D is the dimension of the system considered.

'''

#The main load function: enter new datasets here with an appropriate name
def load_raw(name, useall=False):

    datadir = headfolder + '/Data/Raw' 

    if name == 'Lorenz63':
        fname = '%s/Lorenz63/lorenz_10k_raw.npy' % datadir
        if useall:
            data = (np.load(fname).T)
        else:
            data = (np.load(fname).T)[:,10000:10000+100000]
    
    elif name == 'Lorenz96':
        fname = '%s/Lorenz96/lorenz96_raw.txt' % datadir
        if useall:
            data= np.loadtxt(fname)[:,:200000]
        else:
            data= np.loadtxt(fname)[:,:20000]

    elif name == 'CDV':
        fname = '%s/CDV/cdv_deterministic.np' % datadir
        if useall:
            data= np.fromfile(fname).reshape([6,-1])
        else:
            data= np.fromfile(fname).reshape([6,-1])[:,10000:10000+40000]
    
    elif name == 'JetLat':
        lats = np.loadtxt('%s/JetLat/jetlatitude_ERA20C_djf_1901-2010.txt' % datadir)
        pc0 = np.loadtxt('%s/JetLat/pc0_ua850_ERA20C_1901-2010.txt' % datadir)
        pc1 = np.loadtxt('%s/JetLat/pc1_ua850_ERA20C_1901-2010.txt' % datadir)
        lats /= lats.std()
        pc0 /= pc0.std()
        pc1 /= pc1.std()
        data = np.vstack([lats, pc0, pc1])
   
    elif name == 'Gaussian':
        fname = '%s/Gaussian/gaussian_raw.txt' % datadir
        data = np.loadtxt(fname)

    else:
        sys.exit("(load_raw) ERROR: unrecognized dataset %s" % name)

    return data



#For loading filtered data
def load_filt(name, method='GaussianKDE', identifier=None,\
              num_bins=160, perc_to_keep=70,\
              reverse_order=False, eofs=False, normalise=False):
        

    if method == 'timeseries':
        data = load_raw(name)
        from eofs.standard import Eof
        EOFs=Eof(data.T)
        pcdata=EOFs.pcs(npcs=num_pcs).T
        data = pcdata[:,:10000]

    else:
        if 'Perc' in name:
            orig_name, percstr = name.split('-')
            perc_to_keep = int(percstr.split('Perc')[-1])

        else:
            orig_name = name

        if eofs:
            method += '_EOFs'

        if perc_to_keep == 'XXX':
            data = load_raw(orig_name)
        else:
            datadir = headfolder + './Data/Filtered'
            if identifier is None:
                if reverse_order:
                    fname = '%s/%s/%s_filtered_leastdense_%s_percent_%sbins_%s.txt' % (datadir, orig_name, orig_name, perc_to_keep, num_bins, method)
                else:
                    fname = '%s/%s/%s_filtered_densest_%s_percent_%sbins_%s.txt' % (datadir, orig_name, orig_name, perc_to_keep, num_bins, method)
            else:
                if reverse_order:
                    fname = '%s/%s/%s_filtered_densest_%s_percent_%sbins_%s_%s.txt' % (datadir, orig_name, orig_name, perc_to_keep, num_bins, method, identifier)
                else:
                    fname = '%s/%s/%s_filtered_leastdense_%s_percent_%sbins_%s_%s.txt' % (datadir, orig_name, orig_name, perc_to_keep, num_bins, method, identifier)

            data = np.loadtxt(fname)

    if normalise:
        num_dims = data.shape[0]
        for n in range(num_dims):
            data[n,:] = data[n,:]/data[n,:].std()
  

    return data



#Loading the additionally sparsified data using gudhi's pre_sparse argument
def load_from_pers(name, max_edge=None, min_pers=None, sparse=0.7, pre_sparse=None,\
                   method='GaussianKDE', num_bins=160, reverse_order=False, eofs=False):

   
    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    else:
        basename = name
        perc_to_keep = None


    if perc_to_keep is None:
        if reverse_order:
            outpath_root = headfolder + './Data/Processed/Loops/%s/reversed' % basename
        else:
            outpath_root = headfolder + './Data/Processed/Loops/%s/standard' % basename
    else:
        if reverse_order:
             outpath_root = headfolder + './Data/Processed/Loops/%s/reversed/%s-Perc%s' % (basename, basename, perc_to_keep)
        else:
             outpath_root = headfolder + './Data/Processed/Loops/%s/standard/%s-Perc%s' % (basename, basename, perc_to_keep)



    filtdata = load_filt(name, num_bins=num_bins, method=method, reverse_order=reverse_order, eofs=eofs)
    
    if (max_edge is None):
        max_edge = default_maxedge(filtdata)

    if min_pers is None:
        min_pers = default_minpers(filtdata, name)

    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)

    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)

    
    if eofs:
        fname = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_f.txt' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
    else:
        fname = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_f.txt' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)

    data = []
    with open(fname,'r') as f:
        for i,l in enumerate(f.readlines()):
            if i == 0:
                nele = ([int(j) for j in l.split(" ")][-1])
            if i > 0 and i <=nele:
                data.append([float(j) for j in l.split(" ")])
    data = np.array(data).T

    return data



#For loading loops produced with PersLoop
def load_loops(name, max_edge=None, min_pers=None, sparse=0.7, pre_sparse=None,\
               method='GaussianKDE', reverse_order=False, eofs=False):
    
    import phtools as ph
    
    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    else:
        basename = name
        perc_to_keep = None


    if perc_to_keep is None:
        if reverse_order:
            outpath_root = headfolder + './Data/Processed/Loops/%s/reversed' % basename
        else:
            outpath_root = headfolder + './Data/Processed/Loops/%s/standard' % basename
    else:
        if reverse_order:
             outpath_root = headfolder + './Data/Processed/Loops/%s/reversed/%s-Perc%s' % (basename, basename, perc_to_keep)
        else:
             outpath_root = headfolder + './Data/Processed/Loops/%s/standard/%s-Perc%s' % (basename, basename, perc_to_keep)
 

    filtdata = load_filt(name, method=method, reverse_order=reverse_order, eofs=eofs)

    if (max_edge is None):
        max_edge = default_maxedge(filtdata)

    if min_pers is None:
        min_pers = default_minpers(filtdata, name)


    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)

    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)

    if eofs:
        fname = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_floops' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
    else:
        fname = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_floops' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)

    loop_numbers = ph.loopnums(fname)
    loop_list = []
    for l in loop_numbers:
        tmp=ph.getloop(l, fname)
        loop_list.append(tmp)
        
    return loop_list



#Loading components
def load_comps(name, num2keep=5,\
                     max_edge=None, min_pers=None,\
                     sparse=0.7, pre_sparse=None,\
                     method='GaussianKDE', reverse_order=False,\
                     eofs=False):

    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    else:
        basename = name
        perc_to_keep = None


    if perc_to_keep is None:
        if reverse_order:
            outpath_root = headfolder + './Data/Processed/Loops/%s/reversed' % basename
        else:
            outpath_root = headfolder + './Data/Processed/Loops/%s/standard' % basename
    else:
        if reverse_order:
             outpath_root = headfolder + './Data/Processed/Loops/%s/reversed/%s-Perc%s' % (basename, basename, perc_to_keep)
        else:
             outpath_root = headfolder + './Data/Processed/Loops/%s/standard/%s-Perc%s' % (basename, basename, perc_to_keep)
   

    filtdata = load_filt(name, method=method, reverse_order=reverse_order, eofs=eofs)
  
    if (max_edge is None):
        max_edge = default_maxedge(filtdata)

    if min_pers is None:
        min_pers = default_minpers(filtdata, name)


    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)

    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)

    if eofs:
        fname = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-comps.npy' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str,\
                                                                                      presparse_str, num2keep)
    else:
        fname = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-comps.npy' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str,\
                                                                                 presparse_str, num2keep)


    components = np.load(fname, allow_pickle=True)

    return components


#Loading deathtimes of components
def load_comp_deaths(name, num2keep=5,\
                     max_edge=None, min_pers=None,\
                     sparse=0.7, pre_sparse=None,\
                     method='GaussianKDE', reverse_order=False,\
                     eofs=False):

    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    else:
        basename = name
        perc_to_keep = None


    if perc_to_keep is None:
        if reverse_order:
            outpath_root = headfolder + './Data/Processed/Loops/%s/reversed' % basename
        else:
            outpath_root = headfolder + './Data/Processed/Loops/%s/standard' % basename
    else:
        if reverse_order:
             outpath_root = headfolder + './Data/Processed/Loops/%s/reversed/%s-Perc%s' % (basename, basename, perc_to_keep)
        else:
             outpath_root = headfolder + './Data/Processed/Loops/%s/standard/%s-Perc%s' % (basename, basename, perc_to_keep)
    
    
    filtdata = load_filt(name, method=method, reverse_order=reverse_order, eofs=eofs)
  
    if (max_edge is None):
        max_edge = default_maxedge(filtdata)

    if min_pers is None:
        min_pers = default_minpers(filtdata, name)


    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)

    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)

    if eofs:
        fname = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-comps_deathtimes.txt' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                                 sparse_str, presparse_str, num2keep)
        fname_sizes = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-comps_sizes.npy' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                                  sparse_str, presparse_str, num2keep)
    else:
        fname = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-comps_deathtimes.txt' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                            sparse_str, presparse_str, num2keep)
        fname_sizes = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-comps_sizes.npy' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                            sparse_str, presparse_str, num2keep)


    death_times = np.loadtxt(fname)[:num2keep]

    #Replace 0's with min_pers
    death_times[death_times == 0.] = min_pers

    #Get sizes of components, padding out with zeros
    sizes = np.zeros(num2keep)
    sizes_real = np.load(fname_sizes, allow_pickle=True)
    for n in range(len(sizes_real)):
        sizes[n] = sizes_real[n]

    return death_times, sizes


#Loading birth/death-times of loops
def load_loop_lives(name, num2keep=5,\
                     max_edge=None, min_pers=None,\
                     sparse=0.7, pre_sparse=None,\
                     method='GaussianKDE', reverse_order=False,
                     eofs=False):

    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    else:
        basename = name
        perc_to_keep = None


    if perc_to_keep is None:
        if reverse_order:
            outpath_root = headfolder + './Data/Processed/Loops/%s/reversed' % basename
        else:
            outpath_root = headfolder + './Data/Processed/Loops/%s/standard' % basename
    else:
        if reverse_order:
             outpath_root = headfolder + './Data/Processed/Loops/%s/reversed/%s-Perc%s' % (basename, basename, perc_to_keep)
        else:
             outpath_root = headfolder + './Data/Processed/Loops/%s/standard/%s-Perc%s' % (basename, basename, perc_to_keep)
   

    filtdata = load_filt(name, method=method, reverse_order=reverse_order, eofs=eofs)
  
    if (max_edge is None):
        max_edge = default_maxedge(filtdata)

    if min_pers is None:
        min_pers = default_minpers(filtdata, name)


    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)

    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)
   

    if eofs:
        fname_births = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-loops_birthtimes.txt' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                                        sparse_str, presparse_str, num2keep)
        fname_deaths = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-loops_deathtimes.txt' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                                        sparse_str, presparse_str, num2keep)
    else:
        fname_births = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-loops_birthtimes.txt' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                                   sparse_str, presparse_str, num2keep)
        fname_deaths = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_%s-oldest-loops_deathtimes.txt' % (outpath_root, name, method, max_edge_str, min_pers_str,\
                                                                                                   sparse_str, presparse_str, num2keep)
   

    birth_times = np.loadtxt(fname_births)
    death_times = np.loadtxt(fname_deaths)

    if birth_times.shape:
        births = birth_times[:num2keep]
    else:
        births = np.array([birth_times])

    if death_times.shape:
        deaths = death_times[:num2keep]
    else:
        deaths = np.array([death_times])

    return births, deaths



