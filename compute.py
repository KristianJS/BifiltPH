import sys
import os
import argparse
import numpy as np
import phtools as ph
from shutil import copyfile
from loader import load_raw, load_filt
from misc import make_path, default_maxedge, default_minpers,headfolder




'''
Script to do all the topological/geometric computations. Includes barcodes, loops
and any additional local density/dimension computations we might want to add.
'''



#Function to compute barcodes as well as loops computed with PersLoop
#Assumed to work on filtered data. The optional args are:
#
# * max_edge   - length up-to-which you want to build the filtration
#                (if None, a sensible choice is made)
# * min_pers   - only calculate barcodes etc for objects surviving a distance greater than this
#                (if None, a sensible choice is made)
# * sparse     - See https://gudhi.inria.fr/python/latest/rips_complex_user.html,.
#                If you don't use this we end up with too much data. (edited) 
#                (default is 0.9 based on experience)
#
# * pre_sparse - Another way to sparsify dataset without changing the topology.
#
# * remaining  - The three remaining optargs are unlikely to change at the moment.
#
def compute_homology(name,\
                     max_edge=None,\
                     min_pers=None,\
                     sparse=0.7,\
                     pre_sparse=None,\
                     method='GaussianKDE',\
                     identifier=None,\
                     num_bins=120,\
                     reverse_order=False,\
                     eofs=False,\
                     perc_to_keep=None,\
                     num_comp_to_keep=5,\
                     use_pl=True):
    
    
    print(100*'=')
    print("(compute_homology) Computing topological cocycles for %s" % name)
    if eofs:
        print("(compute_homology) Using EOFs: Yes")
    else:
        print("(compute_homology) Using EOFs: No")
    if not(use_pl):
        print("(compute_homology) NB: Omitting PersLoop computations. Working with Gudhi only.")


    #Load filtered data
    print("(compute_homology) Name = %s" % name)
    data = load_filt(name, method=method, identifier=identifier, num_bins=num_bins,\
                           perc_to_keep=perc_to_keep, reverse_order=reverse_order, eofs=eofs)

    data = data[:3,:] #restrict to 3D in case it's higher

    if 'Perc' in name:
            name, percstr = name.split('-')
            perc_to_keep = percstr.split('Perc')[-1]

    print("(compute_homology) Data shape is %s" % (data.shape,))
    
    #Set the output path
    if perc_to_keep is None:
        if reverse_order:
            outpath_root = headfolder + './Data/Processed/Loops/%s/reversed' % name
            outpath_root_figs = headfolder + './Figures/threeD/%s/reversed' % name
        else:
            outpath_root = headfolder + './Data/Processed/Loops/%s/standard' % name
            outpath_root_figs = headfolder + './Figures/threeD/%s/standard' % name
    else:
        if reverse_order:
             outpath_root = headfolder + './Data/Processed/Loops/%s/reversed/%s-Perc%s' % (name, name, perc_to_keep)
             outpath_root_figs = headfolder + './Figures/threeD/%s/reversed/%s-Perc%s' % (name, name, perc_to_keep)
        else:
             outpath_root = headfolder + './Data/Processed/Loops/%s/standard/%s-Perc%s' % (name, name, perc_to_keep)
             outpath_root_figs = headfolder + './Figures/threeD/%s/standard/%s-Perc%s' % (name, name, perc_to_keep)
    make_path(outpath_root)
    make_path(outpath_root_figs)
    
    
    #Work out sensible default options of max_edge and min_pers if one/both is not specified
    if (max_edge is None):
        
        print("(compute_homology) max_edge not specified. Generating sensible defaults...")
        max_edge = default_maxedge(data)

    if min_pers is None:
            
        print("(compute_homology) min_pers not specified. Generating sensible defaults...")
        min_pers = default_minpers(data, name)
      
    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)
    
    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)

    #Determination of the filename
    if perc_to_keep is None:
        if eofs:
            fname = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
            fname2 = '%s/%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root_figs, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
        else:
            fname = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
            fname2 = '%s/%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root_figs, name, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
    else:
        if eofs:
            fname = '%s/%s-Perc%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root, name, perc_to_keep, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
            fname2 = '%s/%s-Perc%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root_figs, name, perc_to_keep, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
        else:
            fname = '%s/%s-Perc%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root, name, perc_to_keep, method, max_edge_str, min_pers_str, sparse_str, presparse_str)
            fname2 = '%s/%s-Perc%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_' % (outpath_root_figs, name, perc_to_keep, method, max_edge_str, min_pers_str, sparse_str, presparse_str)


    #Do the computations
    print("(compute_homology) Asking for the %s longest-lived components" % num_comp_to_keep)
    print("(compute_homology) Using (max_edge, min_pers, sparse, pre_sparse) = (%s, %s, %s, %s)" % (max_edge, min_pers, sparse, pre_sparse))
    print("(compute_homology) Entering phtools...")
    print(100*'-')
    ph.persloop(data=data.T, nme=fname, num_comp_to_keep=num_comp_to_keep, max_edge=max_edge, min_pers=min_pers, sparse=sparse,
                pre_sparse=pre_sparse,tkbackend=tkback,use_pl=use_pl)
    print(100*'-')
    
    #Signal successful computation!
    print("(compute_homology) Success! Output stored in %s" % outpath_root)
    print(100*'=')

    #Copy the birthdeath plot to Figures folder
    copyfile('%sbirthdeath.png' % fname, '%sbirthdeath-gudhi.png' % fname2)



#Main function
def compute_main(args=None):
    """
    Function to compute persistent topological loops of a dataset using Gudhi + PersLoop".
    Usage:

    ./compute.py DATANAME --max_edge=X --min_pers=X --sparse=X --pre_sparse=X
    
    DATANAME must correspond to a dataset that has been defined in loader.py

    Optional arguments:

        * max_edge   : float (default=None)
        * min_pers   : float (default=None)
        * sparse     : float (default=0.7)
        * pre_sparse : float (default=None)
    """
    

    if args is None:
        pass

    parser = argparse.ArgumentParser(description="Computation of persistent loops using Gudhi + PersLoop",\
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("name", metavar="DATANAME", type=str, help="Name of your dataset (as defined in loader.py)")
    parser.add_argument("--num_pcs", metavar="X", type=int, default=3, help="Dimensionality of dataset, counted in PC space")
    parser.add_argument("--max_edge", metavar="X", type=float, default=None, help="The max edgelength used in the filtration")
    parser.add_argument("--min_pers", metavar="X", type=float, default=None, help="The minimal lifespan for which a cocycle is retained")
    parser.add_argument("--sparse", metavar="X", type=float, default=0.7, help="A Gudhi way to sparsify the computed filtration")
    parser.add_argument("--pre_sparse", metavar="X", type=float, default=0.005, help="A way to sparsify data before computing filtrations")
    parser.add_argument("--method", metavar="FILTER METHOD", type=str, default='binmeans', help="How filtering of data has been achieved.")
    parser.add_argument("--id", metavar="IDENTIFIER", type=str, default=None, help="Identifier to mark any special choices made for filtering")

    args = parser.parse_args()

    compute_homology(args.name, args.max_edge, args.min_pers, args.sparse, args.pre_sparse,\
                  method=args.method, identifier=args.id, num_pcs=args.num_pcs, num_bins=120, num_to_keep=20000) 




#Running as main:
if __name__ == '__main__':
    tkback = True
    compute_main()
else:
    tkback = False
    


