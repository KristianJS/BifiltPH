
'''
The main script to compute everything.
'''


#Bifiltration main
def bifiltration(name, percentages, method,\
                 eofs=False,\
                 reverse=False,\
                 overwrite=False,\
                 nohist=False,\
                 use_pl=True,\
                 plots=False):

    import matplotlib
    matplotlib.use('Agg')
    import sys
    import os
    import numpy as np
    import plots as p 
    from scipy import stats
    from loader import load_raw, load_filt
    from misc import make_path, headfolder, default_maxedge, default_minpers, fmt, dpi, fsize
    from filter import filt_raw_kde, filt_raw_speed, filt_raw_eofbins
    from compute import compute_homology


    #Only support these 3 options for now
    assert method in ['GaussianKDE', 'EOFbinmeans', 'PhaseSpeed']
    

    #Fix some basics
    num_bins = 160


    #Print some preamble
    print(100*'=')
    print("                            GENERATING BIFILTRATION")
    print(100*'-')
    print("  Dataset: %s" % name)
    print("  Density percentages: %s" % percentages)
    print("  Filtering method: %s" % method)
    if reverse:
        print("  Filtering order: reversed")
    else:
        print("  Filtering order: standard")
    print(100*'-')


    #Parameter choices determined through trial/error
    if name == 'Lorenz63':
        normalise = True
        if eofs:
            if normalise:
                max_edge = 5.0
                min_pers = 0.25
                pre_sparse = 0.05
                sparse = 0.7
            else:
                max_edge = 20.0
                min_pers = 0.6
                pre_sparse = 0.3
                sparse = 0.7
        else:
            if normalise:
                max_edge = 5.0
                min_pers = 0.25
                pre_sparse = 0.05
                sparse = 0.7
            else:
                max_edge = 20.0
                min_pers = 0.6
                sparse = 0.7
                pre_sparse = 0.5
        
        if method == 'EOFbinmeans':
            pre_sparse = 0.05
        plot_loops = True
        plot_comps = True
        num2show = 2
        num2keep = 5
        identifier = None
    
    elif name == 'Lorenz96':
        if eofs:
            max_edge = 5.0
            min_pers = 0.25
            pre_sparse = 0.05
            sparse = 0.7
        else:
            max_edge = 5.0
            min_pers = 1.4
            pre_sparse = 0.05
            sparse = 0.7
        
        if method == 'EOFbinmeans':
            pre_sparse = 0.05
        plot_loops = True
        plot_comps = True
        num2show = 3
        num2keep = 5
        normalise = True
        identifier = None
    

    elif 'CDV' in name:
        max_edge = 5.0
        min_pers = 0.25
        min_pers = 0.5
        if method == 'PhaseSpeed':
            if reverse:
                pre_sparse = 0.015
            else:
                pre_sparse = 0.05
        else:
            pre_sparse = 0.005
        sparse = 0.7
        plot_loops = True
        plot_comps = True
        num2show = 3
        num2keep = 5
        eofs = True
        normalise = True
        identifier = None

    elif name == 'JetLat':
        max_edge = 3.50
        min_pers = 0.15
        pre_sparse = None
        sparse = 0.7
        plot_loops = True
        plot_comps = True
        num2show = 3
        num2keep = 5
        eofs = False #we're already in EOF space
        normalise = True #technically redundant
        identifier = None
    
    elif name == 'Gaussian':
        max_edge = 5.0
        min_pers = 0.15 
        pre_sparse = 0.05
        sparse = 0.7
        plot_loops = True
        plot_comps = True
        num2show = 2
        num2keep = 5
        eofs = False
        normalise = False
        identifier = None

    else:
        sys.exit("(bifiltration) ERROR: you need to specify parameter choices up front for new datasets!")


    kwargs = {'max_edge':max_edge, 'min_pers':min_pers, 'sparse':sparse, 'pre_sparse':pre_sparse,\
            'num_bins':num_bins, 'reverse_order':reverse}

    for perc in percentages:
        print(100*'-')
        print("                              PERCENTAGE = %s" % perc)
 
        if not(plots):
            if method == 'GaussianKDE':
                filt_raw_kde(name, perc_to_keep=perc, num_bins=num_bins, reverse_order=reverse,\
                             normalise=normalise, use_eofs=eofs, identifier=identifier, overwrite=overwrite)
            elif method == 'EOFbinmeans':
                filt_raw_eofbins(name, perc_to_keep=perc, num_bins=num_bins, num_hist_dims=3,\
                                 max_num=20000, batch_max=2000, subsample_rate=2,\
                                 mode='bin_means', identifier=identifier,\
                                 use_eofs=eofs, reverse_order=reverse, overwrite=overwrite)
            elif method == 'PhaseSpeed':
                filt_raw_speed(name, perc_to_keep=perc, num_bins=num_bins, reverse_order=reverse,\
                             normalise=normalise, use_eofs=eofs, identifier=identifier, overwrite=overwrite)

            else:
                sys.exit("ERROR: unrecognized method %s" % method)

            compute_homology('%s-Perc%s' % (name, perc), perc_to_keep=perc, method=method, eofs=eofs,\
                             num_comp_to_keep=num2keep, use_pl=use_pl, **kwargs)
       
        if plot_loops and use_pl:
            p.phasespace_loops('%s-Perc%s' % (name, perc), dims=[0,1,2], num2show=num2show, method=method, eofs=eofs,\
                               save=True, **kwargs)
        

        if plot_comps:
            p.phasespace_comps('%s-Perc%s' % (name, perc), dims=[0,1,2], num2keep=num2keep, num2show=num2show, method=method, eofs=eofs,\
                               save=True, **kwargs)
   
    if len(percentages)<4:
        pass
    else:
        p.density_vs_lifetime(name, percentages, num2show=num2keep, identifier=None, method=method, eofs=eofs, save=True, **kwargs)
   
    if not(nohist):
        p.density_histogram(name, percentages, identifier=None, method=method, save=True, **kwargs)
    
    
    print(100*'=')
    print("                            BIFILTRATION COMPLETE")
    print(100*'-')



##########################################################################
#                          RUNNING AS MAIN
##########################################################################
#Main function
def main(args=None):
    """
    Function to perform a full bifiltration across density/length space, using
    Gudhi and Persloop to compute the homology and find long-lived cycles.
    Produces standard output of plots for a given dataset.

    Usage:

    ./bifiltration.py DATANAME --percs=<VALS> --method=[Kernel,Bins] [--eofs, --binmeans, --reverse, --overwrite, --nohist]
    
    DATANAME must correspond to a dataset that has been defined in loader.py
    It is recommended to specify default Gudhi/Persloop parameters for your dataset
    (max_edge, min_pers, sparse and pre_sparse) in the bifiltration function to
    speed up computations.

    Optional arguments:

        Positional:
          --percs      : Percentages to filter density over: 'default' or a list of ints (e.g. [10,50,90])
          --method     : Filtering method. 'Kernel' (Gaussian KDE; default) or 'Binmeans' (direct binning).
        
        Boolean:
          --eofs      : Transform data to EOF space before filtering/computing (default=False)
          --reverse   : Filter by decreasing density as opposed to default of increasing density.
          --overwrite : Force overwrite of any existing filtered data. Otherwise existing files are used.
          --nohist    : Don't produce a density histogram (useful because this can be a slow operation).
    """

    #Module import
    import sys
    import os
    import argparse


    if args is None:
        pass

    parser = argparse.ArgumentParser(description="Full bifiltration of homology using Gudhi + PersLoop",\
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("name", metavar="DATANAME", type=str, help="Name of your dataset (as defined in loader.py)")
    parser.add_argument("--percs", metavar="Filtvals", default='default', choices=['default', 'minimal', 'low', 'fine', 'manual'], help="Choice of density filtration percentages.")
    parser.add_argument("--method", type=str, default='kernel', choices=['kernel', 'bins', 'speed'], help="How to do the density filtering.")
    parser.add_argument("--eofs", action='store_true', help="Change of basis to EOF space.")
    parser.add_argument("--reverse", action='store_true', help="Reverse order of density filtration.")
    parser.add_argument("--overwrite", action='store_true', help="Force overwrite any existing filtered data.")
    parser.add_argument("--nohist", action='store_true', help="Prevent creation of density histogram plot.")
    parser.add_argument("--nopl", action='store_true', help="Don't use PersLoop, only do stuff which uses Gudhi.")
    parser.add_argument("--plots", action='store_true', help="Only create plots, don't compute anything.")

    args = parser.parse_args()

    if args.percs == 'default':
        percentages = [10,20,30,40,50,60,70,80,90,100]
    elif args.percs == 'fine':
        percentages = list(range(5,100,5))
    elif args.percs == 'minimal':
        percentages = [10,80]
    elif args.percs == 'low':
        percentages = [10,15,20,25]
    elif args.percs == 'manual':
        percentages = [30]
    else:
        print("ERROR: --percs must be 'default' or 'minimal'")

    if args.method == 'kernel':
        method = 'GaussianKDE'
    elif args.method == 'bins':
        method = 'EOFbinmeans'
    else:
        method = 'PhaseSpeed'

    bifiltration(args.name, percentages, method, args.eofs, args.reverse, args.overwrite, args.nohist, not(args.nopl), args.plots) 


if __name__ == '__main__':
   main() 


