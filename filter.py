import sys
import os
import argparse
import numpy as np
from eofs.standard import Eof
from scipy import stats
from loader import load_raw
from misc import make_path,headfolder



#Filter phase space using a Gaussian kernel density estimator
def filt_raw_kde(name, perc_to_keep, identifier=None,\
                 num_bins=120, max_num=20000, use_eofs=False,\
                 normalise=True, reverse_order=False, overwrite=False):


    #Load the raw data: note that load_raw will return data in
    #the correct shape
    X = load_raw(name)
    shape_before = X.shape


    print(100*'=')
    print("(filt_raw) Dataset = %s" % name)
    print("(filt_raw) Current shape = %s" % (shape_before,))
    print("(filt_raw) Keeping densest %s percent of data" % perc_to_keep)

    #Determine filename and check if the filtered data already exists.
    #If it does, return, unless overwrite=True
    outdir = headfolder+'./Data/Filtered/%s' % name
    make_path(outdir)

    method = 'GaussianKDE'
    if use_eofs:
        method += '_EOFs'

    if identifier is None:
        if reverse_order:
            outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_%s.txt' % (outdir, name, perc_to_keep, num_bins, method)
        else:
            outname = '%s/%s_filtered_densest_%s_percent_%sbins_%s.txt' % (outdir, name, perc_to_keep, num_bins, method)
    else:
        if reverse_order:
            outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_%s_%s.txt' % (outdir, name, perc_to_keep, num_bins, method, identifier)
        else:
            outname = '%s/%s_filtered_densest_%s_percent_%sbins_%s_%s.txt' % (outdir, name, perc_to_keep, num_bins, method, identifier)
    
    
    if os.path.isfile(outname):
        if overwrite:
            print("(filt_raw) Filtered data already exists. Overwriting...")
        else:
            print("(filt_raw) Filtered data already exists. Returning!")
            return
    else:
        pass
    
   
    #Switch to an EOF basis?
    if use_eofs:
        print("(filt_raw) Converting data to EOF space...")
        from eofs.standard import Eof
        if name == 'Lorenz96':
            num_pcs = 4
            EOFs=Eof(X.T)
            X=EOFs.pcs(npcs=num_pcs).T
        else:
            num_pcs = 3
            EOFs=Eof(X.T)
            X=EOFs.pcs(npcs=num_pcs).T
        print("(filt_raw) PCs computed. Proceeding")


    #Normalise dimensions
    if normalise:
        print("(filt_raw) Normalising dimensions.")
        if name == 'Lorenz96':
            X0 = X[0,:]
            X1 = X[1,:]
            X2 = X[2,:]
            X3 = X[3,:]
            X0 /= X0.std()
            X1 /= X1.std()
            X2 /= X2.std()
            X3 /= X3.std()
            X = np.vstack([X0,X1,X2,X3])
        else:
            X0 = X[0,:]
            X1 = X[1,:]
            X2 = X[2,:]
            X0 /= X0.std()
            X1 /= X1.std()
            X2 /= X2.std()
            X = np.vstack([X0,X1,X2])


    #If you've asked for all the data then that's easy
    if (perc_to_keep == 100) or (perc_to_keep is None):
        Xfiltered = X

    else:

        L=X.shape[1]

        #Kernel density estimate
        kde = stats.gaussian_kde(X)
        density = kde(X)

        #Take the perc_to_keep percentile of the density
        percentile = np.percentile(density, 100-perc_to_keep)

        if reverse_order:
            percentile = np.percentile(density, perc_to_keep)
            print("(filt_raw) WARNING: reverse_order = True")
            inds = np.argwhere(density < percentile)[:,0] #<----This gives you the (100-perc_to_keep) LEAST dense points, i.e. outliers etc.
        else:
            inds = np.argwhere(density > percentile)[:,0] #<----This gives you the (100-perc_to_keep) MOST dense points

        Xfiltered = X[:,inds]


    #Saving output
    shape_after = Xfiltered.shape
    np.savetxt(outname, Xfiltered)
    print("(filt_raw) New data shape = %s" % (shape_after,))
    print("(filt_raw) Output saved to: %s" % outname)
    print(100*'=')





#Filter by direct binning in EOF space - after a script written
#by J. Dorrington (2020)
def filt_raw_eofbins(name, perc_to_keep, num_bins=120, num_hist_dims=3,\
                     num_pcs=3,max_num=20000, batch_max=2000, subsample_rate=2,\
                     mode='bin_means', identifier=None,\
                     use_eofs=False, reverse_order=False,\
                     useall=True, normalise=True, overwrite=False):

    #Load the raw data: note that load_raw will return data in
    #the correct shape
    X = load_raw(name, useall)
    shape_before = X.shape

    print(100*'=')
    print("(filt_raw) Dataset = %s" % name)
    print("(filt_raw) Current shape = %s" % (shape_before,))
    if reverse_order:
        print("(filt_raw) Keeping %s percent least dense datapoints" % perc_to_keep)
    else:
        print("(filt_raw) Keeping densest %s percent of data" % perc_to_keep)


    #Determine filename and check if the filtered data already exists.
    #If it does, return, unless overwrite=True
    outdir = headfolder+'./Data/Filtered/%s' % name
    make_path(outdir)
   
    if mode == 'bin_means':
        tag = 'binmeans'
    else:
        tag = 'bincentres'

    if identifier is None:
        if use_eofs:
            if reverse_order:
                outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_EOF%s_EOFs.txt' % (outdir, name, perc_to_keep, num_bins, tag)
            else:
                outname = '%s/%s_filtered_densest_%s_percent_%sbins_EOF%s_EOFs.txt' % (outdir, name, perc_to_keep, num_bins, tag)
        else:
            if reverse_order: 
                outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_EOF%s.txt' % (outdir, name, perc_to_keep, num_bins, tag)
            else:
                outname = '%s/%s_filtered_densest_%s_percent_%sbins_EOF%s.txt' % (outdir, name, perc_to_keep, num_bins, tag)
    else:
        if use_eofs:
            if reverse_order:
                outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_EOF%s_EOFs_%s.txt' % (outdir, name, perc_to_keep, num_bins, tag, identifier)
            else:
                outname = '%s/%s_filtered_densest_%s_percent_%sbins_EOF%s_EOFs_%s.txt' % (outdir, name, perc_to_keep, num_bins, tag, identifier)
        else:
            if reverse_order:
                outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_EOF%s_%s.txt' % (outdir, name, perc_to_keep, num_bins, tag, identifier)
            else:
                outname = '%s/%s_filtered_densest_%s_percent_%sbins_EOF%s_%s.txt' % (outdir, name, perc_to_keep, num_bins, tag, identifier)

    if os.path.isfile(outname):
        if overwrite:
            print("(filt_raw) Filtered data already exists. Overwriting...")
        else:
            print("(filt_raw) Filtered data already exists. Returning!")
            return
    else:
        pass
   

    #Swap to an EOF basis?
    if use_eofs:
        #Compute EOFs of the data
        EOFs=Eof(X.T)
    
        #Get PCs in different dimensions:
        Y=EOFs.pcs(npcs=num_pcs).T
        Yhist=EOFs.pcs(npcs=num_hist_dims).T
    else:
        Y = X
        Yhist = X


    #Normalise dimensions
    if normalise:
        Y0 = Y[0,:]
        Y1 = Y[1,:]
        Y2 = Y[2,:]
        Y0 /= Y0.std()
        Y1 /= Y1.std()
        Y2 /= Y2.std()
        Y = np.vstack([Y0,Y1,Y2])
    
    
    if (perc_to_keep == 100) or perc_to_keep is None:
        Yfiltered = Y
    else:

        L=Y.shape[1]
    
        #From the perc_to_keep, we work out the number of points to keep
        #based on the length of the time series
        percnum=int(perc_to_keep*1e-2*L)
  
        num2keep = percnum

        #We take a histogram in num_hist_dims dimensions of the EOF data
        #Requires a lot of memory fyi
        density,bins=np.histogramdd(Yhist.T,bins=num_bins, normed=True)
        bins=np.array(bins)
        density = np.array(density)
        density_shape = density.shape


        #find the densest bins containing up to perc_to_keep % of the data
        #and store them in 'inds'

        #We find out which bin each data point is in
        which_bin=np.zeros_like(Y)
        for d in range(num_hist_dims):
            which_bin[d]=np.digitize(Yhist[d],bins[d],right=True)-1
        which_bin[which_bin==-1]=0
    
        #We get all the nonempty bins and we sort them in decreasing order of density
        nonempty_bins,non_empty_counts=np.unique(which_bin,axis=1,return_counts=True)
        nonempty_bins=nonempty_bins[:,np.argsort(non_empty_counts)[::-1]]
        non_empty_counts=np.sort(non_empty_counts)[::-1]

        #inds is a set of bins that contain the num_to_keep data points most
        #densely grouped in num_hist_dims dimensions
        inds=nonempty_bins[:,np.cumsum(non_empty_counts)<num2keep].astype(int)
    

        #Take binmeans if asked for
        if mode == "bin_means":

            print("(filt_raw) Method = bin means (%s bins)" % num_bins)
            tag = 'binmeans'
        
            batch_num=L//batch_max +int(bool(L%batch_max))

            tot_sum=np.zeros([num_pcs,inds.shape[1]])
            tot_counts=np.zeros([inds.shape[1]])
        
            #Split the computation up into smaller batches

            for i in range(batch_num):

                #Progress bar
                j = (i + 1) / float(batch_num)
                sys.stdout.write('\r')
                sys.stdout.write("[%-20s] %d%%" % ('='*int(20*j), 100*j))
                sys.stdout.flush()
            
                #The actual computation
            
                #Find which bins the points are in
                batch_bin=which_bin[:,i*batch_max:(i+1)*batch_max]
            
                #Get the corresponding data values in num_pcs dimensions
                #the corresponding data values
                batch_Y=Y[:,i*batch_max:(i+1)*batch_max]
            
                #Find which points are in one of the densest bins stored in inds,
                #and sum them
                bin_to_data=np.all(batch_bin[:,:,None]==inds[:,None,:],axis=0)
                batch_sum=np.sum(batch_Y[:,:,None]*bin_to_data[None,:,:],axis=1)
                tot_sum+=batch_sum
                tot_counts+=bin_to_data.sum(axis=0)
        
            print('') # blank line due to progress bar
            print("(filt_raw) Filtering complete!")
        
            #Turn the summed data into averaged data in each bin
            Yfiltered=tot_sum/tot_counts
        
            Yfiltered=Yfiltered[:,::subsample_rate][:,:max_num]


        #Otherwise just take the centres of these bins
        elif mode == "bin_centres":
            if num_hist_dims!=num_pcs:
                raise(ValueError("bin_centres makes no sense for num_hist_dims!=num_pcs"))
            
            print("(filt_raw) Method = bin centres (%s bins)" % num_bins)

            tag = 'bincentres'
            bin_widths=np.array((bins[:,-1]-bins[:,-2])/2)
            Yfiltered = np.array([bins[i][inds[i,:]] + bin_widths[i] for i in range(num_pcs)])

            Yfiltered=Yfiltered[:,::subsample_rate][:,:max_num]
            print("(filt_raw) Filtering complete!")


        else:
            raise(IOError(f"Unrecognised filtering mode {mode}"))
    

    #Saving output
    shape_after = Yfiltered.shape
    np.savetxt(outname, Yfiltered)
    print("(filt_raw) New data shape = %s" % (shape_after,))
    print("(filt_raw) Output saved to: %s" % outname)
    print(100*'=')




#Filter by phase space velocity
def filt_raw_speed(name, perc_to_keep, identifier=None,\
                   num_bins=120, max_num=20000, use_eofs=False,\
                   normalise=True, reverse_order=False, overwrite=False):


    #Load the raw data: note that load_raw will return data in
    #the correct shape
    X = load_raw(name)
    shape_before = X.shape


    print(100*'=')
    print("(filt_raw) Dataset = %s" % name)
    print("(filt_raw) Current shape = %s" % (shape_before,))
    print("(filt_raw) Keeping slowest %s percent of data" % perc_to_keep)

    #Determine filename and check if the filtered data already exists.
    #If it does, return, unless overwrite=True
    outdir = headfolder+'./Data/Filtered/%s' % name
    make_path(outdir)

    method = 'PhaseSpeed'
    if use_eofs:
        method += '_EOFs'

    if identifier is None:
        if reverse_order:
            outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_%s.txt' % (outdir, name, perc_to_keep, num_bins, method)
        else:
            outname = '%s/%s_filtered_densest_%s_percent_%sbins_%s.txt' % (outdir, name, perc_to_keep, num_bins, method)
    else:
        if reverse_order:
            outname = '%s/%s_filtered_leastdense_%s_percent_%sbins_%s_%s.txt' % (outdir, name, perc_to_keep, num_bins, method, identifier)
        else:
            outname = '%s/%s_filtered_densest_%s_percent_%sbins_%s_%s.txt' % (outdir, name, perc_to_keep, num_bins, method, identifier)
    
    
    if os.path.isfile(outname):
        if overwrite:
            print("(filt_raw) Filtered data already exists. Overwriting...")
        else:
            print("(filt_raw) Filtered data already exists. Returning!")
            return
    else:
        pass
    
   
    #Switch to an EOF basis?
    if use_eofs:
        print("(filt_raw) Converting data to EOF space...")
        from eofs.standard import Eof
        num_pcs = 3
        EOFs=Eof(X.T)
        X=EOFs.pcs(npcs=num_pcs).T
        print("(filt_raw) PCs computed. Proceeding")


    #Normalise dimensions
    if normalise:
        X0 = X[0,:]
        X1 = X[1,:]
        X2 = X[2,:]
        X0 /= X0.std()
        X1 /= X1.std()
        X2 /= X2.std()
        X = np.vstack([X0,X1,X2])


    #If you've asked for all the data then that's easy
    if (perc_to_keep == 100) or (perc_to_keep is None):
        Xfiltered = X

    else:

        L=X.shape[1]

        #Compute phase speeds: because the timestep between each point is always 1 day,
        #it is equivalent to just computing the Euclidean distance between successive points.
        speeds = np.sqrt(np.sum((X[:,1:]-X[:,:-1])**2,0)) 
        speeds = np.insert(speeds,0,speeds[0])



        #Take the perc_to_keep percentile of the density
        percentile = np.percentile(speeds, 100-perc_to_keep)

        if reverse_order:
            percentile = np.percentile(speeds, perc_to_keep)
            print("(filt_raw) WARNING: reverse_order = True")
            inds = np.argwhere(speeds < percentile)[:,0] #<----This gives you the (100-perc_to_keep) LEAST dense points, i.e. outliers etc.
        else:
            inds = np.argwhere(speeds > percentile)[:,0] #<----This gives you the (100-perc_to_keep) MOST dense points

        Xfiltered = X[:,inds]


    #Saving output
    shape_after = Xfiltered.shape
    np.savetxt(outname, Xfiltered)
    print("(filt_raw) New data shape = %s" % (shape_after,))
    print("(filt_raw) Output saved to: %s" % outname)
    print(100*'=')

