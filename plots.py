import sys
import os
import argparse
import glob
import shutil
import numpy as np
import numpy.ma as ma
from scipy import stats
from eofs.standard import Eof
import matplotlib
if __name__ == '__main__':
    matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from misc import make_path, default_maxedge, default_minpers, headfolder, fmt, dpi, fsize
import loader


'''
Script for creating generic suite of plots.
'''



#Where figures end up
figpath = headfolder + '/Figures'




##########################################################################################


#Add a generic 3D plot to an existing ax instance(to be called within other functions)
def plot_3d_generic(x, y, z,\
                    ax,\
                    label=None,\
                    title=None,\
                    color='k',\
                    marker='o',\
                    cs=None,\
                    s=0.1,\
                    alpha=1.,\
                    fontsize=12.,\
                    labelpad=-10,\
                    tickpad=-5,\
                    eofs=False,\
                    axnames=['x', 'y', 'z'],\
                    xlims=(None,None),\
                    ylims=(None,None),\
                    zlims=(None,None),\
                    elev=10, azim=-40):

   
    if cs is None:
        p=ax.scatter(x, y, z, s=s, color=color, alpha=alpha, marker=marker, label=label)
    else:
        p=ax.scatter(x, y, z, s=s, c=cs, cmap='RdBu_r', alpha=alpha, marker=marker, label=label)

    if not(xlims[0] is None and xlims[1] is None):
        ax.set_xlim(xmin=xlims[0],xmax=xlims[1])
    if not(ylims[0] is None and ylims[1] is None):
        ax.set_ylim(ymin=ylims[0],ymax=ylims[1])
    if not(zlims[0] is None and zlims[1] is None):
        ax.set_zlim(zmin=zlims[0],zmax=zlims[1])
   
    if eofs:
        axnames=['PC1', 'PC2', 'PC3'] 

    
    for t in ax.xaxis.get_major_ticks(): 
        t.label.set_fontsize(6)
        t.set_pad(tickpad)
    for t in ax.yaxis.get_major_ticks():
        t.label.set_fontsize(6)
        t.set_pad(tickpad)
    for t in ax.zaxis.get_major_ticks():
        t.label.set_fontsize(6)
        t.set_pad(tickpad)

    if not(axnames is None):
        ax.set_xlabel(axnames[0], fontsize=7.5, labelpad=labelpad)
        ax.set_ylabel(axnames[1], fontsize=7.5, labelpad=labelpad)
        ax.set_zlabel(axnames[2], fontsize=7.5, labelpad=labelpad)
    
    ax.view_init(elev=elev, azim=azim)

    if not(title is None):
        plt.title(title, fontsize=fontsize)
    
    if not(label is None):
        plt.legend(loc='best')

    return p



#Plot of a dataset using increasingly many points
def phasespace_points(name, eofs=False, dims=[0,1,2], cb=True, save=True):

    #Points to plot
    Ns = [10000, 20000, 50000]

    print("(phasespace_points) Generating raw phase space plot...")
    print("(phasespace_points) Dataset = %s" % name)
    print("(phasespace_points) Number of points to be used: %s" % Ns)


    #Load the raw data
    data = loader.load_raw(name)

    #Optionally transform to EOFs
    if eofs:
        EOFs=Eof(data.T)
        pcdata=EOFs.pcs(npcs=3).T
        data = pcdata

    #Normalise
    X0 = data[0,:]
    X1 = data[1,:]
    X2 = data[2:,:]
    X0 /= X0.std()
    X1 /= X1.std()
    X2 /= X2.std()
    data = np.vstack([X0,X1,X2])
    
    #Grab first N timesteps for each N
    data0 = data[:,:Ns[0]]
    data1 = data[:,:Ns[1]]
    data2 = data[:,:Ns[2]]

    #Compute density for the biggest one
    print("(phasespace_raw) Computing density...")
    kde = stats.gaussian_kde(data2)
    density = kde(data2)


    #Some sensible choices for our datasets
    if ('Lorenz63' in name):
        ss = 0.1
        alpha = 0.01
        if eofs:
            elev = 10
            azim = -80
            xlims = (-3,3)
            ylims = (-3,3)
            zlims = (-1,1)
        else:
            elev = 10
            azim = -40
            xlims = (-2,2)
            ylims = (-2,2)
            zlims = (0,4)
            
    elif 'Lorenz96' in name:
        ss = 0.1
        alpha = 0.01
        if eofs:
            elev = 50
            azim = 20
            xlims = (-2,2)
            ylims = (-2,2)
            zlims = (-2,2)
        else:
            elev = 50
            azim = 20
            xlims = (-2,2)
            ylims = (-2,2)
            zlims = (-2,2)

    elif 'CDV' in name:     
        ss = 0.1
        alpha = 0.01
        elev = 10
        azim = -40
        xlims = (-3,3)
        ylims = (-2,4)
        zlims = (-4,4)
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)

    elif 'JetLat' in name:
        ss = 0.1
        alpha = 0.1
        elev = 10
        azim = -40
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)

    else:
        ss = 0.1
        alpha = 0.1
        elev = 10
        azim = -40
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)


    if name == 'Lorenz':
        tname = 'Lorenz63'
    else:
        tname = name

    fontsize=10
    s=0.2

    #Begin plotting
    print("(phasespace_raw) Plotting...")
    fig = plt.figure(figsize=fsize)
    
    if eofs:
        axnames=['PC1', 'PC2', 'PC3']
    else:
        axnames=['x', 'y', 'z']
    
    d0 = dims[0]
    d1 = dims[1]
    d2 = dims[2]
    
    ax = plt.subplot(221, projection = '3d')
    p=plot_3d_generic(data0[d0,:], data0[d1,:], data0[d2,:],\
                    ax,\
                    label=None,\
                    title='(a) %s: first %s timesteps' % (tname, Ns[0]),\
                    cs=None,\
                    s=s,\
                    fontsize=fontsize,\
                    marker='o',\
                    alpha=0.6,\
                    axnames=axnames,\
                    elev=elev, azim=azim)

    ax = plt.subplot(222, projection = '3d')
    p=plot_3d_generic(data1[d0,:], data1[d1,:], data1[d2,:],\
                    ax,\
                    label=None,\
                    title='(b) %s: first %s timesteps' % (tname, Ns[1]),\
                    cs=None,\
                    s=s,\
                    fontsize=fontsize,\
                    marker='o',\
                    alpha=0.6,\
                    axnames=axnames,\
                    elev=elev, azim=azim)
    
    ax = plt.subplot(223, projection = '3d')
    p=plot_3d_generic(data2[d0,:], data2[d1,:], data2[d2,:],\
                    ax,\
                    label=None,\
                    title='(c) %s: first %s timesteps' % (tname, Ns[2]),\
                    cs=None,\
                    s=s,\
                    fontsize=fontsize,\
                    marker='o',\
                    alpha=0.6,\
                    axnames=axnames,\
                    elev=elev, azim=azim)
    
    ax = plt.subplot(224, projection = '3d')
    
    if cb:
        axnames=[axnames[0],axnames[1],None]
    p=plot_3d_generic(data2[d0,:], data2[d1,:], data2[d2,:],\
                    ax,\
                    label=None,\
                    title='(d) %s: density estimate' % (tname),\
                    cs=density,\
                    s=s,\
                    fontsize=fontsize,\
                    marker='o',\
                    alpha=0.6,\
                    axnames=axnames,\
                    elev=elev, azim=azim)
    
    if cb:
        axins = inset_axes(ax, width="5%",height="60%",loc='lower left',bbox_to_anchor=(1, 0.2, 0.6, 1),bbox_transform=ax.transAxes,borderpad=0)
        fig.colorbar(p, cax=axins, orientation='vertical')

    #Saving/showing
    if save:
        if eofs:
            fname = '%s_rawphasespace_EOFs_density_evolution.%s' % (name, fmt)
        else:
            fname = '%s_rawphasespace_density_evolution.%s' % (name, fmt)
        fpath = '%s' % (figpath)
        figname = '%s/%s' % (fpath, fname)
        make_path(fpath)
        plt.savefig(figname, format=fmt, dpi=dpi)
        plt.close()
        print("(phasespace_raw) Figure saved as %s" % figname)
        print(100*'-')
    else:
        plt.show()



#Plot of filtered phase space with and without loops overlaid on top
def phasespace_loops(name, dims=[0,1,2], num2show=2, max_edge=None, min_pers=None, sparse=0.7, pre_sparse=None,\
                          method='GaussianKDE', num_bins=120, eofs=False,\
                          identifier=None, reverse_order=False,\
                          save=True):
    
    print("(phasespace_loops) Generating phase space plot with loops...")
   
    #Get basename
    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    
    
    #Get filtered data and normalise axes
    filtdata = loader.load_from_pers(name, num_bins=num_bins, max_edge=max_edge, min_pers=min_pers,\
                                     sparse=sparse, pre_sparse=pre_sparse, method=method,\
                                     reverse_order=reverse_order, eofs=eofs)
    

    normdata = filtdata

    #Get the loops
    loop_list= loader.load_loops(name, max_edge=max_edge, min_pers=min_pers,\
                                     sparse=sparse, pre_sparse=pre_sparse,\
                                     method=method, reverse_order=reverse_order, eofs=eofs)
    


    #If max_edge/min_pers = None, determine what these actually are
    if (max_edge is None):
        max_edge = default_maxedge(loader.load_filt(name, method, identifier, num_bins, None, eofs=eofs, reverse_order=reverse_order))

    if min_pers is None:
        min_pers = default_minpers(loader.load_filt(name, method, identifier, num_bins, None, eofs=eofs, reverse_order=reverse_order), name)

    max_edge_str = "{0:.2f}".format(max_edge) #strip away decimals to avoid overly long filenames
    min_pers_str = "{0:.2f}".format(min_pers)
    sparse_str = "{0:.2f}".format(sparse)

    if pre_sparse is None:
        presparse_str = '0'
    else:
        presparse_str = "{0:.3f}".format(pre_sparse)    


    #Some sensible choices for our datasets
    if ('Lorenz63' in name):
        ss = 0.1
        alpha = 0.01
        if eofs:
            elev = 10
            azim = -80
            xlims = (-3,3)
            ylims = (-3,3)
            zlims = (-3,3)
        else:
            elev = 10
            azim = -40
            xlims = (-2,2)
            ylims = (-2,2)
            zlims = (0,4)
            
    elif 'Lorenz96' in name:
        ss = 0.1
        alpha = 0.01
        if eofs:
            elev = 50
            azim = 20
            xlims = (-2,2)
            ylims = (-2,2)
            zlims = (-2,2)
        else:
            elev = 50
            azim = 20
            xlims = (-2,2)
            ylims = (-2,2)
            zlims = (-2,2)

    elif 'CDV' in name:     
        ss = 0.1
        alpha = 0.01
        elev = 10
        azim = -40
        xlims = (-3,3)
        ylims = (-2,4)
        zlims = (-4,4)
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)

    elif 'JetLat' in name:
        ss = 0.1
        alpha = 0.1
        elev = 10
        azim = -40
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)

    else:
        ss = 0.1
        alpha = 0.1
        elev = 10
        azim = -40
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)



    #Begin plotting
    fig = plt.figure(figsize=fsize)
    axnames=['x', 'y', 'z']
    
    d0 = dims[0]
    d1 = dims[1]
    d2 = dims[2]
    
    if reverse_order:
        tit_a = '(a) %s: %s (rev.)' % (name, method)
    else:
        tit_a = '(a) %s: %s' % (name, method)
    tit_b = '(b) With persistent loops'


    #Filtered phase space by itself
    ax1 = plt.subplot(121, projection = '3d')
    plot_3d_generic(normdata[d0,:], normdata[d1,:], normdata[d2,:],\
                    ax1,\
                    label=None,\
                    title=tit_a,\
                    color='k',\
                    marker='o',\
                    s=0.5,\
                    xlims=xlims,\
                    ylims=ylims,\
                    zlims=zlims,\
                    alpha=0.8,\
                    axnames=axnames,\
                    elev=elev, azim=azim)
    
    
    
    #Then with loops overlaid
    ax2 = plt.subplot(122, projection='3d') 
    plot_3d_generic(normdata[d0,:], normdata[d1,:], normdata[d2,:],\
                    ax2,\
                    label=None,\
                    title=tit_b,\
                    color='k',\
                    marker='o',\
                    s=ss,\
                    xlims=xlims,\
                    ylims=ylims,\
                    zlims=zlims,\
                    alpha=0.3,\
                    axnames=axnames,\
                    elev=elev, azim=azim)


    
    loop_colors = ['r', 'b', 'g', 'y', 'm', 'tab:brown', 'tab:pink', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red', 'tab:olive', 'tab:cyan',
                   'xkcd:cloudy blue', 'xkcd:dark pastel green', 'xkcd:dust', 'xkcd:electric lime','xkcd:tea']

    for i in range(min(len(loop_list),num2show)):
        col = loop_colors[i]
        l = loop_list[i]
        ax2.plot(l[:,d0], l[:,d1], l[:,d2], color=col, alpha=1., linewidth=3.0)



    #Saving/showing
    plt.tight_layout()
    if save:
        if eofs:
            fname = '%s_%s_EOFs_maxe_%s_minp_%s_sp_%s_presp_%s_phasespace-loops.%s' % \
                   (name, method, max_edge_str, min_pers_str, sparse_str, presparse_str, fmt)
        else:
            fname = '%s_%s_maxe_%s_minp_%s_sp_%s_presp_%s_phasespace-loops.%s' % \
                   (name, method, max_edge_str, min_pers_str, sparse_str, presparse_str, fmt)
        if reverse_order:
            fpath = '%s/%s/reversed/%s' % (figpath, basename, name)
        else:
            fpath = '%s/%s/standard/%s' % (figpath, basename, name)
        figname = '%s/%s' % (fpath, fname)
        make_path(fpath)
        plt.savefig(figname, format=fmt, dpi=dpi)
        plt.close()
        print("(phasespace_loops) Figure saved as %s" % figname)
        print(100*'-')
    else:
        plt.show()




#Bespoke plot for each dataset
def phasespace_bespoke(name, save=True):
    
    print("(phasespace_bespoke) Generating phase space plot with long-lived cycles...")

    if name == 'Lorenz63':
        tname = name
        max_edge = 5.0
        min_pers = 0.25
        pre_sparse = 0.05
        sparse = 0.7
        percs = [80,60,20]
        cycles = ['loops', 'loops', 'comps']
        ss = 0.1
        alpha = 0.05
        elev = 10
        azim = -40
        xlims = (-2.5,2.5)
        ylims = (-2.5,2.5)
        zlims = (0,4.5)
        eofs = False
        num2show = 1
        num_bins = 160
        method = 'GaussianKDE'

    elif name == 'Lorenz96':
        tname = name
        max_edge = 5.0
        min_pers = 0.25
        pre_sparse = 0.05
        sparse = 0.7
        percs = [50,20,10]
        cycles = ['comps', 'loops', 'comps']
        ss = 0.1
        alpha = 0.05
        elev = 50
        azim = 20
        xlims = (-2,2)
        ylims = (-2,2)
        zlims = (-2,2)
        eofs = True
        num2show = 3
        num_bins = 160
        method = 'GaussianKDE'

    elif 'CDV' in name:     
        tname = name
        max_edge = 5.0
        min_pers = 0.50
        pre_sparse = 0.005
        sparse = 0.7
        percs = [70,50,20]
        cycles = ['loops', 'loops', 'comps']
        ss = 0.4
        alpha = 0.05
        elev = 10
        azim = -40
        xlims = (-3,3)
        ylims = (-2,4)
        zlims = (-4,4)
        eofs = True
        num2show = 3
        num_bins = 160
        method = 'EOFbinmeans'

    elif 'JetLat' in name:
        tname = name
        max_edge = 3.50
        min_pers = 0.15
        pre_sparse = None
        sparse = 0.7
        percs = [70,50,10]
        cycles = ['comps', 'comps', 'comps']
        ss = 1.0
        alpha = 0.2
        elev = 10
        azim = -40
        xlims = (4.5,7.5)
        ylims = (-2,2)
        zlims = (-2,2)
        eofs = False
        num2show = 3
        num_bins = 160
        method = 'GaussianKDE'

    elif name == 'Gaussian':
        tname = name
        max_edge = 2.50
        min_pers = 0.15
        pre_sparse = 0.050
        sparse = 0.7
        percs = [70,40,10]
        cycles = ['comps', 'comps', 'comps']
        ss = 1.0
        alpha = 0.2
        elev = 10
        azim = -40
        xlims = (-2,2)
        ylims = (-2,2)
        zlims = (-2,2)
        eofs = False
        num2show = 3
        num_bins = 160
        method = 'GaussianKDE'

    
    else:
        ss = 0.1
        alpha = 0.1
        elev = 10
        azim = -40
        xlims = (None,None)
        ylims = (None,None)
        zlims = (None,None)
    
    #Get basename
    if 'Perc' in name:
        basename, percstr = name.split('-')
        perc_to_keep = percstr.split('Perc')[-1]
    

    filt_dict = {}
    cycle_dict = {}

    filt_dict[100] = loader.load_filt('%s-Perc%s' % (name, 100), method, None,\
                                                  num_bins=160, perc_to_keep=100,\
                                                  reverse_order=False, eofs=eofs, normalise=False)[:,:100000]

    

    for n in range(len(percs)):

        filtdata = loader.load_from_pers('%s-Perc%s' % (name, percs[n]), num_bins=num_bins, max_edge=max_edge, min_pers=min_pers,\
                                         sparse=sparse, pre_sparse=pre_sparse, method=method,\
                                         reverse_order=False, eofs=eofs)
        
        if cycles[n] == 'comps':
            cycle = loader.load_comps('%s-Perc%s' % (name, percs[n]), 5,\
                                      max_edge, min_pers,\
                                      sparse, pre_sparse,\
                                      method, eofs=eofs, reverse_order=False)
        else:
            cycle = loader.load_loops('%s-Perc%s' % (name, percs[n]), max_edge=max_edge, min_pers=min_pers,\
                                     sparse=sparse, pre_sparse=pre_sparse,\
                                     method=method, reverse_order=False, eofs=eofs)

        filt_dict[percs[n]] = filtdata
        cycle_dict[percs[n]] = cycle


    
    #Plotting
    colors = ['r', 'b', 'g']
    tits = ['(b)', '(c)', '(d)', '(e)', '(f)']
    d0 = 0
    d1 = 1
    d2 = 2

    if eofs:
        axnames = ['PC1', 'PC2', 'PC3']
    else:
        axnames=['x', 'y', 'z']
   

    fig = plt.figure(figsize=fsize)
    
    #First the 100% data
    ax = plt.subplot(221, projection = '3d')
    plot_3d_generic(filt_dict[100][d0,:], filt_dict[100][d1,:], filt_dict[100][d2,:],\
                    ax,\
                    label=None,\
                    title='(a) %s: 100%% densest' % tname,\
                    color='k',\
                    marker='o',\
                    s=0.8,\
                    xlims=xlims,\
                    ylims=ylims,\
                    zlims=zlims,\
                    alpha=alpha,\
                    axnames=axnames,\
                    elev=elev, azim=azim)
  

    #Then the specifically chosen thresholds
    for n in range(len(percs)):
        normdata = filt_dict[percs[n]]

        ax = plt.subplot('22%s' % (n+2), projection = '3d')
        if cycles[n] == 'loops':
            if name == 'Lorenz96':
                loop_list = cycle_dict[percs[n]][:3]
            else:
                loop_list = cycle_dict[percs[n]][:2]

            if name == 'CDV':
                lw = 2.5
                alpha1 = 0.4
            else:
                lw = 3.0
                alpha1 = 0.4

            plot_3d_generic(normdata[d0,:], normdata[d1,:], normdata[d2,:],\
                            ax,\
                            label=None,\
                            title='%s %s: %s%% densest' % (tits[n], tname, percs[n]),\
                            color='k',\
                            marker='o',\
                            s=ss,\
                            xlims=xlims,\
                            ylims=ylims,\
                            zlims=zlims,\
                            alpha=alpha1,\
                            axnames=axnames,\
                            elev=elev, azim=azim)

            for i in range(len(loop_list)):
                col = colors[i]
                l = loop_list[i]
                ax.plot(l[:,d0], l[:,d1], l[:,d2], color=col, alpha=1., linewidth=lw)

        elif cycles[n] == 'comps':
            components = cycle_dict[percs[n]]

            for k in range(min(len(components), num2show)):
                points = components[k]
                if name in ['Lorenz96', 'CDV', 'Gaussian', 'JetLat']:
                    alpha2 = 0.4
                else:
                    alpha2 = 0.2
                if len(points) < 4:
                    s = 20.
                    alpha2 = 0.8 
                else:
                    s = 0.5
                if k == 0:
                    titx = '%s %s: %s%% densest' % (tits[n], tname, percs[n])
                else:
                    titx = None

                plot_3d_generic(normdata[d0,points], normdata[d1,points], normdata[d2,points],\
                                ax,\
                                label=None,\
                                title=titx,\
                                color=colors[k],\
                                marker='o',\
                                s=s,\
                                xlims=xlims,\
                                ylims=ylims,\
                                zlims=zlims,\
                                alpha=alpha2,\
                                axnames=axnames,\
                                elev=elev, azim=azim)

    
    
    #Saving/showing
    if save:
        fname = '%s_chosen_thresholds.png' % name
        fpath = '%s' % figpath
        figname = '%s/%s' % (fpath, fname)
        make_path(fpath)
        plt.savefig(figname, format=fmt, dpi=dpi)
        plt.close()
        print("(phasespace_bespoke) Figure saved as %s" % figname)
        print(100*'-')
    else:
        plt.show()




#Plot the bifiltration lifetime plot: density vs lifetime of comps/cycles
def density_vs_lifetime(name, percentages='default', max_edge=None, min_pers=None, sparse=0.7, pre_sparse=None,\
                          method='GaussianKDE', num_bins=120, eofs=False,\
                          num2show=10, identifier=None,\
                          reverse_order=False, save=True):

    print("(density_vs_lifetime) Generating plots of density vs lifetimes...")


    if percentages == 'default':
        percentages = [10,20,30,40,50,60,70,80,90,100]


    #Load the death_times of components and birth/death_times of loops, for each percentage
    comp_lives_dict = {}
    comp_sizes_dict = {}
    loop_lives_dict = {}
    max_life_arr = []
    for perc in percentages:
        pname = '%s-Perc%s' % (name, perc)
        comp_deaths, comp_sizes = loader.load_comp_deaths(pname, num2show, max_edge, min_pers,\
                                                          sparse, pre_sparse, method,\
                                                          reverse_order, eofs=eofs)

        loop_births, loop_deaths = loader.load_loop_lives(pname, num2show, max_edge, min_pers,\
                                                          sparse, pre_sparse, method,\
                                                          reverse_order, eofs=eofs)
        


        #If these are empty (no loops found) then add in the min_pers as a placeholder
        if len(comp_deaths) == 0:
            comp_deaths = np.array(num2show*[min_pers])
            comp_sizes = np.array(num2show*[0])
        if len(loop_births) == 0:
            loop_deaths = np.array(num2show*[min_pers])
            loop_births = np.array(num2show*[0.])
        
        loop_lives = loop_deaths - loop_births
    
        comp_lives_dict[perc] = comp_deaths
        comp_sizes_dict[perc] = comp_sizes
        loop_lives_dict[perc] = loop_lives
        max_life_arr.append(ma.masked_invalid(comp_deaths).max())
        max_life_arr.append(ma.masked_invalid(loop_lives).max())

    #Find the maximum lifespan across all cycles
    max_life = ma.array(max_life_arr).max()
    inf_val = max_life*1.2 #Place infinity line just above this


    #Plot percentages vs lifetimes
    if name == 'Lorenz':
        tname = 'Lorenz63'
    else:
        tname = name

    fig = plt.figure(figsize=fsize)
    ax = fig.add_subplot(111)
    s = 30
    s_out = 30

    for perc in percentages:
        c_lives = comp_lives_dict[perc]
        c_sizes = comp_sizes_dict[perc]
        l_lives = loop_lives_dict[perc]
        M = len(c_lives)
        K = len(l_lives)

        c_lives[c_lives==np.inf] = inf_val
        l_lives[l_lives==np.inf] = inf_val

        comp_pairs = [(perc, c_lives[n]) for n in range(min(M,num2show)) if c_sizes[n]>3]
        loop_pairs = [(perc, l_lives[n]) for n in range(min(K,num2show))]

        comp_pairs_outliers = [(perc, c_lives[n]) for n in range(min(M,num2show)) if c_sizes[n]<4]

        x_comps, y_comps = zip(*comp_pairs)
        if not(comp_pairs_outliers == []):
            x_comps_outliers, y_comps_outliers = zip(*comp_pairs_outliers)
        x_loops, y_loops = zip(*loop_pairs)
        x_loops = np.array(x_loops) + 1.0

        if percentages.index(perc) == 0:
            plt.scatter(x_comps, y_comps, marker='o', color='r', s=s, label='Components')
            if not(comp_pairs_outliers == []):
                plt.scatter(x_comps_outliers, y_comps_outliers, marker='*', color='r', s=s_out, alpha=0.4)
                #pass
            plt.scatter(x_loops, y_loops, marker='^', color='b', s=s, label='Loops')
        else:
            plt.scatter(x_comps, y_comps, marker='o', color='r', s=s)
            if not(comp_pairs_outliers == []):
                plt.scatter(x_comps_outliers, y_comps_outliers, marker='*', color='r', s=s_out, alpha=0.4)
               #pass
            plt.scatter(x_loops, y_loops, marker='^', color='b', s=s)

    
    ax.axhline(y=min_pers, color='k', linewidth=1)
    ax.axhline(y=0.35, color='k', linestyle='--', linewidth=1)
    ax.axhline(y=inf_val, color='k', linewidth=1)
    trans = transforms.blended_transform_factory(ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(0,inf_val, "inf", color="black", transform=trans, 
        ha="right", va="center")
    ax.text(0,min_pers, "min_pers", color="black", transform=trans, 
        ha="right", va="center", fontsize=8)
    ax.text(0,0.35, "gaussian", color="black", transform=trans, 
        ha="right", va="center", fontsize=8)
    
    old_ticks = ax.get_yticks()
    if name in ['CDV', 'Lorenz96']:
        new_ticks = old_ticks[2:-2]
    else:
        new_ticks = old_ticks[1:-1]

    ax.set_yticks(new_ticks)

    for t in ax.xaxis.get_major_ticks():
        t.label.set_fontsize(9)
    for t in ax.yaxis.get_major_ticks():
        t.label.set_fontsize(9)

    ax.margins(0.05)
    plt.xticks(list(range(10,100+10,10)))
    plt.ylim(ymax=1.05*inf_val)
    plt.legend(loc='best')
    if reverse_order:
        plt.title('%s bifiltration: %s longest lived (rev.)' % (tname, num2show))
    else:
        plt.title('%s bifiltration: %s longest lived cycles' % (tname, num2show))
    plt.xlabel('Density threshold (%)')
    plt.ylabel('Lifespan')
   

    #Saving/showing
    plt.tight_layout()
    if save:
        if reverse_order:
            if eofs:
                fname = '%s_EOFs_bifiltration_reversed_lifespans_%s-oldest_%s.%s' % \
                       (name, num2show, method, fmt)
            else:
                fname = '%s_bifiltration_reversed_lifespans_%s-oldest_%s.%s' % \
                       (name, num2show, method, fmt)
        else:
            if eofs:
                fname = '%s_EOFs_bifiltration_lifespans_%s-oldest_%s.%s' % \
                       (name, num2show, method, fmt)
            else:
                fname = '%s_bifiltration_lifespans_%s-oldest_%s.%s' % \
                       (name, num2show, method, fmt)
        fpath = '%s' % (figpath)
        figname = '%s/%s' % (fpath, fname)
        make_path(fpath)
        plt.savefig(figname, format=fmt, dpi=dpi)
        plt.close()
        print("(density_vs_lifetime) Figure saved as %s" % figname)
        print(100*'-')
    else:
        plt.show()
    



#Plot of density histogram with percentages marked
def density_histogram(name, percentages, max_edge=None, min_pers=None, sparse=0.7, pre_sparse=None,\
                      method='GaussianKDE', num_bins=120,\
                      identifier=None, reverse_order=None,\
                      save=True):

    print("(density_histogram) Generating histogram of density...")

    #Get the raw data and compute the density
    data = loader.load_raw(name)
    kde = stats.gaussian_kde(data)
    density = kde(data)


    #Prepare figure and plot
    percentiles = [np.percentile(density, percentages[n]) for n in range(len(percentages))]
    percentile_labels = [str(percentages[n]) for n in range(len(percentages))]

    fig = plt.figure(figsize=fsize)
    ax = fig.add_subplot(111)
    plt.hist(density, facecolor='g', density=True)
    plt.xlabel('Percentiles')
    plt.ylabel('Density')
    plt.title('%s density histogram' % name)
    plt.xticks(percentiles, percentile_labels)
    plt.tight_layout()

    if save:
        fname = '%s_density_histogram_%s.%s' % \
                   (name, method, fmt)
        fpath = '%s/%s' % (figpath, name)
        figname = '%s/%s' % (fpath, fname)
        make_path(fpath)
        plt.savefig(figname, format=fmt, dpi=dpi)
        plt.close()
        print("(density_histogram) Figure saved as %s" % figname)
        print(100*'-')
    else:
        plt.show()






