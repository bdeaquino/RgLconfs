#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:46:13 2019

@author: belisa
"""

import numpy as np
import pylab as pl

class RgL(object):
    def __init__(self, fin, lim, dt=1, 
                 idx=0, skip=0, sigma_fac=2, thermo_fac=1,
                 energy=False, contacts=False,
                 xmin=None, xmax=None, len_xbin=None, 
                 ymin=None, ymax=None, len_ybin=None,
                 criterion=None, nhubs=None, hubs=None, 
                 selHubIDs=None, traj=None, nnat=None,
                 compute_distribution=True, compute_histogram=False,
                 compute_average=False, compute_thermo=False, compute_tmedian=False,
                 compute_occupancy=True, compute_transitions=True, 
                 compute_energy=False, compute_contacts=False, only_hubs=False,
                 first_interval=False, only_X=False, only_Y=False):
        
        super(RgL, self).__init__()
        
        self.fin=fin
        self.lim=lim
        self.dt=dt
        self.idx=idx
        self.skip=skip
        self.sigma_fac=sigma_fac
        self.thermo_fac=thermo_fac
        self.traj=traj  
        self.nnat=nnat
        
        self.only_X=only_X
        self.only_Y=only_Y
        
        self.xmin=xmin
        self.xmax=xmax
        self.len_xbin=len_xbin
        self.ymin=ymin
        self.ymax=ymax
        self.len_ybin=len_ybin
        
        if None in [xmin, xmax, len_xbin, ymin, ymax, len_ybin] and compute_distribution:
            self.compute_data(compute_distribution=True,
                              compute_occupancy=False, compute_transitions=False,
                              sigma_fac=sigma_fac, energy=energy, contacts=contacts)

        if None not in [self.xmin, self.xmax, self.len_xbin, self.ymin, self.ymax, self.len_ybin]:   
            self.bins_x, self.xnbins = bin_params(self.xmin, self.xmax, self.len_xbin)
            self.bins_y, self.ynbins = bin_params(self.ymin, self.ymax, self.len_ybin)
        
        self.compute_occupancy=compute_occupancy
        self.compute_transitions=compute_transitions
        self.compute_distribution=compute_distribution
        self.compute_histogram=compute_histogram
        self.compute_average=compute_average
        self.compute_thermo=compute_thermo
        self.compute_tmedian=compute_tmedian
        self.compute_energy=compute_energy
        self.compute_contacts=compute_contacts
        
        if not criterion and not compute_energy and not compute_contacts:
            self.compute_data(compute_histogram=compute_histogram,
                              compute_average=compute_average,
                              compute_thermo=compute_thermo,
                              compute_tmedian=compute_tmedian,
                              compute_occupancy=compute_occupancy,
                              compute_transitions=compute_transitions,
                              sigma_fac=sigma_fac, 
                              energy=energy, contacts=contacts)
            
        elif criterion or compute_energy or compute_contacts:    
            self.criterion=criterion
            self.only_hubs=only_hubs
            self.first_interval=first_interval
            
            self.hubs=hubs
            if not nhubs:
                if hasattr(self.hubs, '__len__'):
                    self.nhubs=len(hubs)
            else:
                self.nhubs=nhubs
                
            self.selHubIDs=selHubIDs
            self.selHubIDs = self.check_selHubs()
            if hasattr(self.selHubIDs, '__len__'):
                self.nselHubs=len(self.selHubIDs)
            else:
                self.nselHubs=None
                
            self.compute_data(compute_transitions=False,
                              compute_tmedian=compute_tmedian,
                              sigma_fac=sigma_fac, 
                              energy=energy, contacts=contacts)
            
            [self.hubs, self.idHubs, self.occHubs] = self.determine_hubs()  
            print('Hubs: '+str(self.hubs))
            
            self.selHubs = [self.hubs[k] for k in self.selHubIDs]
            
            self.compute_data(criterion=criterion, only_hubs=only_hubs,
                              traj=traj, first_interval=first_interval,
                              compute_occupancy=compute_occupancy,
                              compute_transitions=compute_transitions,
                              sigma_fac=sigma_fac, 
                              compute_energy=compute_energy,
                              compute_contacts=compute_contacts,
                              energy=energy, contacts=contacts)
        
        if compute_occupancy:            
            self.occupancyMatrix = self.occupancy_matrix()
    
    def compute_data(self, criterion=None, traj=None, 
                     compute_distribution=False, sigma_fac=2,
                     compute_occupancy=True, compute_transitions=True,
                     only_hubs=False, first_interval=False,
                     compute_average=False, compute_thermo=False,
                     compute_histogram=False, compute_tmedian=False,
                     compute_energy=False, compute_contacts=False,
                     energy=False, contacts=False):
        """
        possible criterion values: 'hub_data', 'between_hubs', 'between_hubs_ordered', 'from_hub', 'to_hub', 'folding', 'folding_hub'
        """
        
        if compute_average:
            print("Calculating averages.")
            avX = []
            avY = []
            
        if compute_thermo:
            print("Computing thermodynamic properties.")
            avQ = []
            avP0 = []
            avnc = []
            
        if compute_tmedian:
            if not self.nnat:
                print("Error: compute_tmedian requires number of native contacts (nnat).")
            else:
                print("Computing median folding time.")
                tf = []
            
        if compute_histogram:
            print("Computing histogram.")
            H = []
        
        if compute_distribution:
            print("Calculating binning parameters.")
            difX = []
            difY = []
            xmin = None
            xmax = None
            ymin = None
            ymax = None
            
        interval = None

        occupancy = []
        scanned = []
        transitions = []
        f = []
            
        if criterion=='hub_data':
            print('Computing data information for '+str(self.nhubs)+' hubs.')
            output = []
        
            if not hasattr(self.hubs, '__len__'):
                print("Error: missing declaration of 'hubs'.")
        
        elif criterion=='folding':
            self.ntraj_fold = 0
            self.ntraj_notfold = 0
        
        elif criterion=='to_hub':
            ntraj_nottohub = 0
            th = []
        
        elif criterion=='between_hubs_ordered':
            ntraj_ord = 0
            self.direct_trans = 0
            self.notdirect_trans = 0
            self.ntrans = 0
            self.t_trans_dir = []
            self.t_trans_notdir = []
                    
        if traj is not None:
            itraj = traj[0]
            ittraj = 0
        else:
            itraj = 0
            
        i = 0
        
        if compute_energy:
            self.Elist = [[] for y in range(self.nselHubs)]
        
        if compute_contacts:
            self.nclist = [[] for y in range(self.nselHubs)]
        
        with open(self.fin) as fp:
            for j, line in enumerate(fp):            
                if i == 0:
                    i+=1
                    start = self.idx+(self.lim)*itraj
                    stop = start+self.lim

                    t = []
                    x = []
                    y = []
                    
                    if compute_energy:
                        Etraj = []
                        
                    if compute_contacts or criterion == 'folding' or criterion == 'folding_hub' or compute_thermo or (compute_tmedian and self.nnat):
                        nc = []
            
                if j >= start:
                    if (j-start)%self.dt == 0:
                        if energy == True and contacts == True:
                            tt, Et, nct, xt, yt = line.split()
                        elif energy == True and contacts == False:
                            tt, Et, xt, yt = line.split()
                        elif energy == False and contacts == True:
                            tt, nct, xt, yt = line.split()
                        else:
                            tt, xt, yt = line.split()
                        t.append(float(tt))
                        x.append(float(xt))
                        y.append(float(yt))
                        
                        if compute_contacts or criterion == 'folding' or criterion == 'folding_hub' or compute_thermo or (compute_tmedian and self.nnat):
                            nc.append(float(nct))
                        if compute_energy:
                            Etraj.append(float(Et))
                        
                if j == stop-1:
                    if self.only_X:
                        y = np.zeros(len(y))
                    elif self.only_Y:
                        x = np.zeros(len(x))
                        
                    if compute_histogram:
                        h, xedges, yedges = self.create_histogram(x, y)
                        if len(H) == 0:
                            H = np.array(h)
                        else:
                            H += np.array(h)
                        
                    if compute_average:
                        avX.append(np.mean(x))
                        avY.append(np.mean(y))
                        
                    if compute_tmedian and self.nnat:
                        tft, inc = tstep(self.nnat, nc, t)
                        if tft != -1:
                            tf.append(tft)
                        
                    if compute_thermo:
                        ncmax = max(nc)
                        Q, p0, ncnts = calc_QnP0(nc, ncmax, skip=self.skip, fac=self.thermo_fac)
                        avQ.append(Q)
                        avP0.append(p0)                      
                        avnc.append(ncnts)
                    
                    if compute_distribution:
                        difX.extend(np.diff(x))
                        difY.extend(np.diff(y))
                        
                        xmin, xmax = get_minmax(x, xmin, xmax)
                        ymin, ymax = get_minmax(y, ymin, ymax)
                                            
                    compute = True
                    
                    if compute_occupancy or compute_transitions:
                        xinds, bins_x, xnbinspri = digitize_vec(x, self.len_xbin, self.xmin, self.xmax)
                        yinds, bins_y, ynbins = digitize_vec(y, self.len_ybin, self.ymin, self.ymax)
                        inds = list(zip(xinds,yinds))                    
                        bins = list(bins_x) + list(bins_y)                            
                    
                    if criterion or compute_energy:
                        for hID in self.selHubIDs:
                            if self.hubs[hID] not in inds:
                                print('Error: Some of the selected hubs were not reached. Skipping trajectory.')
                                compute=False
                                
                    if compute_energy:
                        for k1,x1 in enumerate(self.selHubs):
                           self.Elist[k1].extend([(itraj, k2, Etraj[k2]) for k2,x2 in enumerate(inds) if x2 == x1])
                    
                    if compute_contacts:
                        for k1,x1 in enumerate(self.selHubs):
                           self.nclist[k1].extend([(itraj, k2, nc[k2]) for k2,x2 in enumerate(inds) if x2 == x1])
                            
                    if criterion:          
                        if compute:
                            if criterion == 'hub_data':
                                self.dataHubs = self.t_hub_list(output, itraj, inds, t)
                                self.f_most_probable_path, self.hubs_most_probable_path = self.most_probable_path()
                            
                            elif criterion == 'between_hubs_ordered':
                                if len(self.selHubIDs) != 2:
                                    print("Error: criterion demands len(selHubIDs)=2")
                                    break
                                    
                                else:
                                    interval = self.tscan_hubs_ordered(inds,
                                                                       first_interval=first_interval)   

                                    if len(interval) < 1:
                                        ntraj_ord+=1
                                        print('Inexistent route between ordered hubs. Skipping '+str(ntraj_ord)+' trajectories.')
                                        compute=False
                                    
                                    else:
                                        direct, notdirect, t_trans_dir, t_trans_notdir = count_direct_transitions(interval, inds, t)
                                        self.direct_trans+=direct
                                        self.notdirect_trans+=notdirect
                                        
                                        self.t_trans_dir.extend(t_trans_dir)
                                        self.t_trans_notdir.extend(t_trans_notdir)
                            
                            elif criterion == 'between_hubs':
                                if len(self.selHubIDs) < 2:
                                    print("Error: criterion demands len(selHubIDs)>=2")
                                    break
                                
                                else:   
                                    interval = self.tscan_hubs(inds, 
                                                               first_interval=first_interval)
                                    
                            elif criterion == 'from_hub' or criterion == 'to_hub' or criterion == 'folding_hub':                            
                                if len(self.selHubIDs) != 1:
                                    print("Error: criterion demands len(selHubIDs)=1.")
                                    break
                                    
                                else:
                                    k = self.selHubIDs[0]
                                    if criterion == 'from_hub':                    
                                        interval = self.tscan_hubs(inds, hbeg=self.hubs[k],
                                                                   first_interval=first_interval)
                                    elif criterion == 'to_hub':
                                        if self.only_hubs:
                                            print("Error: only_hubs cannot be chosen with this criterion.")
                                        else:
                                            ibeg = 0
                                            th_traj, iend = tstep(self.hubs[k], inds, t)
                                            if th == -1:
                                                compute=False
                                                ntraj_nottohub+=1
                                                print("Hub not reached. Skipping trajectory ("+str(ntraj_nottohub)+").")
                                            else:
                                                th.append(th_traj)
                                                interval = [(ibeg, iend)]   
                                            
                                    elif criterion == 'folding_hub':
                                        if not self.nnat or self.nnat <= 0:
                                            print('Error: criterion demands a positive value for nnat.')
                                            break                                    
                                        
                                        else:
                                            if self.nnat not in nc:
                                                print("Number of native contacts not reached. Skipping trajectory.")
                                                compute=False
                                            else:
                                                ibeg = 0
                                                self.tf, inc = tstep(self.nnat, nc, t)
                                                hubt = inds[inc]
                                                self.th, iend = tstep(hubt, inds, t)
                                                interval = [(ibeg, iend)]                                    
                                                
                            elif criterion == 'folding':
                                if not self.nnat:
                                    print("Error: compute_tmedian requires number of native contacts (nnat).")
                                else:
                                    tf_traj, inc = tstep(self.nnat, nc, t)
                                    if tf_traj != -1:                                
                                        ibeg = 0                                    
                                        iend = inc
                                        interval = [(ibeg, iend)]
                                        self.ntraj_fold+=1
                                    
                                    else:
                                        compute=False
                                        self.ntraj_notfold+=1
                                        print("Folding not reached. Skipping trajectory ("+str(self.ntraj_notfold)+").")
                                

                    if compute:
                        if compute_occupancy:                            
                            occupancy, scanned = self.compute_property(inds, interval=interval, 
                                                                           prop=occupancy, 
                                                                           prop2=scanned,
                                                                           only_hubs=only_hubs,
                                                                           propname='occupancy')
                        
                        if compute_transitions:                            
                            f, transitions = self.compute_property(inds, interval=interval, 
                                                                       prop=f, 
                                                                       prop2=transitions, 
                                                                       only_hubs=only_hubs,
                                                                       propname='transitions')
                    i=0
                    
                    if traj is not None:
                        ittraj+=1
                        if ittraj < len(traj):
                            itraj = traj[ittraj]
                        else:
                            break
                        
                    else:
                        itraj+=1
        
        self.ntraj = itraj
        
        if criterion == 'hub_data':
            self.t_hubs, self.frac_traj_hubs, self.frac_traj_fh =\
            self.time_data_hubs()
            
        if criterion == 'between_hubs_ordered':
            self.ntrans = self.direct_trans + self.notdirect_trans
            
            self.t_trans = self.t_trans_dir+self.t_trans_notdir
            self.tmedian_trans = np.median(self.t_trans)
            self.tmedian_trans_dir = np.median(self.t_trans_dir)
            self.tmedian_trans_notdir = np.median(self.t_trans_notdir)
        
        if criterion == 'to_hub':
            if len(th):
                self.ntraj_tohub = len(th)
                self.th = np.median(th)
            else:
                self.ntraj_tohub = 0
                self.th = self.lim
        
        if compute_histogram:
            self.H = H
            self.xedges = xedges
            self.yedges = yedges
            
        if compute_average:
            self.avX = np.mean(avX)
            self.avY = np.mean(avY)
            
        if compute_tmedian and self.nnat:
            if len(tf):
                self.ntraj_fold = len(tf)
                self.tf = np.median(tf)
            else:
                self.ntraj_fold = 0
                self.tf = self.lim
            
        if compute_thermo:
            self.Q = np.mean(avQ)
            self.P0 = np.mean(avP0)
            self.nc = np.mean(avnc)
            
        if compute_distribution:            
            self.stdX, self.len_xbin = displ_std(difX, sigma_fac)
            self.stdY, self.len_ybin = displ_std(difY, sigma_fac)
            self.xmin = xmin
            self.xmax = xmax
            self.difX = difX
            self.ymin = ymin
            self.ymax = ymax
            self.difY = difY
            print('xmin ='+str(self.xmin)+', xmax='+str(self.xmax)+', len_xbin='+str(self.len_xbin))
            print('ymin ='+str(self.ymin)+', ymax='+str(self.ymax)+', len_ybin='+str(self.len_ybin))
            
        if compute_occupancy:
            self.occupancyBins = occupancy
            self.scannedBins = scanned
            self.posBins = bins
            
        if compute_transitions:
            self.transitionsBins = transitions
            self.fTransBins = f
            
            if not compute_occupancy:
                self.posBins = bins
    
    def compute_property(self, inds, interval=None, prop=None, prop2=None, 
                         only_hubs=False, propname='None'):                
        if prop==None or prop2==None:
            print('Please specify scanned and occupancy arrays.')
        
        if interval == None:
            interval = []
            ibeg = 0
            
            if propname == 'transitions':
                iend = len(inds)-1
                
            else:
                iend = len(inds)
            interval.append((ibeg, iend))
            
        for (ibeg, iend) in interval:
            if iend != -1:
                if only_hubs:
                    ihubs = []
                    for n in inds[ibeg:iend+1]:
                        if n in self.hubs:
                            ihubs.append(inds.index(n))
                            
                    print(ihubs)
            
                    if propname == 'transitions':
                        nvec = range(len(ihubs)-1)
                    else:
                        nvec = range(len(ihubs))
                else:
                    nvec = range(ibeg,iend)
                    
                for n in nvec:
                    if only_hubs:
                        k1 = ihubs[n]
                    else:
                        k1 = n
                    
                    if propname == 'transitions':
                        if only_hubs:
                            k2 = ihubs[n+1]
                        else:
                            k2 = n+1
                        (a,b) = inds[k1]
                        (c,d) = inds[k2]
                         
                        if a != c or b != d:        
                             update_vectors_n_freqs((a, b, c, d), prop, prop2)
                         
                    else:
                        update_vectors_n_freqs(inds[k1], prop, prop2)
                
        return prop, prop2

    def occupancy_matrix(self):
        H = np.zeros((self.ynbins,self.xnbins))
        
        for (n,m) in self.scannedBins:
            idx = self.scannedBins.index((n,m))
            H[m-1][n-1] = self.occupancyBins[idx]
            
        return H

    def plot_occ_matrix(self):
        self.occupancy_matrix()
            
        fig,sbpl=pl.subplots(1,1)
        H = np.ma.array(self.occupancyMatrix, mask=np.isinf(np.log10(self.occupancyMatrix)))
        X, Y = np.meshgrid(self.bins_y, self.bins_x)
        im = sbpl.pcolormesh(Y, X, H.T, cmap='nipy_spectral_r')
        
        return im
    
    def create_histogram(self, x, y): 
        bins = bins_histogram(self.xmin, self.xmax, self.ymin, self.ymax)
            
        H, xedges, yedges = np.histogram2d(x, y, bins=bins, density=1)
        H = H.T#/np.max(H)
        
        return H, xedges, yedges
    
    def check_selHubs(self):
        try:
            self.selHubIDs
        except NameError:
            self.selHubIDs=None
        
        if self.selHubIDs:      
            for sh in self.selHubIDs:
                if sh > self.nhubs-1:
                    print('Error: selHubIDs values must be between 0 and nhubs-1.')
                    if self.criterion == 'from_hub' or self.criterion == 'to_hub':
                        print('Setting selHubIDs = [0].')
                        self.selHubIDs = [0]
                    if self.criterion != 'from_hub' or self.criterion != 'to_hub':
                        print('Setting selHubIDs = '+str(np.arange(0,self.nhubs))+'.')
                        self.selHubIDs = np.arange(0,self.nhubs)       
        
        if self.criterion == 'between_hubs':
            if self.selHubIDs:
                if not hasattr(self.selHubIDs, '__len__') or len(self.selHubIDs) < 2:
                    print('Error: criterion demands len(selHubIDs) >= 2.')
            else:
                print('Setting selHubIDs = np.arange(0,nhubs).')
                self.selHubIDs = np.arange(0,self.nhubs)
                
        elif self.criterion == 'between_hubs_ordered':
            if self.selHubIDs:
                if not hasattr(self.selHubIDs, '__len__') or len(self.selHubIDs) != 2:
                    print('Error: criterion demands len(selHubIDs) = 2.')
            else:
                print('Setting selHubIDs = [0,nhubs-1].')
                self.selHubIDs = [0,self.nhubs-1]
            
        elif self.criterion == 'from_hub' or self.criterion == 'to_hub' or self.criterion == 'folding_hub':
            if self.selHubIDs:
                if not hasattr(self.selHubIDs, '__len__') or len(self.selHubIDs) != 1:
                    print('Error: criterion demands a single value for selHubIDs.')
            else:
                print('Setting selHubIDs = [0].')
                self.selHubIDs = [0]
                
        elif self.criterion == 'hub_data':
            if not self.selHubIDs:
                print('Setting selHubIDs = np.arange(0,nhubs).')
                self.selHubIDs = np.arange(0,self.nhubs)
                
        else:
             if not self.selHubIDs:
                print('Setting selHubIDs = np.arange(0,nhubs).')
                self.selHubIDs = np.arange(0,self.nhubs)
            
                
        return self.selHubIDs
    
    def determine_hubs(self):
        if not self.hubs:
            occ = np.sort(self.occupancyBins)    
            occ = occ[::-1]
            idHubs = [self.occupancyBins.index(k) for k in occ[0:self.nhubs]]
            occHubs = occ[0:self.nhubs]/sum(occ)
            hubs = [self.scannedBins[k] for k in idHubs]
            
        else:
            hubs = self.hubs
            idHubs = [k for k,x in enumerate(self.scannedBins) if x in hubs]
            occHubs = [x for k,x in enumerate(self.occupancyBins) if k in idHubs]
            
        return [hubs, idHubs, occHubs]
    
    def tscan_hubs(self, inds, hbeg=None, first_interval=False):        
        interval = []
        
        if hbeg:
            ibeg = inds.index(hbeg)
        else:
            ibeg = 0
            
        if self.criterion == 'from_hub':
            hubs = self.hubs
        else:
            hubs = self.selHubs
        
        hubs_temp = []
        i = 0
        for k,x in enumerate(inds[ibeg:]):
            if x in self.selHubs and i == 0:
                ibeg_t = k
                i = k
                hubs_temp.append(x)
            elif x in hubs and x not in hubs_temp and k >= i and i > 0:
                hubs_temp.append(x)
                i = k
                
                if len(hubs_temp) == len(hubs):
                    iend_t = k+1
                    interval.append((ibeg_t, iend_t))
                    i = 0
                    if first_interval:
                        break
            elif x not in hubs and k == len(inds[ibeg:])-1 and not len(interval):
                print("Route of selected hubs not achieved by trajectory. Skipping trajectory.")
                interval.append((0, -1))
                break
        
        return interval
    
    def tscan_hubs_once(self, inds, ibeg=None):        
        if ibeg:
            inds = inds[ibeg:]
        else:
            ibeg = 0
            
        ik = []
        for k,x in enumerate(self.hubs):
            if x in inds:
                ik.append(inds.index(x))
            else:
                ik.append(-1)
    
        iksel = [x for x in ik if x != -1]
        ibeg += min(iksel)
        iend = ibeg+max(iksel)
        
        return ibeg, iend
    
    def tscan_hubs_ordered(self, inds, ibeg=None, first_interval=False):
        interval = []
        
        if not ibeg:           
            ibeg = 0

        if first_interval:
            if self.selHubs[0] in inds[ibeg:]:
                k0 = inds[ibeg:].index(self.selHubs[0])
            if self.selHubs[-1] in inds[ibeg:]:
                k1 = inds[ibeg:].index(self.selHubs[-1])
                
            if k1 > k0:
                ibeg_t = k0
                iend_t = k1
                interval.append((ibeg_t, iend_t))
                
        else:
            i=0
            for k,x in enumerate(inds[ibeg:]):
                if x == self.selHubs[0] and i == 0:
                    ibeg_t = k
                    i = k
                elif x == self.selHubs[-1] and k >= i and i > 0:
                    iend_t = k
                    interval.append((ibeg_t, iend_t))
                    i = 0
            
        return interval
    
    def t_hub_list(self, output, itraj, inds, t):
        for n in range(self.nhubs):
            if self.hubs[n] in inds:                                
                ilim = inds.index(self.hubs[n])
                output.append([itraj, n, self.hubs[n], t[ilim]])
            else:
                output.append([itraj, n, self.hubs[n], -1])
        
        return output
        
    def time_data_hubs(self):
        t_hubs = [[-1 for x in range(self.nhubs)] for y in range(self.nhubs+1)] 
        frac_traj_hubs = []
        ntraj = int(len(self.dataHubs)/self.nhubs)
        thubs = []
        for n in range(self.nhubs):
            tt = [d for [a,b,c,d] in self.dataHubs if b == n and d != -1]
            ntrajt = len(tt)
            
            thubs.append(np.median(tt))
            frac_traj_hubs.append(ntrajt/ntraj)
        
        t_hubs[0] = thubs 
        
        for n1 in range(self.nhubs):        
            tthubs = [[-1 for x in range(ntraj)] for y in range(self.nhubs)]
            for n2 in range(ntraj):
                tt = [d for [a,b,c,d] in self.dataHubs if a == n2]
                ito = list(np.argsort(tt))
                
                tthubs[ito.index(n1)][n2] = tt[n1]
                
            for n2 in range(self.nhubs):        
                tsel = [a for a in tthubs[n2] if a != -1]
                    
                if len(tsel):
                    tt2 = np.median(tsel)
                    
                else:
                    tt2 = -1
                
                t_hubs[n2+1][n1] = tt2
            
                #break
        
        ntraj_fh = np.zeros(self.nhubs)
        
        for n1 in range(ntraj):
            tt1 = []
            for [a,b,c,d] in self.dataHubs:
                if a == n1:
                    tt1.append(d)
            
            for n2 in range(self.nhubs):
                if tt1[n2] == min(tt1):
                    ntraj_fh[n2] += 1
                    
        frac_traj_fh = ntraj_fh/ntraj
        
        return t_hubs, frac_traj_hubs, frac_traj_fh
    
    def traj_hub_list(self):
        traj_list = []
        for n1 in range(self.ntraj):
            tt1 = []
            for [a,b,c,d] in self.dataHubs:
                if a == n1:
                    tt1.append(d)
                    
            hub_list = list(np.argsort(tt1))
            traj_list_t = [n1]
            traj_list_t.extend(hub_list)
            
            traj_list.append(traj_list_t)
        
        return traj_list
    
    def most_probable_path(self, index=None):
        traj_hubs = np.array(self.traj_hub_list()).T
        idcs = list(range(self.ntraj))
        hub_traj = []
        
        for n in range(self.nhubs):
            if index != None and n == 1:
                hub = np.bincount(traj_hubs[n+1][idcs]).argsort()[::-1][index]
            else:
                hub = np.bincount(traj_hubs[n+1][idcs]).argmax()
                
            hub_traj.append(hub)
            idcst = [i for i,j in enumerate(traj_hubs[n+1][idcs]) if j == hub]
            idcs = np.array(idcs)[idcst]
            
        rate = len(idcs)/self.ntraj
        return rate, hub_traj, idcs

    def prop_to_file(self, propname, files):
        """
        #occ_var = [occupancy, scanned, xmin, xmax, len_xbin, ymin, ymax, len_ymin]
        occ_var = H
        bins_var = bins
        trans_var = [transitions, f]
        hubs_var = [nhubs, hubs, occ, t_hubs, frac_traj_hubs, frac_traj_fh]
        """
        length = len(propname)
        
        print('Generating output files:')
        
        for i in range(length):
            pname = propname[i]
    
            if pname == 'occupancy':
                #prop, bins_x, bins_y = occupancy_matrix(*occ_var)
                prop = self.occupancyMatrix
                fmt="%.18e"
                    
            elif pname == 'bins':
                prop = self.posBins
                fmt="%.18e"
                
            elif pname == 'transitions':
                prop = [(int(a), int(b), int(c), int(d), e) for (a, b, c, d), e in zip(self.transitionsBins, self.fTransBins)]
                fmt="%d"
                
            elif pname == 'hub_data':
                if self.compute_occupancy:
                    hubs_var = [self.nhubs, self.hubs, self.occHubs, self.t_hubs, self.frac_traj_hubs, self.frac_traj_fh]
                else:
                    hubs_var = [self.nhubs, self.hubs, self.t_hubs, self.frac_traj_hubs, self.frac_traj_fh]
                prop = []
                for n in range(len(hubs_var)):
                    if n == 0:
                        prop.append(list(np.arange(hubs_var[n])))
                    elif n == 3:
                        for line in hubs_var[n]:
                            prop.append(list(line))
                    else:
                        prop.append(list(hubs_var[n]))
                fmt="%s"
    
            else:
                print("Error: propname variables must be 'occupancy', 'transitions', 'bins' or 'hub_data'.")
            
            txtname = files[i]
            
            print(txtname)
            np.savetxt(txtname, prop, fmt=fmt)

###############################################################################

def get_minmax(x, xmin=None, xmax=None):    
    if None in [xmin, xmax]:
        xmin = min(x)
        xmax = max(x)
    else:
        xmin = min(xmin, min(x))
        xmax = max(xmax, max(x))
        
    xmin = int(xmin//1)
    xmax = int(-(-xmax//1))
    
    return xmin, xmax

def displ_std(difX, sigma_fac):
    stdX = np.std(difX)        
    len_xbin = sigma_fac*stdX
    
    return stdX, len_xbin

def bin_params(xmin, xmax, len_xbin):
    xnbins = int(-(-(xmax-xmin)//len_xbin))
    bins_x = np.arange(xmin, xmax, len_xbin)
    
    return bins_x, xnbins

def bins_histogram(xmin, xmax, ymin, ymax):
    xlen = xmax-xmin
    ylen = ymax-ymin
    
    if xlen < ylen:
        rat = int(float(ylen)/xlen)
        xedges = np.linspace(xmin,xmax,xlen*rat)
        yedges = np.linspace(ymin,ymax,ylen)
    else:
        rat = int(float(xlen)/ylen)
        xedges = np.linspace(xmin,xmax,xlen)
        yedges = np.linspace(ymin,ymax,ylen*rat)
    bins=(xedges, yedges)
    
    return bins
    
def digitize_vec(x, len_xbin, xmin=None, xmax=None):
    x = np.array(x)
    
    if xmin == None:
        xmin = min(x)
    if xmax == None:
        xmax = max(x)
    
    bins_x, xnbins = bin_params(xmin, xmax, len_xbin)
    xinds = np.digitize(x, bins_x)
    
    return xinds, bins_x, xnbins

def update_vectors_n_freqs(a, fvec, vec):
        if a not in vec:
            vec.append(a)
            fvec.append(1)
        else:
            idx = vec.index(a)
            fvec[idx] += 1
            
def calc_QnP0(nc, ncmax, skip=0, fac=1):    
    nc = list(map(int, nc[skip:-1]))
        
    if fac == 1:
        P0 = float(list(nc).count(ncmax))/len(nc)
            
    else:
        p0 = [x for x in nc if x >= fac*ncmax]
        P0 = float(len(p0))/len(nc)
        
    Q = np.mean(nc)/ncmax
    ncnts = np.mean(nc)
        
    return Q, P0, ncnts

def tstep(ntgt, nc, t):
    if ntgt in nc:
        inc = nc.index(ntgt)
        return t[inc], inc
    else:
        return -1, len(t)-1
    
def count_direct_transitions(interval, inds, t):
        direct = 0
        notdirect = 0
        t_dir = []
        t_notdir = []
        
        for (a,b) in interval:
            for i in range(a,b+1):
                if i == a:
                    hub_init = inds[i]
                    
                if inds[i] != hub_init:
                    if i != b:
                        notdirect+=1
                        t_notdir.append(t[b]-t[a])
                        break
                    else:
                        direct+=1
                        t_dir.append(t[b]-t[a])
        
        return direct, notdirect, t_dir, t_notdir
