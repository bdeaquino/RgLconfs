#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:44:48 2019

@author: belisa
"""

import os
#from Property import occupancy_matrix

def get_files_with(files, secnames, paths):
    for i in paths:
        for j in os.listdir(i):
            fname = os.path.join(i,j)
            
            count = 0
            nsecname = len(secnames)
            for k in range(nsecname):
                if os.path.isfile(fname):
                    if secnames[k] in j:
                        count+=1
                    else:
                        break
                    
                    if count == nsecname:
                        files.append(fname)   
                        
def set_output_name(label, parameter, dt, tp=None, path=None, criterion=None, 
                    nhubs=None, selHubIDs=None, traj=None, only_hubs=False,
                    first_interval=False, fmt='txt', finfo=False, ffac=0.5, 
                    facbothdir=0.8, only_X=False, only_Y=False):
        if criterion:
            if not nhubs:
                print('Error: missing nhubs value.')
                
            if criterion == 'hub_data' or criterion == 'between_hubs' or criterion == 'between_hubs_ordered':
                sellen = len(selHubIDs)
                if sellen == nhubs and (criterion == 'between_hubs' or criterion == 'hub_data'):
                    hubname='_all'
                else:
                    hubname=''
                    for n in range(sellen):
                        hubname+='_'+str(selHubIDs[n])
                        
            elif criterion == 'from_hub' or criterion == 'to_hub' or criterion == 'folding_hub':
                hubname=str(selHubIDs[0])
        
        if path:
            txtname=path+label+'_'+parameter+'_dt'+str(dt)
        else:
            txtname=label+'_'+parameter+'_dt'+str(dt)
        
        if tp:
            txtname+='_tp'+str(tp)
        
        if traj is not None:
            txtname+='_traj'+str(traj)
        
        if criterion:
            txtname+='_nhubs'+'_'+str(nhubs)+'_'+criterion+hubname
            
        if only_hubs:
            txtname+='_only_hubs'
            
        if first_interval:
            txtname+='_first_interval'
            
        if finfo and None not in [ffac, facbothdir]:
           txtname+='_finfo_'+str(ffac)+'_'+str(facbothdir)
           
        if only_X:
            txtname+='_Rg'
        
        if only_Y:
            txtname+='_L'
            
        txtname+='.'+fmt
    
        return txtname