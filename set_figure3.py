#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:05:23 2019

@author: belisa
"""

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl

def init(ncols, nrows, hspc=None, wspc=None, aspect=1, d=1, hsize=None, wsize=None, legend=None, lfac=None,
         bottom=None, left=None, right=None, top=None):    
    if hspc != None:
        hfac = hspc*(nrows-1)
    else:
        hfac = 0
    if wspc != None:
        wfac = wspc*(ncols-1)
    else:
        wfac = 0
        
    x = float(d*ncols + wfac)
    y = float(d*aspect*nrows + hfac)
    
    if legend != None:
        if lfac == None:
            print('Please specify legend size.')
        else:
            if legend == 'top' or legend == 'bottom':
                y+=lfac
            elif legend == 'left' or legend == 'right':
                x+=lfac
    
    if None not in (hsize, wsize):
        figsize = (wsize, hsize)
    else:
        figsize = plt.figaspect(y/x)
        print(figsize)
    fontsize=15*d
        
    fig, sbpl = pl.subplots(nrows, ncols, figsize=figsize)
    
    if bottom != None:
        fig.subplots_adjust(bottom=bottom)
    if left != None:
        fig.subplots_adjust(left=left)
    if right != None:
        fig.subplots_adjust(right=right)
    if top != None:
        fig.subplots_adjust(top=top)
    if hspc != None:
        fig.subplots_adjust(hspace=hspc)
    if wspc != None:
        fig.subplots_adjust(wspace=wspc)
    
    return fig, sbpl, fontsize

def set_subplot(sbpl, nx, ny, ix, iy):
    if nx < 1 or ny < 1:
        print('Wrong assignment of number of rows or columns.')
    elif nx == 1 and ny == 1:
        ax = sbpl
    elif nx == 1 and ny > 1:
        ax = sbpl[iy]
    elif nx > 1 and ny == 1:
        ax = sbpl[ix]
    else:
        ax = sbpl[iy, ix]
        
    return ax
    
def custom_subplot(sbpl, xmin=None, xmax=None, ymin=None, ymax=None, aspect=None,
                   xmin_tick=None, xmax_tick=None, xstep_tick=None, 
                   ymin_tick=None, ymax_tick=None, ystep_tick=None, 
                   top=False, right=False, left=True, bottom=True,
                   labelbottom=True, labelleft=True, labelright=False, labeltop=False,
                   set_spines=False, axes=['top', 'bottom', 'left', 'right'], 
                   fontsize=12, linewidth=2, ticklength=5,
                   xlabel=None, ylabel=None, labelpad=None):
    
    if xmin != None and xmin_tick == None:
        xmin_tick=xmin
    if xmax != None and xmax_tick == None:
        xmax_tick=xmax
    if xmin != None and xmax != None and xstep_tick == None:
        xstep_tick=5
    if ymin != None and ymin_tick == None:
        ymin_tick=ymin
    if ymax != None and ymax_tick == None:
        ymax_tick=ymax
    if ymin != None and ymax != None and ystep_tick == None:
        ystep_tick=5
    
    if aspect != None:
        sbpl.set_aspect(aspect)
    
    if None not in [xmin, xmax]:
        sbpl.set_xlim(xmin,xmax)
    if None not in [ymin, ymax]:
        sbpl.set_ylim(ymin,ymax)
    
    if xlabel != None:
        sbpl.set_xlabel(xlabel, fontsize=fontsize, labelpad=labelpad)
    if ylabel != None:
        sbpl.set_ylabel(ylabel, fontsize=fontsize, labelpad=labelpad)
    
    if None not in [xmin_tick, xmax_tick, xstep_tick]:
        sbpl.set_xticks(np.arange(xmin_tick, xmax_tick, xstep_tick))
    if None not in [ymin_tick, ymax_tick, ystep_tick]:
        sbpl.set_yticks(np.arange(ymin_tick, ymax_tick, ystep_tick))
    
    sbpl.tick_params(axis='both', top=top, right=right, left=left, bottom=bottom, 
                     labeltop=labeltop, labelleft=labelleft, labelright=labelright, labelbottom=labelbottom, 
                     direction='in', width=linewidth, length = ticklength, labelsize = fontsize)
            
    for axis in axes:
        sbpl.spines[axis].set_linewidth(linewidth)

def color_axes(sbpl, axis=None, naxis=None, color='k', label=False, labelname=None, fontsize=None): 
    sbpl.tick_params(axis=axis, colors=color)
    sbpl.spines[naxis].set_color(color)
    
    if label == True and labelname != None:
            if axis == 'y':
                if fontsize == None:
                    sbpl.set_ylabel(labelname, color=color)
                else:
                    sbpl.set_ylabel(labelname, color=color, fontsize=fontsize)
            else:
                if fontsize == None:
                    sbpl.set_xlabel(labelname, color=color)
                else:
                    sbpl.set_xlabel(labelname, color=color, fontsize=fontsize)
        
def custom_colorbar(fig, im, x, y, xlen, ylen, orientation='horizontal',
                    min_tick=None, max_tick=None, step_tick=None,
                    fontsize=12, linewidth=2, ticklength=7):
    cb_sbpl = fig.add_axes([x, y, xlen, ylen])
    cbar = fig.colorbar(im, cax=cb_sbpl, orientation=orientation)   
    
    if None in [min_tick, max_tick, step_tick]:
        cbar.set_ticks([0, 0.5, 1])
        cbar.set_ticklabels(['low', 'medium', 'high'])
    else:
        cbar.set_ticks(np.arange(min_tick, max_tick*1.000001, step_tick))       
        
    cbar.ax.tick_params(labelsize=fontsize, width=linewidth, length = ticklength)
    
    for axis in ['top', 'bottom', 'left', 'right']:
        cbar.ax.spines[axis].set_linewidth(linewidth)