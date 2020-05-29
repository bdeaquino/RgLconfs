#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 12:55:48 2019

@author: belisa
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pygraphviz as pgv

class graph(object):
    def __init__(self, fin_occ, fin_trans, hubs=None, ffac=0.05, facbothdir=0.8,
                 dpi=100, cmap='nipy_spectral_r', 
                 directed=True, nodesep=0.5, splines="false", 
                 outputorder="edgesfirst", overlap="scale", shift_color=0.8,
                 min_size_nodes=0, fac_size_nodes=0.5, nodeshape='circle',
                 nodecolor='black', nodefillcolor='black', nodestyle='filled', 
                 nodelabel='', nodefontsize='10', hubnodeshape='doublecircle',
                 hubnodecolor='black', hubnodefillcolor='black',
                 hubnodestyle='filled', hubnodelabel='', hubnodefontsize='10',
                 node_pos_hfac=1, node_pos_wfac=1,
                 min_size_edges=2.5, fac_size_edges=5,
                 style_sym_edges=None, style_asym_edges=["dashed", "dotted"]):
        
        super(graph, self).__init__()
        
        self.fin_occ=fin_occ
        self.fin_trans=fin_trans
        self.hubs=hubs
        self.ffac=ffac
        self.facbothdir=facbothdir
        self.cmap=plt.get_cmap(cmap)
        self.shift_color=shift_color
        self.min_size_nodes=min_size_nodes
        self.fac_size_nodes=fac_size_nodes
        
        self.nodeshape=nodeshape
        self.nodecolor=nodecolor
        self.nodefillcolor=nodefillcolor
        self.nodestyle=nodestyle
        self.nodelabel=nodelabel
        self.nodefontsize=nodefontsize
        if hubs:
            self.hubnodeshape=hubnodeshape
            self.hubnodecolor=hubnodecolor
            self.hubnodefillcolor=hubnodefillcolor
            self.hubnodestyle=hubnodestyle
            self.hubnodelabel=hubnodelabel
            self.hubnodefontsize=hubnodefontsize
            
        self.node_pos_hfac=node_pos_hfac
        self.node_pos_wfac=node_pos_wfac
        
        self.min_size_edges=min_size_edges
        self.fac_size_edges=fac_size_edges
        self.style_sym_edges=style_sym_edges
        self.style_asym_edges=style_asym_edges       
        
        #initializing graph
        self.G = pgv.AGraph(directed=directed, nodesep=nodesep, 
                       splines=splines, outputorder=outputorder, 
                       overlap=overlap, dpi=dpi)
        
        self.h = self.occupancy_matrix_from_file()
        self.nodes, self.edges, self.source, self.target, self.f, self.fmin, self.fmax, self.flen = self.sep_transitions_data_from_file()        
        
    def occupancy_matrix_from_file(self):
        h = np.loadtxt(self.fin_occ)
    
        maxh = 0
        for j in h:
            maxh = max(maxh, max(j))
            
        return h.T/maxh
    
    def sep_transitions_data_from_file(self):
        source = []
        target = []
        f = []
        with open(self.fin_trans) as fp:
            for j, line in enumerate(fp):  
                a,b,c,d,e = line.split()
                
                source.append((int(a),int(b)))
                target.append((int(c),int(d)))
                f.append(float(e))
                
        f = np.array(f)/sum(f)
        fmax = max(f)
        f /= fmax
        fmin = min(f)
        fmax2 = max(f)
        flen = fmax2-fmin
        
        idx_good = [k for k,x in enumerate(f) if x > self.ffac*fmax2]
        f = f[idx_good]
        
        source = [x for k,x in enumerate(source) if k in idx_good]
        target = [x for k,x in enumerate(target) if k in idx_good]
        nodes = list(set(source+target))
        edges = list(zip(source,target))
        
        return nodes, edges, source, target, f, fmin, fmax, flen
    
    def dir_edges(self):
        edges_dec = []
        edges_inc = []
        
        f_dec = []
        f_inc = []
        
        for ((a,b),(c,d),e) in zip(self.source,self.target,self.f):
            if c > a:
                direc = 'inc'
            elif c < a:
                direc = 'dec'
            elif c == a:
                if d > b:
                    direc = 'inc'
                elif d < b:
                    direc = 'dec'
                    
            if direc == 'dec':
                edges_dec.append(((a,b),(c,d)))
                f_dec.append(e)
            elif direc == 'inc':
                edges_inc.append(((a,b),(c,d)))
                f_inc.append(e)
                
        return edges_dec, f_dec, edges_inc, f_inc
    
    def add_nodes(self):
        for (a,b) in self.nodes:
            pos=str(self.node_pos_wfac*a)+','+str(self.node_pos_hfac*b)+'!'
            color, width = color_n_width(self.h[a-1][b-1], self.cmap, self.shift_color, 
                                              self.min_size_nodes, 
                                              self.fac_size_nodes, 
                                              self.fmin, self.flen)
            
            if self.hubs and (a,b) in self.hubs:
                print(a,b)
                self.G.add_node((int(a),int(b)), width=width, shape=self.hubnodeshape,
                           color=self.hubnodecolor, fillcolor=self.hubnodefillcolor, 
                           fixedsize='true', fontsize=self.hubnodefontsize, 
                           style=self.hubnodestyle, label=self.hubnodelabel, pos=pos)
            else:
                self.G.add_node((int(a),int(b)), width=width, shape=self.nodeshape,
                           color=self.nodecolor, fillcolor=self.nodefillcolor, 
                           fixedsize='true', fontsize=self.nodefontsize, 
                           style=self.nodestyle, label=self.nodelabel, pos=pos)
            
    def add_edges(self):
        edges = []
        
        edges_dec, f_dec, edges_inc, f_inc = self.dir_edges()  
        
        fac1 = self.min_size_edges
        fac2 = self.fac_size_edges
        j = 0
        for (a,b) in edges_dec:
            color1, ewidth1 = color_n_width(f_dec[j], self.cmap, self.shift_color, fac1, fac2, 
                                            self.fmin, self.flen)
            if (b,a) in edges_inc:
                k = edges_inc.index((b,a))
                fmx = max(f_dec[j],f_inc[k])
                fmn = min(f_dec[j],f_inc[k])
                
                color, ewidth = color_n_width(fmx, self.cmap, self.shift_color, fac1, fac2, 
                                              self.fmin, self.flen)
                color2, ewidth2 = color_n_width(f_inc[k], self.cmap, self.shift_color, fac1, fac2, 
                                                self.fmin, self.flen)
                   
                if fmx and fmn:
                    fmx = max(f_dec[j],f_inc[k])
                    fmn = min(f_dec[j],f_inc[k])
                              
                    if fmn/fmx > self.facbothdir:
                        if self.style_sym_edges:
                            self.G.add_edge((a,b), color=color, penwidth=ewidth, 
                                       style=self.style_sym_edges, dir="both")
                        else:
                            self.G.add_edge((a,b), color=color, 
                                            penwidth=ewidth, dir="both")
                        edges.append((a,b))
                        edges.append((b,a))
                    else:      
                        if f_dec[j] > f_inc[k]:
                            style1 = self.style_asym_edges[0]
                            style2 = self.style_asym_edges[1]
                        else:
                            style1 = self.style_asym_edges[1]
                            style2 = self.style_asym_edges[0]
                        self.G.add_edge(a,b, color=color1, 
                                        penwidth=ewidth1, style=style1)
                        self.G.add_edge(b,a, color=color2, 
                                        penwidth=ewidth2, style=style2)
                        edges.append((a,b))
                        edges.append((b,a))
                        
                elif fmn is None:
                    if f_dec[j]:
                        self.G.add_edge((a,b), color=color1, penwidth=ewidth1, 
                                   style=self.style_asym_edges[0])
                        edges.append((a,b))
                    elif f_inc[j]:
                        self.G.add_edge((b,a), color=color2, penwidth=ewidth2, 
                                   style=self.style_asym_edges[0])                                  
                        edges.append((b,a))
            else:
                if f_dec[j]:
                    self.G.add_edge((a,b), color=color1, penwidth=ewidth1, 
                               style=self.style_asym_edges[0])
                    edges.append((a,b))
            j+=1
            
        j = 0
        for (a,b) in edges_inc:
            if (a,b) not in edges:
                if f_inc[j]:
                    color2, ewidth2 = color_n_width(f_inc[j], self.cmap, self.shift_color, fac1, fac2, 
                                                    self.fmin, self.flen)
                    self.G.add_edge((a,b), color=color2, penwidth=ewidth2, 
                               style=self.style_asym_edges[0])
                
            j+=1
    
    def draw_graph(self, graphname, fmt='eps'):
        try:
            self.add_nodes()
            self.add_edges()
            
            print('Drawing '+graphname)
            self.G.draw(graphname, prog='neato', format=fmt)
            
        except NameError:
            print("Error: it was not possible to draw the network.")
            print("Please try choosing a different value for ffac or facbothdir.")
            print("Alternatively, set spline='false'.")

###############################################################################

def color_n_width(a, cmap, m, fac1, fac2, vmin, vlen):
        color = cmap((a-m*vmin)/vlen)
        color = mpl.colors.rgb2hex(color)
        width = fac1+fac2*a
        
        return color, width
