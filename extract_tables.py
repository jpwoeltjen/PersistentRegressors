#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:11:37 2019

@author: jan
"""

import pandas as pd
df = pd.read_csv("tabula-efficienttests_appendix032005.csv", header=None, index_col=0)

directory = 'ci'

df_gls95 = df.iloc[:,[0,1]]
df_gls95.columns = ['cl','cu']
df_gls95.to_csv(f'{directory}/df_gls95.csv',index=True)


df_gls90 = df.iloc[:,[2,3]]
df_gls90.columns = ['cl','cu']
df_gls90.to_csv(f'{directory}/df_gls90.csv',index=True)


df_gls80 = df.iloc[:,[4,5]]
df_gls80.columns = ['cl','cu']
df_gls80.to_csv(f'{directory}/df_gls80.csv',index=True)

df = pd.read_csv("ci_rho_q_test.csv", header=None, index_col=0)


n_rows = 61
n_tables = 10
delta_start = -1
delta_diff = 0.025
for i in range(n_tables):
    
    table = df.iloc[i*n_rows:(i+1)*n_rows]
    if i==0:
        columns = [-0.999, delta_start+delta_diff*(4*i+1), delta_start+delta_diff*(4*i+2),
                   delta_start+delta_diff*(4*i+3)]
    else:
        columns = [round(delta_start+delta_diff*(4*i),5),
                   round( delta_start+delta_diff*(4*i+1),5),
                   round(delta_start+delta_diff*(4*i+2),5),
                   round(delta_start+delta_diff*(4*i+3),5)]
    
    print(columns)
    
    for j in range(len(columns)):
        df_gls_delta = table.iloc[:,j*2:(j+1)*2]
        df_gls_delta.columns = ['cl','cu']
        df_gls_delta.to_csv(f'{directory}/df_gls_delta{columns[j]}.csv',index=True)
    
    
