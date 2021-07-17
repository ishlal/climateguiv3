#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[1]:


import pandas as pd
import numpy as np
#import seaborn as sns
import math
#import networkx as nx

#from pymcdm import methods as mcdm_methods
#from pymcdm import weights as mcdm_weights
#from pymcdm import normalizations as norm
#from pymcdm import correlations as corr
#from pymcdm.helpers import rankdata, rrankdata

#%matplotlib inline
import matplotlib.pyplot as plt

from sklearn import preprocessing

# reading in data and pre-processing
data_sheet = 'SWAGv3.xlsx'

acsdf_orig = pd.read_excel('SWAGv3.xlsx', header = 1, index_col = 0, nrows = 12)
#acsdf_orig.head()
acsdf_orig = acsdf_orig.drop(axis=1, columns=['Unnamed: 8'])
COLNAMES = acsdf_orig.columns
#acsdf_orig



# normalize data. Specify the max score that should be given:
cols = acsdf_orig.columns
normalize_max = 100
x = acsdf_orig.values
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)*normalize_max
acsdf = pd.DataFrame(x_scaled)
# reinstate the column names
acsdf.columns = cols
# create a copy of acsdf
acsdf2 = acsdf.copy(deep = True)

# Read in the weights data. May have to change parameters of pd.read_excel
# depending on the shape of the ACS table.
weights=pd.read_excel('SWAGv3.xlsx',header=1,index_col=0,nrows=1,skiprows=18)
# Drop the extra column
weights=weights.drop(axis=1,columns=['Unnamed: 8'])
# Set the column names
weights.columns = acsdf.columns
weights2 = weights.copy(deep = True)


types = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
acsdf_np = acsdf2.to_numpy()
acsdf_np_f = acsdf_np.astype(float)
weights_np = weights.to_numpy()[0]


# In[2]:


num_criteria = len(acsdf.columns)
num_alternatives = (acsdf.size)//(num_criteria)

min_val = 0
max_val = normalize_max


barplot_x = []
for i in range(num_alternatives):
    barplot_x.append(i)
    

def get_deviation(alt1, alt2):
    # obtain the two rows of the dataframe corresponding to parameters
    alternative1 = acsdf.iloc[alt1]
    alternative2 = acsdf.iloc[alt2]
    deviations = []
    #loop through every entry of each row and append their deviation to a list.
    for i in range(len(alternative1)):
        deviations.append(alternative1[i] - alternative2[i])
    return deviations

def get_vshape_value(deviation, p_param):
    if(deviation<=0):
        return 0
    elif (deviation >= 0 and deviation <= p_param):
        return (deviation/p_param)
    else:
        return 1
    
def get_gaussian_value(deviation, s_param):
    # By definition, if deviation <=0, then the preference value will be 0 because the second
    # alternative is clearly better
    if (deviation <= 0):
        return 0
    else:
        # computes the gaussian value using the exponential equation from above.
        exponent = - ((math.pow(deviation, 2))/(2*math.pow(s_param, 2)))
        return (1 - math.exp(exponent))
    
def get_adj_list(p_q, p_p, p_s):
    adj_list = []
    for i in range (num_criteria):
        preferences = []
        # for each criterion, loop through every pair of alternatives
        for j in range(num_alternatives):
            for k in range(num_alternatives):
                # appen the preference value to a sub-list by getting the gaussian value
                # the preference value is rounded to 3 decimal places
                
                
                #preferences.append(round(get_gaussian_value(get_deviation(j, k)[i], p_s[i]),3))
                preferences.append(round(get_vshape_value(get_deviation(j,k)[i], p_p[i]), 3))
                
        #append the sublist to the larger adjaceny list
        adj_list.append(preferences)
    return adj_list


def get_adj_list_v2():
    adj_list = []
    for i in range (num_criteria):
        preferences = []
        # for each criterion, loop through every pair of alternatives
        if i in ord_cols:
            for j in range(num_alternatives):
                for k in range(num_alternatives):
                    preferences.append(round(get_gaussian_value(get_deviation(j,k)[i], param_s_ord),3))
        else:
            for j in range(num_alternatives):
                for k in range(num_alternatives):
                    # appen the preference value to a sub-list by getting the gaussian value
                    # the preference value is rounded to 3 decimal places
                    preferences.append(round(get_gaussian_value(get_deviation(j, k)[i], param_s),3))
            #append the sublist to the larger adjaceny list
        adj_list.append(preferences)
    return adj_list

def get_adj_list_v3(p_q, p_p, p_s, stress):
    adj_list = []
    for i in range (num_criteria):
        preferences = []
        for j in range(num_alternatives):
            for k in range(num_alternatives):
                if stress[i] == 1:
                    preferences.append(round(get_vshape_value(get_deviation(j,k)[i], p_p[i]),3))
                else:
                    preferences.append(round(get_vshape_value_min(get_deviation(j,k)[i], p_p[i]),3))
        adj_list.append(preferences)
    return adj_list


def get_pi(adj_list, weights):
    pi_sum_list = []
    # loop through every pair of alternatives
    for a in range(num_alternatives):
        for b in range(num_alternatives):
            pi_sum = 0
            # next, obtain the aggregate preference index
            for j in range(num_criteria):
                # Pj(a,b) * wj
                temp_sum = adj_list[j][num_alternatives*a + b] * weights[j]
                pi_sum += temp_sum
            # append this aggregate index to a list (rounded to 3 decimal places)
            pi_sum_list.append(round(pi_sum, 3))
    return pi_sum_list

def get_pos_flow(list_pi, a):
    flow_sum = 0
    for x in range(num_alternatives):
        # flow_sum(a) = sum(pi(a, x))
        flow_sum += list_pi[num_alternatives*a + x]
    # phi_pos(a) = 1/(n-1) * flow_sum(a)
    return (1/(num_alternatives - 1))*(flow_sum)

def get_neg_flow(list_pi, a):
    neg_flow_sum = 0
    for x in range(num_alternatives):
        # neg_flow_sum(a) = sum(pi(x, ))
        neg_flow_sum += list_pi[num_alternatives*x + a]
    # phi_neg(a) = 1/(n-1) * neg_flow_sum(a)
    return (1/(num_alternatives -1))*(neg_flow_sum)

def get_all_pos_flows(pl):
    all_pos_flows = []
    for a in range(num_alternatives):
        # loops through all alternatives and appends their pos_flow to a list
        all_pos_flows.append(get_pos_flow(pl, a))
    return all_pos_flows

def get_all_neg_flows(pl):
    all_neg_flows = []
    for a in range(num_alternatives):
        # loops through all alternatives and appends their neg_flow to a list
        all_neg_flows.append(get_neg_flow(pl, a))
    return all_neg_flows

def get_net_flow(pl):
    pos_flows = get_all_pos_flows(pl)
    neg_flows = get_all_neg_flows(pl)
    net_flows = []
    for a in range(num_alternatives):
        #net flow(a) = pos_flow(a) - neg_flow(a)
        # append net flow for each alternative to a list
        net_flows.append(pos_flows[a] - neg_flows[a])
    return net_flows


def specialized_visual(pos_outflows, neg_outflows, x_low_lim, x_lim, y_low_lim, y_lim):
    #plt.figure(1)
    fig, ax = plt.subplots()
    ax.set_xlim(x_low_lim, x_lim)
    ax.set_ylim(y_low_lim, y_lim)
    for i in np.arange(x_low_lim, x_lim+0.01, 0.1):
        ax.plot([0, 0.01], [i,i], 'k')
    for i in np.arange(y_low_lim, y_lim+0.01, 0.1):
        ax.plot([i,i], [0, 0.01], 'k')
    #x_coor = [0,x_lim]
    #y_coor = [0,y_lim]
    x_coor = [0,1]
    y_coor = [0,1]
    plt.plot(x_coor, y_coor, linewidth = 3)
    
    #ax.text(-0.08, y_low_lim + 0.2, 'Phi-', fontsize='large')
    #ax.text(0.2,y_low_lim-0.08, 'Phi+', fontsize='large')
    plt.xlabel('Phi+')
    plt.ylabel('Phi-')
    ax.axis('on')
    
    for i in range(len(pos_outflows)):
        phi_pos = pos_outflows[i]
        phi_neg = neg_outflows[i]
        plt.scatter(phi_pos, 1-phi_neg, label = f'A{i}')
        ax.text(phi_pos, 1-phi_neg+0.01, f'A{i}')
        start_coors = [phi_pos, phi_pos]
        end_coors = [0, phi_neg]
        #plt.plot(start_coors, end_coors, linewidth = 1, linestyle = '--', color = 'black')
        start_coors = [0, phi_pos]
        end_coors = [phi_neg, phi_neg]
        #plt.plot(start_coors, end_coors, linewidth = 1, linestyle = '--', color = 'black')
    plt.legend()


# In[17]:


#GUI code
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import *

stress = [1, 1, 1, 1, 1, 1, 1]
q_list_gui = [1, 1, 1, 1, 1, 1, 1]
p_list_gui = [6, 6, 6, 6, 6, 6, 6]
s_list_gui = [2, 2, 2, 2, 2, 2, 2]


class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self,master)
        self.createWidgets()

    def createWidgets(self):
        fig=plt.figure(figsize=(8,8))
        ax=fig.add_axes([0.1, 0.1, 0.8, 0.8])
        #ax = fig.add_axes([0, 0, 1, 1])
        

        canvas=FigureCanvasTkAgg(fig,master=root)
        canvas.get_tk_widget().grid(row=1,column=3, rowspan = 20)
        canvas.draw()

        self.plotbutton=tk.Button(master=root, text="plot", command=lambda: self.plot(canvas,ax))
        self.plotbutton.grid(row=15,column=1)
        
        self.labelw = Label(root, text = 'Weights', padx = 10, pady = 5)
        self.labelw.grid(row = 0, column = 0)
        self.labelq = Label(root, text = 'Threshold of Indifference', padx = 10, pady = 5)
        self.labelq.grid(row = 0, column = 1)
        self.labelp = Label(root, text = 'Threshold of Strict Preference', padx = 10, pady = 5)
        self.labelp.grid(row = 0, column = 2)
        
        self.label1 = Label(root, text = 'Budget Expense Cost', padx = 0, pady = 0)
        self.label1.grid(row = 2, column = 0)
        self.horizontal1 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal1.grid(row = 1, column = 0)
        self.label2 = Label(root, text = 'Social Expense to Stakeholders', padx = 0, pady = 0)
        self.label2.grid(row = 4, column = 0)
        self.horizontal2 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal2.grid(row = 3, column = 0)
        self.label3 = Label(root, text = 'Greenhouse Gas Reductions', padx = 0, pady = 0)
        self.label3.grid(row = 6, column=0)
        self.horizontal3 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal3.grid(row = 5, column = 0)
        self.label4 = Label(root, text = 'Health and Safety Effects', padx = 0, pady = 0)
        self.label4.grid(row= 8, column = 0)
        self.horizontal4 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal4.grid(row = 7, column = 0)
        self.label5 = Label(root, text= 'Co-benefits', padx = 0, pady = 0)
        self.label5.grid(row = 10, column = 0)
        self.horizontal5 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal5.grid(row = 9, column = 0)
        self.label6 = Label(root, text = 'Disbenefits', padx = 0, pady =0)
        self.label6.grid(row = 12, column = 0)
        self.horizontal6 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal6.grid(row = 11, column = 0)
        self.label7 = Label(root, text = 'Transition Effects', padx = 0, pady = 0)
        self.label7.grid(row = 14, column = 0)
        self.horizontal7 = Scale(root, from_=0, to = 200, orient = HORIZONTAL)
        self.horizontal7.grid(row = 13, column = 0)

        self.horizontal11 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal11.grid(row = 1, column = 1)
        self.horizontal12 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal12.grid(row = 3, column = 1)
        self.horizontal13 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal13.grid(row = 5, column = 1)
        self.horizontal14 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal14.grid(row = 7, column = 1)
        self.horizontal15 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal15.grid(row = 9, column = 1)
        self.horizontal16 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal16.grid(row = 11, column = 1)
        self.horizontal17 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal17.grid(row = 13, column = 1)
        
        self.horizontal21 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal21.grid(row = 1, column = 2)
        self.horizontal22 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal22.grid(row = 3, column = 2)
        self.horizontal23 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal23.grid(row = 5, column = 2)
        self.horizontal24 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal24.grid(row = 7, column = 2)
        self.horizontal25 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal25.grid(row = 9, column = 2)
        self.horizontal26 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal26.grid(row = 11, column = 2)
        self.horizontal27 = Scale(root, from_=0, to = normalize_max, orient = HORIZONTAL)
        self.horizontal27.grid(row = 13, column = 2)

    def plot(self,canvas,ax):
        #c = ['r','b','g']  # plot marker colors
        ax.clear()         # clear axes from previous plot
        #for i in range(3):
            #theta = np.random.uniform(0,360,10)
            #r = np.random.uniform(0,1,10)
            #ax.plot(theta,r,linestyle="None",marker='o', color=c[i])
            #canvas.draw()
            
            
            
        arr = [self.horizontal1.get(), self.horizontal2.get(), self.horizontal3.get(), self.horizontal4.get(), self.horizontal5.get(), 
          self.horizontal6.get(), self.horizontal7.get()]
        arr = [float(i)/sum(arr) for i in arr]
        arr = np.array(arr)
        
        q_list_gui = [self.horizontal11.get(), self.horizontal12.get(), self.horizontal13.get(), self.horizontal14.get(), 
                      self.horizontal15.get(), self.horizontal16.get(), self.horizontal17.get()]
        p_list_gui = [self.horizontal21.get(), self.horizontal22.get(), self.horizontal23.get(), self.horizontal24.get(), 
                      self.horizontal25.get(), self.horizontal26.get(), self.horizontal27.get()]
        for i in range(len(q_list_gui)):
            s_list_gui[i] = (q_list_gui[i] + p_list_gui[i])/2
        
        adjlistgui = get_adj_list_v3(p_list_gui, q_list_gui, s_list_gui, stress)
        pi_list_gui = get_pi(adjlistgui, arr)
        net_flows_gui = get_net_flow(pi_list_gui)
        ax.bar(barplot_x, net_flows_gui)
        ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        ax.set_xticklabels(['LECC', 'CSCom', 'GreF', 'EBus', 'InCP', 'GasB', 'EVCharge', 'ELawn', 'FoodR', 'BnP', 'WorkR', 'GreElec'])
        plt.xticks(rotation = 45)
        canvas.draw()

root=tk.Tk()
app=Application(master=root)
app.mainloop()


# In[ ]:




