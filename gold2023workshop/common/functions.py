import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt
import re 
import os

import difflib
import random
import requests
import ipywidgets as widgets

from collections import OrderedDict
from itertools import cycle
from collections import OrderedDict
from matplotlib.lines import Line2D
from ipywidgets import Layout, Label

from IPython.display import HTML
from IPython.display import display
import base64

def create_csv_download_link(filename):
    df = pd.read_csv(filename)
    title = "Download the " + filename + " file"
    csv = df.to_csv()
    b64 = base64.b64encode(csv.encode())
    payload = b64.decode()
    html = '<a download="{filename}" style="font-size: 20px; text-decoration: none" href="data:text/csv;base64,{payload}" target="_blank">{title}</a>'
    html = html.format(payload=payload,title=title,filename=filename)
    return HTML(html)


def plot_substances_properties_vs_temperature(results_csv_file, substances, pressure=0):

    # plot settings
    mpl.rcParams['lines.linewidth']=2
    mpl.rcParams['axes.labelsize']=20
    mpl.rcParams['axes.linewidth']=2
    mpl.rcParams['font.size']=18
    mpl.rcParams['figure.figsize']=[9,9]

#    fig = plt.figure()


    plt.rc('grid', linestyle="--", color='gray')
    plt.grid(True)
    markers = [m for m in Line2D.markers]
    data = pd.read_csv(results_csv_file)
    if pressure != 0.0:
    	data = data.loc[data.iloc[:, 2] == pressure]
    c_cycle = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    for s, substance in enumerate(substances):
        data_select = data.loc[data['Symbol']==substance]
        c = next(c_cycle)
        m_cycle = cycle(markers)
        next(m_cycle) 
        next(m_cycle)
        for column in data_select.columns[3:]: # loop over properties
            plt.plot(data_select.iloc[:, 1], data_select[column], color=c, 
                     marker=next(m_cycle), label=substance, markersize=12, markeredgecolor="w")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    legend1 = plt.legend(by_label.values(), by_label.keys()) #, loc=9, bbox_to_anchor=(1.15, 1)) #, loc=legend_loc)
    lines = plt.gca().get_lines()
    
    m_cycle = cycle(markers)
    next(m_cycle) 
    next(m_cycle)
    props_indexes= range(len(data.columns[3:]))
    legend2 = plt.legend([lines[i] for i in props_indexes], data_select.columns[3:],  ncol=1)
    for handle, text in zip(legend2.legendHandles, legend2.get_texts()):
        text.set_color('black')
        #handle.set_marker(next(m_cycle))
        #handle.set_color('black')
        
    plt.gca().add_artist(legend1)
    plt.xlabel("Temperature ($^\circ$C)")
    plt.ylabel("ThermoProp")

    return plt
    
def plot_properties_vs_temperature(results_csv_file, substances, pressure=0):

    # plot settings
    mpl.rcParams['lines.linewidth']=2
    mpl.rcParams['axes.labelsize']=20
    mpl.rcParams['axes.linewidth']=2
    mpl.rcParams['font.size']=18
    mpl.rcParams['figure.figsize']=[9,9]

#    fig = plt.figure()


    plt.rc('grid', linestyle="--", color='gray')
    plt.grid(True)
    markers = [m for m in Line2D.markers]
    data = pd.read_csv(results_csv_file)
    if pressure != 0.0:
    	data = data.loc[data.iloc[:, 2] == pressure]
    c_cycle = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    for s, substance in enumerate(substances):
        data_select = data.loc[data['Symbol']==substance]
        c = next(c_cycle)
        m_cycle = cycle(markers)
        next(m_cycle) 
        next(m_cycle)
        for column in data_select.columns[3:]: # loop over properties
            plt.plot(data_select.iloc[:, 1], data_select[column], color=c, 
                     marker=next(m_cycle), label=column + ' '+ substance, markersize=12, markeredgecolor="w")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    legend1 = plt.legend(by_label.values(), by_label.keys()) #, loc=9, bbox_to_anchor=(1.15, 1)) #, loc=legend_loc)
    lines = plt.gca().get_lines()
        
    plt.gca().add_artist(legend1)
    plt.xlabel("Temperature ($^\circ$C)")
    plt.ylabel("ThermoProp")

    return plt


def multi_checkbox_widget(descriptions):
    """ Widget with a search field and lots of checkboxes """
    search_widget = widgets.Text(description="Search")
    options_dict = {description: widgets.Checkbox(description=description, value=False) for description in descriptions}
    options = [options_dict[description] for description in descriptions]
    box_layout = Layout(overflow='scroll',
                    #border='3px solid black',
                    width='max-content',
                    height='300px',
                    flex_flow='column',
                    display = 'flex',
                    align_items = 'stretch',
                    color='blue')

    options_widget = widgets.VBox(options, layout=box_layout)
    multi_select = widgets.VBox([search_widget, options_widget])

    # Wire the search field to the checkboxes
    def on_text_change(change):
        search_input = change['new']
        if search_input == '':
            # Reset search field
            new_options = [options_dict[description] for description in descriptions]
        else:
            # Filter by search field using difflib.
  #          close_matches = difflib.get_close_matches(search_input, descriptions, n=5, cutoff=0.0)
            # close_matches = list(filter(lambda x: search_input in x, descriptions)) 
            close_matches = [x for x in descriptions if re.search(search_input, x, re.IGNORECASE)]
            new_options = [options_dict[description] for description in close_matches]
        options_widget.children = new_options

    search_widget.observe(on_text_change, names='value')
    return multi_select

def tabs(select_subst, select_reactions):
    tab_contents = ['Substances', 'Reactions']
    children = [select_subst , select_reactions]
    tab = widgets.Tab()
    tab.children = children
    for i in range(len(children)):
        tab.set_title(i, tab_contents[i])
    tab
    return tab

'''
For the given path, get the List of all files in the directory tree 
'''
def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
    # Create full path
        fullPath = os.path.join(dirName, entry)
    # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
    return allFiles

def select_databases(dirName) :
    origin = widgets.Dropdown(
    options= getListOfFiles(dirName),
    value= 'databases/aq17-fun.json',
  #  description='Available Databases:',
    layout= Layout(width= 'max-content')
    )

    def on_dropdown_change(change):
        if change['name'] == 'value' and (change['new'] != change['old']):
            print('do something with the change')
    
    origin.observe(on_dropdown_change)


    items = [Label(value='Available Databases:'), origin]

    box_layout = Layout(
                    width='max-content',
                    height='',
                    flex_direction='row',
                    display='flex')
    carousel = widgets.HBox(children=items, layout=box_layout)

    return display(carousel)
    #return origin

           
