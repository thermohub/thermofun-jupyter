import thermofun.PyThermoFun as fun

import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt
import re 
import os
import base64
import difflib
import random
import requests
import ipywidgets as widgets

from collections import OrderedDict
from itertools import cycle
from matplotlib.lines import Line2D
from ipywidgets import Layout, Label
from ipywidgets import interact, interact_manual

from IPython.display import HTML
from IPython.display import display, clear_output
from beakerx import TableDisplay

from common.functions import create_csv_download_link


def plot_substances_properties_vs_temperature(results_csv_file, substances, pressure):

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
    legend1 = plt.legend(by_label.values(), by_label.keys(), loc=9, bbox_to_anchor=(1.15, 1)) #, loc=legend_loc)
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

def multi_checkbox_widget(records, properties):
    """ Widget with a search field and lots of checkboxes """
    search_widget = widgets.Text(description="Search")
    options_dict = {description: widgets.Checkbox(description=description, value=False) for description in records}
   # options_dict.values()[0].value = True
    props_dict = {description: widgets.Checkbox(description=description, value=False) for description in properties}
    options = [options_dict[description] for description in records]
    options_props = [props_dict[description] for description in properties]
    box_layout = Layout(overflow='scroll',
                    #border='3px solid black',
                    width='max-content',
                    height='300px',
                    flex_flow='column',
                    display = 'flex',
                    align_items = 'stretch',
                    color='blue')
    txt = widgets.Label(value='Properties')
    options_widget = widgets.VBox(options, layout=box_layout)
    props_widget = widgets.VBox(options_props, layout=box_layout)
    props_txt = widgets.VBox([txt, props_widget])
    multi_select = widgets.VBox([search_widget, options_widget])
    multi_data = widgets.HBox([multi_select, props_txt])

    # Wire the search field to the checkboxes
    def on_text_change(change):
        search_input = change['new']
        if search_input == '':
            # Reset search field
            new_options = [options_dict[description] for description in records]
        else:
            # Filter by search field using difflib.
  #          close_matches = difflib.get_close_matches(search_input, records, n=5, cutoff=0.0)
            # close_matches = list(filter(lambda x: search_input in x, records)) 
            close_matches = [x for x in records if re.search(search_input, x, re.IGNORECASE)]
            new_options = [options_dict[description] for description in close_matches]
        options_widget.children = new_options

    search_widget.observe(on_text_change, names='value')
    
    return multi_data

def make_tabs(select_subst, select_reactions, names):
    tab_contents = names
    children = [select_subst, select_reactions]
    tab = widgets.Tab()
    tab.children = children
    for i in range(len(children)):
        tab.set_title(i, tab_contents[i])
    tab
    return tab

def ui(databases_w, dropdown_w):
    table_style = {'description_width': 'initial'}
    table_layout = {'width':'150px', 'min_width':'150px', 'height':'28px', 'min_height':'28px'}
    row_layout = {'width':'122px', 'min_width':'122px'}
    row_layout2 = {'width':'250px', 'min_width':'250px'}

    databases_w.layout = row_layout
    dropdown_w.layout = row_layout2

    table_header_0_widget = widgets.Text(
        description='',
        disabled=True,
        button_style='',
        tooltip='',
        icon='',
        layout=row_layout,
        style=table_style
    )

    t0_1_widget = widgets.BoundedIntText( description='Tmin',
                                        value=25,
                                        min=0,
                                        max=1000,
                                        step=1,
                                        disabled=True,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    t0_2_widget = widgets.BoundedIntText( description='Tmax',
                                        value=150,
                                        min=0,
                                        max=1000,
                                        step=1,
                                        disabled=True,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   ) 
    t0_3_widget = widgets.BoundedIntText( description='Tstep',
                                        value=5,
                                        min=0,
                                        max=1000,
                                        step=1,
                                        disabled=True,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    t1_1_widget = widgets.BoundedIntText( description='Pmin',
                                        value=0,
                                        min=0,
                                        max=100000,
                                        step=1,
                                        disabled=True,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    t1_2_widget = widgets.BoundedIntText( description='Pmax',
                                        value=0,
                                        min=0,
                                        max=100000,
                                        step=1,
                                        disabled=True,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   ) 
    t1_3_widget = widgets.BoundedIntText( description='Pstep',
                                        value=0,
                                        min=0,
                                        max=100000,
                                        step=1,
                                        disabled=True,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    

    hbox1 = widgets.HBox([databases_w, dropdown_w, t0_1_widget, t0_2_widget, t0_3_widget])
    hbox2 = widgets.HBox([Label(value='',layout=row_layout), Label(value='', layout=row_layout2), t1_1_widget, t1_2_widget, t1_3_widget])
    ui = widgets.VBox([hbox1, hbox2])

    display(ui)

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

db_file = 'databases/aq17-fun.json'

def load_widgets(dfDatabase, dfSelectSubst, dfSelectReact) :
    dirName = 'databases'
    origin  = widgets.Dropdown(
    options = getListOfFiles(dirName),
    value   = db_file,
  #  description='Available Databases:',
    layout  = Layout(width= 'max-content')
    )
#    database_file = 'databases/aq17-fun.json'
#    global database
#    database = fun.Database(database_file)
#    select_subst = multi_checkbox_widget(database.mapSubstances())
#    select_react = multi_checkbox_widget(database.mapReactions())
#    global db_file
    
    tabs = make_tabs(dfSelectSubst[db_file], dfSelectReact[db_file], ['Substances', 'Reactions'])
    
    plot_out=widgets.Output()
    link_out=widgets.Output()

    box_layout2 = Layout(overflow='scroll',
                    #border='3px solid black',
                    width='max-content',
                    height='500px',
    #                flex_direction='row',
    #                display='flex',
    #                align_items = 'stretch',
                    color='blue')
    
    table_out=widgets.Output()#layout=box_layout2) 

    button = widgets.Button(
    value=False,
    description='Calculate',
    disabled=False,
    button_style='success', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Calculate Properties',
#    icon='check'
    )

    def on_dropdown_change(change):
        if change['name'] == 'value' and (change['new'] != change['old']):
            global db_file 
            db_file = change['new']
#            global database
#            database = fun.Database(database_file)
#            select_subst_n = multi_checkbox_widget(dfDatabase[f].mapSubstances())
#            select_react_n = multi_checkbox_widget(database.mapReactions())
#            select_subst.children = dfSelectSubst[f].children
#            select_react.children = dfSelectReact[f].children
            tabs_new = make_tabs(dfSelectSubst[db_file], dfSelectReact[db_file], ['Substances', 'Reactions'])
            tabs.children = tabs_new.children

            with plot_out:
                clear_output()
            with table_out:
                clear_output()
    
    origin.observe(on_dropdown_change)
    
    items = [Label(value='Available Databases:'), origin]

    box_layout = Layout(
                    width='max-content',
                    height='',
                    flex_direction='row',
                    display='flex')
    carousel = widgets.HBox(children=items, layout=box_layout)
    

    def on_button_clicked(b):
  #      with output:
  #      print("Button clicked.")
        global db_file
        #display(select_subst)
        batch = fun.ThermoBatch(dfDatabase[db_file])
        
        op = fun.BatchPreferences()
        op.isFixed = True # values are written using fixed-point notation
        # if True properties of reactions are calculated from the properties of reactants
        # if False properties of reactions are calculated from the method in the reaction record
        op.calcReactFromSubst   = False 
        # if True properties of substances are calculated from the dependent reactions
        # if False properties of substances are calculated from the method in the substance record
        op.calcSubstFromReact   = False 
        batch.setBatchPreferences(op)
        
        batch.setPropertiesUnits(["temperature", "pressure"],["degC","bar"])
        batch.setPropertiesDigits(["gibbs_energy","entropy", "volume",
                            "enthalpy","logKr", "temperature", "pressure"], [0, 4, 4, 4, 4, 0, 0])
        #properties = ["gibbs_energy", "enthalpy", "entropy"]
        
        temperature_pressure_pairs = [[50,1000],  [150,1000], [200,1000], [250,1000], [300,1000], [350,1000], 
                              [400,1000], [450,1000], [500,1000], [550,1000], [600,1000], [650,1000], 
                              [700,1000], [800,1000], [900,1000], [1000,1000]]
        if tabs.selected_index == 0:
            selected_options = [w.description for w in dfSelectSubst[db_file].children[0].children[1].children if w.value]
            selected_properties = [w.description for w in dfSelectSubst[db_file].children[1].children[1].children if w.value]
            batch.thermoPropertiesSubstance(temperature_pressure_pairs, selected_options, selected_properties).toCSV("results.csv")
        else:
            selected_options = [w.description for w in dfSelectReact[db_file].children[0].children[1].children if w.value]
            selected_properties = [w.description for w in dfSelectReact[db_file].children[1].children[1].children if w.value]
            batch.thermoPropertiesReaction(temperature_pressure_pairs, selected_options, selected_properties).toCSV("results.csv")
        
#        ax=plt.gca()
        with link_out:
            clear_output(wait=True)
            display(create_csv_download_link("results.csv"))
        
        with plot_out:
            clear_output(wait=True)
            plot_substances_properties_vs_temperature('results.csv', selected_options, 1000)
            plt.show()
        with table_out:
            clear_output(wait=True)
            df = pd.read_csv('results.csv')
            #display(df)
            table = TableDisplay(df)
            display(table)
            
    button.on_click(on_button_clicked)

    #display(carousel)
    ui(Label(value='Available Databases:'), origin)
    display(tabs)
    #display(button)
    hbox=widgets.HBox(children=(button, link_out))
    display(hbox)
    
    
    #display(link_out)
    tabs_results = make_tabs(table_out, plot_out, ['Table', 'Plot'])
 #   hbox=widgets.HBox(children=(table_out, plot_out))
    display(tabs_results)

def initialize_widgets() :
    dfDatabase    = {}
    dfSelectSubst = {}
    dfSelectReact = {}

    data_out=widgets.Output()
    
    button = widgets.Button(
    value=False,
    description='Load Data',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Load Data',
#    icon='check'
    )

    files = getListOfFiles('databases')
#    print('files')

    progress = widgets.IntProgress(
    value=0,
    min=0,
    max=len(files),
 #   step=1,
 #   description='Loading:',
    bar_style='', # 'success', 'info', 'warning', 'danger' or ''
    orientation='horizontal'
    )

    hbox=widgets.HBox(children=(button, progress))
    display(hbox)

    def on_button_clicked(b):
        for i, f in enumerate(files):
            progress.value = i+1
            dfDatabase[f] = fun.Database(f)
            dfSelectSubst[f] = multi_checkbox_widget(dfDatabase[f].mapSubstances(), ["gibbs_energy","enthalpy","entropy","heat_capacity_cp","heat_capacity_cv","volume","helmholtz_energy","internal_energy"])
            dfSelectReact[f] = multi_checkbox_widget(dfDatabase[f].mapReactions(), ["reaction_gibbs_energy","reaction_helmholtz_energy","reaction_internal_energy","reaction_enthalpy","reaction_entropy","reaction_volume","reaction_heat_capacity_cp","reaction_heat_capacity_cv","logKr","lnKr"])
            progress.value = i+1
            #time.sleep(0.1)
        with data_out:
            clear_output(wait=True)
            load_widgets(dfDatabase, dfSelectSubst, dfSelectReact)
            
    button.on_click(on_button_clicked)

    display(data_out)

    


