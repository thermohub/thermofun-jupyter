import thermofun as fun

import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt
import re 
import os
import base64
import copy
#import difflib
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
import uuid
from IPython.core.display import display, HTML, Javascript

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import uuid
from IPython.display import display_javascript, display_html, display
import json

import threading
import time

from common.progress2 import log_progress

''' function for rendering JSON record'''
class RenderJSON_(object):
    def __init__(self, json_data):
        if isinstance(json_data, dict):
            self.json_str = json.dumps(json_data)
        else:
            self.json_str = json
        self.uuid = str(uuid.uuid4())

        self._ipython_display_()

    def _ipython_display_(self):
        display_html('<div id="{}" style="height: 600px; width:100%;"></div>'.format(self.uuid),
            raw=True
        )
        display_javascript("""
        require(["https://rawgit.com/caldwell/renderjson/master/renderjson.js"], function() {
          renderjson.set_show_to_level(1)
          document.getElementById('%s').appendChild(renderjson(%s))
        });
        """ % (self.uuid, self.json_str), raw=True)

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
    legend1 = plt.legend(by_label.values(), by_label.keys(), loc=9, bbox_to_anchor=(1.2, 1)) #, loc=legend_loc)
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

def make_options(options_dict, description):
    options_dict.add(description, widgets.Checkbox(description=description, value=False))
    return options_dict

def input_reactions_widget(substances, properties):
    input_reactions_textbox = widgets.Textarea(
    value='# Write reactions using substance \'symbols\'\n# e.g. Calcite (left column) not CaCO3\n# Carefully balance the reactions!\nH2O@ = H+ + OH-',
    description='',
    layout=Layout(max_width='max-content', 
    min_width='360px',
    max_height='max-content', 
    height='340px')
)

    """ Widget with a search field and lots of checkboxes """
    search_widget = widgets.Text(description="Search", style={'description_width': '50px'})

    subst_dict = {description: widgets.Text(layout={'width':'285px'}, style={'description_width': '100px'},
        description=description, value=substances[description].formula(), disabled=True) for description in substances.keys()}
    box_layout = Layout(#overflow_y='scroll',
                    #border='3px solid black',
                    width='max-content',
                    height='320px',
                    flex_flow='column',
                    display = 'flex',
                    align_items = 'stretch',
                    color='blue')
    substs = [subst_dict[description] for description in substances.keys()]#records]

    subst_widget = widgets.VBox(substs, layout=box_layout)
    
    # Wire the search field to the checkboxes
    def on_text_change(change):
        search_input = change['new']
        if search_input == '':
            # Reset search field
            new_options = [subst_dict[description] for description in substances.keys()]
        else:
            # Filter by search field using difflib.
  #          close_matches = difflib.get_close_matches(search_input, records, n=5, cutoff=0.0)
            # close_matches = list(filter(lambda x: search_input in x, records)) 
            close_matches = [x for x in substances.keys() if re.search(search_input, substances[x].formula(), re.IGNORECASE)]
            new_options = [subst_dict[description] for description in close_matches]
        subst_widget.children = new_options

    props_dict = {description: widgets.Checkbox(description=description, value=False, style={'description_width': '10px'}) for description in properties}
    props_dict[properties[0]].value = True
    options_props = [props_dict[description] for description in properties]
    txt = widgets.Label(value='Properties')
    props_widget = widgets.VBox(options_props, layout=box_layout)
    props_txt = widgets.VBox([txt, props_widget])

    search_widget.observe(on_text_change, names="value")
    multi_select = widgets.VBox([search_widget, subst_widget])
    multi_data = widgets.HBox([ multi_select, input_reactions_textbox, props_txt])

    return multi_data

# somehow remember the text area per database, not really important


def multi_checkbox_widget(records, properties, reaction_equations = []):
    """ Widget with a search field and lots of checkboxes """
    search_widget = widgets.Text(description="Search", style={'description_width': '50px'})

    ''' Creates reac_keys_dic to display reaction equations but still be able 
    to call the calculation function using the rection symbol '''
    records_keys=[]
    if len(reaction_equations)>0:
        for i, key_ in  enumerate(records.keys()):
            records_keys.append(reaction_equations[i])
            reac_keys_dic.update( {reaction_equations[i] : key_} )
    else:
        records_keys = records.keys()

    # options_dict = {}
    # for descr in log_progress(iter(records_keys), every=1):
    #     options_dict=make_options(options_dict, descr)

    options_dict = {description: widgets.Checkbox(description=description, value=False,style={'description_width': '10px'}, layout={'width':'285px'}) for description in log_progress(iter(records_keys), every=1, size=len(records_keys), name='Records', progress_out=progress_out)}
   # options_dict.values()[0].value = True
    """ Properties of subst or reactions """
    props_dict = {description: widgets.Checkbox(description=description, value=False, style={'description_width': '10px'}) for description in properties}
    props_dict[properties[0]].value = True
    
    options = [options_dict[description] for description in records_keys]#records]
    options_props = [props_dict[description] for description in properties]
    box_layout = Layout(#overflow_y='scroll',
                    #border='3px solid black',
                    width='max-content',
                    height='320px',
                    flex_flow='column',
                    display = 'flex',
                    align_items = 'stretch',
                    color='blue')

    box_layout2 = Layout(#overflow_y='scroll',
                    #border='3px solid black',
                    width='300px',
                    height='320px',
                    flex_flow='column',
                    display = 'flex',
                    align_items = 'stretch',
                    color='blue')

    button_details_records = widgets.Button(
    value=False,
    description='Show Details of Selected Records',
    disabled=False,
    button_style='info', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Show Details of Selected Records',
    layout  = Layout(width= 'max-content')
#    icon='check'
    )

    def make_accordion(selected_records):
        recs_dic = {rec: widgets.Output() for rec in selected_records}
        recs = [recs_dic[rec] for rec in selected_records]
        accordion = widgets.Accordion(children=recs, selected_index=None)
        for i, r in enumerate(selected_records):
            accordion.set_title(i, r)
            with accordion.children[i]:
                #RenderJSON_(records[r].jsonString())
                clear_output()
                uuid_ = str(uuid.uuid4())
                display(HTML('<div id="{}" style="height: auto; width:100%;"></div>'.format(uuid_)))
                display(HTML("""<script>
                require(["https://rawgit.com/caldwell/renderjson/master/renderjson.js"], function() {
                renderjson.set_show_to_level(1)
                document.getElementById('%s').appendChild(renderjson(%s))
                });</script>
                """ % (uuid_, records[r].jsonString())))

#            with accordion.children[i]:
#                display(display(RenderJSON(records[r].jsonString())))
#                print(i)
        
        return accordion
    
    options_widget = widgets.VBox(options, layout=box_layout)
    ac_ = make_accordion([w.description for w in options_widget.children if w.value])

    txt = widgets.Label(value='Properties')
    props_widget = widgets.VBox(options_props, layout=box_layout)
    props_txt = widgets.VBox([txt, props_widget])
    multi_select = widgets.VBox([search_widget, options_widget])
    records_description = widgets.VBox( [ac_], layout=box_layout)
    records_description_ = widgets.VBox([button_details_records, records_description])
    multi_data = widgets.HBox([multi_select, props_txt, records_description_])

    # Wire the search field to the checkboxes
    def on_text_change(change):
        search_input = change['new']
        if search_input == '':
            # Reset search field
            new_options = [options_dict[description] for description in records_keys]
        else:
            # Filter by search field using difflib.
  #          close_matches = difflib.get_close_matches(search_input, records, n=5, cutoff=0.0)
            # close_matches = list(filter(lambda x: search_input in x, records)) 
            close_matches = [x for x in records_keys if re.search(search_input, x, re.IGNORECASE)]
            new_options = [options_dict[description] for description in close_matches]
        options_widget.children = new_options

    search_widget.observe(on_text_change, names="value")

    def on_button_clicked(b):
        selected_options_ = [w.description for w in options_widget.children if w.value]
        
        reac_keys = []
        global reac_keys_dic
        if len(reaction_equations)>0:
            for k in selected_options_:
                reac_keys.append(reac_keys_dic[k])
            selected_options_= reac_keys

        new_ac =  make_accordion(selected_options_)
        records_description.children = [new_ac]
            
    button_details_records.on_click(on_button_clicked)
    
    return multi_data

def make_tabs_3(select_subst, select_reactions, input_reactions, names):
    tab_contents = names  
    children = [select_subst, select_reactions, input_reactions]
    tab = widgets.Tab(layout=Layout(max_width='max-content', min_width='600px'))
    tab.children = children
    for i in range(len(children)):
        tab.set_title(i, tab_contents[i])
    return tab

def make_tabs(select_subst, select_reactions, names):
    tab_contents = names
    children = [select_subst, select_reactions]
    tab = widgets.Tab(layout=Layout(max_width='max-content', min_width='600px'))
    tab.children = children
    for i in range(len(children)):
        tab.set_title(i, tab_contents[i])
    return tab

''' Global variables '''
temperatures = [25, 150, 5]
pressures    = [0,0,0]

op = fun.BatchPreferences()
op.isFixed = True # values are written using fixed-point notation
reac_keys_dic = {}

db_file = 'databases/0-select-a-database.json'

progress = widgets.IntProgress(min=0, max=100, value=0)
progress.bar_style = 'info'

progress_out=widgets.Output()

style = {'description_width': '40px'}
style2 = {'description_width': '10px'}

toggle_buttons= {
          'SubstFromReact':widgets.Checkbox(description='Substance properties from dependent reaction',  value=False, style=style2),
          'ReactFromSubst':widgets.Checkbox(description='Reaction properties from reactants',  value=True, style=style2)}
        #   'loopTthenP':widgets.Checkbox(description='loop over temperatures then over pressures', value=True),
        #   'loopTPthenRecords':widgets.Checkbox(description='loop over TP pairs followed by selected records', value=True)}

toggle_buttons['SubstFromReact'].layout.width='375px'
# toggle_buttons['loopTthenP'].layout.width='375px'
toggle_buttons['ReactFromSubst'].layout.width='375px'

# toggle_buttons['loopTPthenRecords'].layout.width='375px'


op.substancePropertiesFromReaction = toggle_buttons['SubstFromReact'].value
# op.loopTemperatureThenPressure = toggle_buttons['loopTthenP'].value 
op.reactionPropertiesFromReactants = toggle_buttons['ReactFromSubst'].value
# op.loopOverTPpairsFirst = toggle_buttons['loopTPthenRecords'].value 

def on_toggle(**toggles):
    global op
    op.substancePropertiesFromReaction = toggle_buttons['SubstFromReact'].value
    # op.loopTemperatureThenPressure = toggle_buttons['loopTthenP'].value 
    op.reactionPropertiesFromReactants = toggle_buttons['ReactFromSubst'].value
    # op.loopOverTPpairsFirst = toggle_buttons['loopTPthenRecords'].value 

    #print(toggles)
    # do something with list of selected

interact_out = widgets.interactive_output(on_toggle, toggle_buttons)
#display(interact_out)

def ui(databases_w, dropdown_w):
  #  table_style = {'description_width': 'initial'}
    table_layout = {'width':'120px', 'min_width':'120px', 'height':'28px', 'min_height':'28px'}
    row_layout = {'width':'122px', 'min_width':'122px'}
    row_layout2 = {'width':'254px', 'min_width':'254px'}

    databases_w.layout = row_layout
    dropdown_w.layout = row_layout2

    t0_1_widget = widgets.BoundedIntText( description='Tmin',
                                        value=25,
                                        min=0,
                                        max=1000,
                                        step=1,
                                        disabled=False,
                                        style=style,

                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    def on_value_change_t0_1(change):
        temperatures[0] = change['new']
    temperatures[0] = t0_1_widget.value
    t0_1_widget.observe(on_value_change_t0_1, names='value')
    t0_2_widget = widgets.BoundedIntText( description='Tmax',
                                        value=150,
                                        min=0,
                                        max=1000,
                                        step=1,
                                        disabled=False,
                                        style=style,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    def on_value_change_t0_2(change):
        temperatures[1] = change['new']
    temperatures[1] = t0_2_widget.value
    t0_2_widget.observe(on_value_change_t0_2, names='value')
    t0_3_widget = widgets.BoundedIntText( description='Tstep',
                                        value=5,
                                        min=0,
                                        max=1000,
                                        step=1,
                                        disabled=False,
                                        style=style,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    def on_value_change_t0_3(change):
        temperatures[2] = change['new']
    temperatures[2] = t0_3_widget.value
    t0_3_widget.observe(on_value_change_t0_3, names='value')
    t1_1_widget = widgets.BoundedIntText( description='Pmin',
                                        value=0,
                                        min=0,
                                        max=100000,
                                        step=1,
                                        disabled=False,
                                        style=style,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    def on_value_change_t1_1(change):
        pressures[0] = change['new']
    pressures[0] = t1_1_widget.value
    t1_1_widget.observe(on_value_change_t1_1, names='value')
    t1_2_widget = widgets.BoundedIntText( description='Pmax',
                                        value=0,
                                        min=0,
                                        max=100000,
                                        step=1,
                                        disabled=False,
                                        style=style,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    def on_value_change_t1_2(change):
        pressures[1] = change['new']
    pressures[1] = t1_2_widget.value
    t1_2_widget.observe(on_value_change_t1_2, names='value')
    t1_3_widget = widgets.BoundedIntText( description='Pstep',
                                        value=0,
                                        min=0,
                                        max=100000,
                                        step=1,
                                        disabled=False,
                                        style=style,
                                    #layout=header_layout,
                                    #style=table_style
                                    layout=table_layout
                                   )
    def on_value_change_t1_3(change):
        pressures[2] = change['new']
    pressures[2] = t1_3_widget.value
    t1_3_widget.observe(on_value_change_t1_3, names='value')

    global progress
    hbox=widgets.HBox(children=[progress])
    # with progress_out:
    #     #display(widgets.IntSlider())
    #     display(hbox)
    
    # display(progress_out)

    progress_out.layout={'width':'380px', 'min_width':'380px'}

    hbox1 = widgets.HBox([databases_w, dropdown_w, t0_1_widget, t0_2_widget, t0_3_widget, toggle_buttons['SubstFromReact']])#, toggle_buttons['loopTthenP'] ])
    hbox2 = widgets.HBox([progress_out, t1_1_widget, t1_2_widget, t1_3_widget , toggle_buttons['ReactFromSubst']]) #, toggle_buttons['loopTPthenRecords'] ])
    ui = widgets.VBox([hbox1, hbox2])

    display(ui)

    return t0_1_widget, t0_2_widget, t0_3_widget, t1_1_widget, t1_2_widget, t1_3_widget

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

def stripComments(code):
    code = str(code)
    return re.sub(r'(?m)^ *#.*\n?', '', code)

def parse_reactions(text):
    textn = stripComments(text)
    postString = textn.split("\n")

    return postString

def load_widgets(dfDatabase, dfSelectSubst, dfSelectReact, dfInputReact) :
    dirName = 'databases'
    global db_file
    origin  = widgets.Dropdown(
    options = getListOfFiles(dirName),
    value   = db_file,
  #  description='Available Databases:',
    layout  = Layout(width= 'max-content')
    )
    # database = fun.Database(db_file)
    # select_subst = multi_checkbox_widget(database.mapSubstances())
    # select_react = multi_checkbox_widget(database.mapReactions())

    
    tabs = make_tabs_3( dfSelectSubst[db_file], dfSelectReact[db_file], dfInputReact[db_file], ['Substances', 'Reactions', 'Write Reactions'])
    
    plot_out=widgets.Output()
    link_out=widgets.Output()
    status_out=widgets.Output()

    box_layout2 = Layout(#overflow='scroll',
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
            origin.disabled = True
            button.disabled = True
            db_file = change['new']
            global database
            try:
                database = fun.Database(db_file)
                reaction_equations = []
                for r in database.mapReactions().values():
                     reaction_equations.append(r.equation())
                dfSelectReact[db_file] = multi_checkbox_widget(database.mapReactions(), ["logKr", "reaction_gibbs_energy","reaction_enthalpy","reaction_entropy","reaction_heat_capacity_cp", "reaction_volume"],
                                                             reaction_equations)
                dfSelectSubst[db_file] = multi_checkbox_widget(database.mapSubstances(), ["gibbs_energy","enthalpy","entropy","heat_capacity_cp","volume"])
                dfInputReact[db_file] = input_reactions_widget(database.mapSubstances(), ["logKr", "reaction_gibbs_energy","reaction_enthalpy","reaction_entropy","reaction_heat_capacity_cp", "reaction_volume"])
                # select_subst_n = multi_checkbox_widget(database.mapSubstances())
                # select_react_n = multi_checkbox_widget(database.mapReactions())
                # select_subst.children = dfSelectSubst[db_file].children
                # select_react.children = dfSelectReact[db_file].children
                #ndx = copy.deepcopy(tabs.selected_index)
                tabs_new = make_tabs_3(dfSelectSubst[db_file], dfSelectReact[db_file], dfInputReact[db_file], ['Substances', 'Reactions', 'Write Reactions'])
                tabs.children = tabs_new.children
                #tabs.selected_index=ndx

                with plot_out:
                    clear_output()
                with table_out:
                    clear_output()
            except Exception as e:
                with status_out:
                    clear_output(wait=True)
                    print(e)
            origin.disabled = False
            button.disabled = False
    origin.observe(on_dropdown_change)
    
    items = [Label(value='Available Databases:'), origin]

    # box_layout = Layout(
    #                 width='max-content',
    #                 height='',
    #                 flex_direction='row',
    #                 display='flex')
    # carousel = widgets.HBox(children=items, layout=box_layout)
    

    def on_button_clicked(b):
  #      with output:
  #      print("Button clicked.")
        origin.disabled = True
        button.disabled = True
        global db_file
        #display(select_subst)
        batch = fun.ThermoBatch(database)
        with status_out:
            clear_output(wait=True)
            print('Calculating...')    
        batch.setPropertiesUnits(["temperature", "pressure"],["degC","bar"])
        batch.setPropertiesDigits(["gibbs_energy","entropy", "volume",
                            "enthalpy","logKr", "temperature", "pressure"], [0, 4, 4, 4, 4, 0, 4])
        #properties = ["gibbs_energy", "enthalpy", "entropy"]
        
        temperature_pressure_pairs = [[50,1000],  [150,1000], [200,1000], [250,1000], [300,1000], [350,1000], 
                              [400,1000], [450,1000], [500,1000], [550,1000], [600,1000], [650,1000], 
                              [700,1000], [800,1000], [900,1000], [1000,1000]]
        global pressures
        global temperatures
        global op

        batch.setBatchPreferences(op)
        batch.setTemperatureIncrement(temperatures[0], temperatures[1], temperatures[2])
        batch.setPressureIncrement(pressures[0], pressures[1], pressures[2])

        #print(pressures)
        #print(temperatures)

        try:
            if tabs.selected_index == 0:
                selected_options = [w.description for w in dfSelectSubst[db_file].children[0].children[1].children if w.value]
                selected_properties = [w.description for w in dfSelectSubst[db_file].children[1].children[1].children if w.value]
                if (len(selected_options) == 0):
                    selected_options = list(database.mapSubstances().keys())
                batch.thermoPropertiesSubstance(selected_options, selected_properties).toCSV("results.csv")
            elif tabs.selected_index == 1:
                selected_options = [w.description for w in dfSelectReact[db_file].children[0].children[1].children if w.value]
                new_keys = []
                global reac_keys_dic
                for k in selected_options:
                    new_keys.append(reac_keys_dic[k])
                selected_options= new_keys
                selected_properties = [w.description for w in dfSelectReact[db_file].children[1].children[1].children if w.value]
                if (len(selected_options) == 0):
                    selected_options = list(database.mapReactions().keys())
                batch.thermoPropertiesReaction(selected_options, selected_properties).toCSV("results.csv")
            elif tabs.selected_index == 2:
                selected_properties = [w.description for w in dfInputReact[db_file].children[2].children[1].children if w.value]
                selected_options = parse_reactions(dfInputReact[db_file].children[1].value)
                op2 = fun.BatchPreferences()
                op2.isFixed = op.isFixed
                op2.substancePropertiesFromReaction=False
                op2.reactionPropertiesFromReactants=False
                batch.setBatchPreferences(op2)
                batch.thermoPropertiesReaction(selected_options, selected_properties).toCSV("results.csv")
            
    #        ax=plt.gca()

            with link_out:
                clear_output(wait=True)
                display(create_csv_download_link("results.csv"))
            
            with status_out:
                clear_output(wait=True)
                print('Reading the results')
            
            with plot_out:
                clear_output(wait=True)
                df = pd.read_csv('results.csv')
                plot_substances_properties_vs_temperature('results.csv', selected_options, df['P(bar)'])
                p = pd.unique(df['P(bar)'])
                plt.show()
                print(f'Pressure {p}(bar)')
            with table_out:
                clear_output(wait=True)
                df = pd.read_csv('results.csv')
                #display(df)
                table = TableDisplay(df)
                display(table)
                print("Select no record to calculate properties for all records:). ")
            with status_out:
                clear_output(wait=True)
                print('Finished')
        except Exception as e:
            with status_out:
                clear_output(wait=True)
                print(e)
        origin.disabled = False
        button.disabled = False
            
    button.on_click(on_button_clicked)

    #display(carousel)
    ui(Label(value='Available Databases:'), origin)
    display(tabs)
    #display(button)
    hbox=widgets.HBox(children=(button, link_out, status_out))
    display(hbox)
    
    
    #display(link_out)
    tabs_results = make_tabs(table_out, plot_out, ['Table', 'Plot'])
 #   hbox=widgets.HBox(children=(table_out, plot_out))
    display(tabs_results)

def initialize_widgets() :
    dfDatabase    = {}
    dfSelectSubst = {}
    dfSelectReact = {}
    dfInputReact = {}

    # global progress
    # hbox=widgets.HBox(children=[progress])
    # with progress_out:
    #     #display(widgets.IntSlider())
    #     display(hbox)
    
    # display(progress_out)

#     def on_button_clicked(b):
#         for i, f in enumerate(files):
#             progress.value = i+1
#             dfDatabase[f] = fun.Database(f)
#             reaction_equations = []
#             for r in dfDatabase[f].mapReactions().values():
#                 reaction_equations.append(r.equation())
#             dfSelectSubst[f] = multi_checkbox_widget(dfDatabase[f].mapSubstances(), ["gibbs_energy","enthalpy","entropy","heat_capacity_cp","volume"])
#             dfSelectReact[f] = multi_checkbox_widget(dfDatabase[f].mapReactions(), ["logKr", "reaction_gibbs_energy","reaction_enthalpy","reaction_entropy","reaction_heat_capacity_cp", "reaction_volume"],
#                                                         reaction_equations)
#             progress.value = i+1
#             #time.sleep(0.1)
#         with data_out:
#             clear_output(wait=True)
#             load_widgets(dfDatabase, dfSelectSubst, dfSelectReact)
            
#     button.on_click(on_button_clicked)
    global db_file
    try:
        database = fun.Database(db_file)
        reaction_equations = []
        for r in database.mapReactions().values():
            reaction_equations.append(r.equation())
        dfSelectReact[db_file] = multi_checkbox_widget(database.mapReactions(), ["logKr", "reaction_gibbs_energy","reaction_enthalpy","reaction_entropy","reaction_heat_capacity_cp", "reaction_volume"],
                                                            reaction_equations)
        dfSelectSubst[db_file] = multi_checkbox_widget(database.mapSubstances(), ["gibbs_energy","enthalpy","entropy","heat_capacity_cp","volume"])

        dfInputReact[db_file] = input_reactions_widget(database.mapSubstances(), ["logKr", "reaction_gibbs_energy","reaction_enthalpy","reaction_entropy","reaction_heat_capacity_cp", "reaction_volume"])
    except Exception as e:
        with progress_out:
            clear_output(wait=True)
            print(e)

    load_widgets(dfDatabase, dfSelectSubst, dfSelectReact, dfInputReact)


