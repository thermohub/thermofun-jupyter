import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt

from collections import OrderedDict
from itertools import cycle
from collections import OrderedDict
from matplotlib.lines import Line2D

from IPython.display import HTML
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


def plot_substances_properties_vs_temperature(results_csv_file, substances, pressure):

    # plot settings
    mpl.rcParams['lines.linewidth']=2
    mpl.rcParams['axes.labelsize']=20
    mpl.rcParams['axes.linewidth']=2
    mpl.rcParams['font.size']=18
    mpl.rcParams['figure.figsize']=[9,9]

    fig = plt.figure()


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

    return fig
           