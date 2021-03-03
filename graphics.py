#!/usr/bin/env python3

""" Functions for creating Bokeh plot objects """

from Bio import SeqIO
from Bio.SeqUtils import GC
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, CustomJS, TapTool, HoverTool


def get_gc_from_file(filepath):
    """ Reads a (multi-)FASTA file and returns a dict as {id: GC-percentage, ...} """
    gc_percentage = {}
    for record in SeqIO.parse(filepath, "fasta"):
        gc_percentage[record.id] = GC(record.seq)

    return gc_percentage


def create_gc_plot_old(filepath):
    """ Creates a Bokeh plot object """

    # Calculate and get the GC-percentages
    percentages = get_gc_from_file(filepath)
    ids = list(percentages.keys())
    values = list(percentages.values())

    source = ColumnDataSource(data=dict(ids=ids, values=values))

    plot = figure(x_range=ids, plot_height=250, title="GC Percentages",
                  toolbar_location=None)

    plot.vbar(x=ids, top=values, width=0.9)
    plot.xgrid.grid_line_color = None
    plot.y_range.start = 0

    # TODO: replace x-axis labels with numbers and add a legend

    # TODO: format bar color depending on percentage

    # TODO: add tooltip showing sequence ID, length and GC-percentage

    # Return the plot object (instead of using 'show(plot)')
    return plot


def create_gc_plot(filename):
    # Create the data for the points
    x1 = [0, 1, 2, 3]
    y1 = [0, 1, 0, 1]
    x2 = [0, 1, 2, 3]
    y2 = [0.5, 0, 0.5, 0]
    colors = ['blue', 'red']

    sources = [ColumnDataSource(data=dict(x=x1, y=y1)),
               ColumnDataSource(data=dict(x=x2, y=y2))]

    # Add tools to the plot
    tap = TapTool()
    hover = HoverTool(tooltips=[("X", "@x"),
                                ("Y", "@y")])

    # Create a plotting figure
    p = figure(plot_width=400, plot_height=400, tools=[tap, hover])

    # Plots circles
    for i, source in enumerate(sources):
       p.circle('x', 'y', source=source, size=25, alpha=1, color=colors[i], hover_color='black', hover_alpha=1)

    # Create a CustomJS callback with the code and the data
    tap.callback = CustomJS(args=dict(sources=sources), code="""
        // Get indices array of all selected items
        for (let i = 0; i < sources.length; i++) { 
            let selected = sources[i].selected.indices;
            let data = sources[i].data;
            if (selected.length > 0) {
                console.log(selected + " selected from data source " + i);
                console.log("values: [" + data['x'][selected[0]] + ", " + data['y'][selected[0]] + "]");
                get_new_plot(i, selected);
            }
        };
    """)

    return p


if __name__ == "__main__":
    print(create_gc_plot('./data/rdp.fa'))
