#!/usr/bin/env python3

""" Functions for creating Bokeh plot objects """

from Bio import SeqIO
from Bio.SeqUtils import GC
from bokeh.plotting import figure

def get_gc_from_file(filepath):
    """ Reads a (multi-)FASTA file and returns a dict as {id: GC-percentage, ...} """
    gc_percentage = {}
    for record in SeqIO.parse(filepath, "fasta"):
        gc_percentage[record.id] = GC(record.seq)

    return gc_percentage

def create_gc_plot(filepath):
    """ Creates a Bokeh plot object """

    # Calculate and get the GC-percentages
    percentages = get_gc_from_file(filepath)
    ids = list(percentages.keys())
    values = list(percentages.values())

    plot = figure(x_range=ids, plot_height=250, title="GC Percentages",
                  toolbar_location=None, tools="")

    plot.vbar(x=ids, top=values, width=0.9)
    plot.xgrid.grid_line_color = None
    plot.y_range.start = 0


    # TODO: replace x-axis labels with numbers and add a legend

    # TODO: format bar color depending on percentage

    # TODO: add tooltip showing sequence ID, length and GC-percentage

    # Return the plot object (instead of using 'show(plot)')
    return plot


if __name__ == "__main__":
    print(create_gc_plot('./data/rdp.fa'))
