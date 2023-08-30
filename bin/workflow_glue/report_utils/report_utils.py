#!/usr/bin/env python
"""Create tables for the report."""

import numpy as np

RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def add_taxa(rank_dict, new_element=None):
    """
    Add new elements to a lineages structure dictionary.

    :param rank_dict (dict): Taxa counts in json structure. Nested dictionary.
        `{"Sample": {"Taxon name":{"rank":str, "value":int (opt), "children": dict},}}`.

    :return (dict)
    """
    if len(new_element) > 2:  # last element is species and next is the abundance
        if not rank_dict.get(new_element[0]):
            rank_dict[new_element[0]] = {'name': new_element[0], 'children': {}}
        return add_taxa(rank_dict[new_element[0]]['children'], new_element[1:])
    else:  # last element
        rank_dict[new_element[0]] = {
            'name': new_element[0], 'children': {}, 'value': new_element[1]}
        return rank_dict[new_element[0]]['children']


# update with counts
def itertaxa(d, new_d=[]):
    """Parse lineages structure dictionary to a list with the sunburst data structure.

    :param d (dict): Taxa counts in json structure. Nested dictionary.
            `{"Sample": {"Taxon name":{"rank":str, "value":int, "children": dict},}}`.
    :param new_d (list, optional): list of dictionaries, each one contains a taxon rank
        and their children taxa following the sunburst data. Defaults to [].
    :return (dict): list of dictionaries, in which dictionary represents a
            taxon rank and their children taxa following the sunburst data.
            structure.
    """
    for taxon, taxon_data in d.items():
        if taxon_data.get('value'):
            new_d.append(dict(name=taxon, value=int(taxon_data['value'])))
        else:
            new_d.append(dict(name=taxon))
        if bool(taxon_data["children"]):
            new_d[-1].update(children=[])
            itertaxa(taxon_data["children"], new_d[-1]["children"])
    return new_d


def prepare_data_to_sunburst(df_init):
    """Convert dataframe into a dictionary that can be use to plot a sunburst.

    :param df_init (pandsa Dataframe): It contains the relative abundance of each taxon.
    :return (list): List of dictionaries with a hierarchic structure:
        name of the node, value and children nodes.
    """
    df_ranks = df_init[RANKS]
    df_ranks = df_ranks.replace(np.nan, 'Unclassified', regex=True)
    df_ranks['abundances'] = df_init['estimated counts']
    multiindex = df_ranks.apply(tuple, axis=1)
    lineages = {}
    for r in multiindex.to_list():
        lineages.update(add_taxa(lineages, new_element=r))
    data2sunburst = itertaxa(lineages)
    return data2sunburst


def sum_terminal_nodes_in_list(data2sunburst, total=None):
    """Sum terminal nodes of each of the elements to be plotted.

    :param data2sunburst (list): List of dictionaries with a hierarchic structure.
    :param total (int, optional): Sum of the terminal nodes. Defaults to None.
    :return (int): Sum of the terminal nodes.
    """
    for i in data2sunburst:
        if i.get('children'):  # iterate to get the last level
            sum_terminal_nodes_in_list(i['children'], total)
        if (i.get('value')):
            total.append(i['value'])
    return sum(total)
