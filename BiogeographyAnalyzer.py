from ete3 import Tree, NodeStyle, TreeStyle

# color for each location
NODE_COLORS = {
    'SA': "green",
    'NA': "blue",
    'SEA': "red",
    'OC': "orange",
    'AF': "brown",
    'MAD': "chocolate",
    'americas': '#1E88E5', 'asian-pacific':'#D81B60', 'africa': '#FFC107',
}

def get_biogeography(node):
    '''
    Find biogeographic information of ete3 node and places it in 'biogeo' feature
    :param node: the ete3 node to get the biogeography information of
    :return: a dictionary with format {location: chance of originating at location}
    '''

    # if a leaf, find location and return 100% chance to be at said location
    if node.is_leaf():
        name = node.name
        # leaves are named <Genus>_<species>_<accession>_<location>
        index = len(name) - 1
        while name[index] != '_':
            index -= 1

        # add 100% chance to originate at location
        node.add_feature('biogeo', {name[index + 1:]: 1})

        # change name to be 'Genus species'
        space = node.name.find('_')
        cutoff = node.name[space + 1:].find('_')
        node.name = node.name[:space] + ' ' + node.name[space + 1:][:cutoff]

        return {name[index + 1:]: 1}
    # if not leaf then base biogeography on child nodes
    else:
        biogeo = {}
        num_child = len(node.children)

        # go through each child and add biogeography info into dictionary. Each child is weighted the same.
        for child in node.children:
            child_biogeo = get_biogeography(child)
            for location in child_biogeo:
                if location in biogeo:
                    biogeo[location] += child_biogeo[location] / num_child
                else:
                    biogeo[location] = child_biogeo[location] / num_child

        # round for easy parsing of information
        for location in biogeo:
            biogeo[location] = round(biogeo[location], 2)

        # sort dictionary by highest probability and add biogeo feature to node
        node.add_feature('biogeo', dict(sorted(biogeo.items(), key=lambda x:x[1], reverse=True)))
        return biogeo

def get_color(biogeo):
    '''
    Gets color of node based on biogeo information
    :param biogeo: dictionary holding locations as keys and probabilities as values
    :return: Color corresponding to most probable location
    '''
    max = 0;
    loc = ''
    for location in biogeo:
        if biogeo[location] > max:
            loc = location
            max = biogeo[location]
    return NODE_COLORS[loc]

def get_color_broad(biogeo):
    '''
    Gets color of node based on a broader category of biogeo information
    :param biogeo: dictionary holding locations as keys and probabilities as values
    :return: Color corresponding to most probable location
    '''

    # get probabilities of each specific region
    sea = biogeo['SEA'] if 'SEA' in biogeo else 0
    oc = biogeo['OC'] if 'OC' in biogeo else 0
    sa = biogeo['SA'] if 'SA' in biogeo else 0
    na = biogeo['NA'] if 'NA' in biogeo else 0
    af = biogeo['AF'] if 'AF' in biogeo else 0
    mad = biogeo['MAD'] if 'MAD' in biogeo else 0

    # get combined scores of broader regions
    asian_pacific = sea + oc
    americas = sa + na
    africa = af + mad

    # get highest probability and return color
    score_dict = {asian_pacific: 'asian-pacific', americas: 'americas', africa: 'africa'}
    winner = score_dict[max(asian_pacific, americas, africa)]
    return NODE_COLORS[winner]

def color_nodes(tree, specific):
    '''
    Color the nodes based on biogeographic probabilities
    :param tree: The ete3 tree to color
    :param specific: True to color by specific regions and False to do the broader classification
    '''

    # goes through each  node in tree
    for node in tree.traverse():
        # make node style object
        nstyle = NodeStyle()
        if specific:
            nstyle['fgcolor'] = get_color(node.biogeo)
        else:
            nstyle['fgcolor'] = get_color_broad(node.biogeo)
        # make leaves bigger than internal nodes
        nstyle['size'] = 15 if node.is_leaf() else 10
        node.set_style(nstyle)

if __name__ == '__main__':
    # load tree from file containing newick
    tree = Tree('tree.nwk')

    # Calculate biogeography chance and color nodes
    tree.ladderize(direction=1)
    get_biogeography(tree)
    color_nodes(tree, False)

    # extra tree styling
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -180
    ts.arc_span = 180
    tree.img_style["size"] = 30

    # show tree
    tree.show(tree_style=ts)
