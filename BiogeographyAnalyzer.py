from ete3 import Tree, NodeStyle, TreeStyle, faces

# color for each location
NODE_COLORS = {
    'SA': "#009E73", 'NA': "#0072B2",
    'SEA': "#D55E00", 'OC': "#E69F00",
    'AF': "#F0E442", 'MAD': "#CC79A7",
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

def add_pie_chart_node(node, specific):
    '''
    Adds pie chart face to internal node
    :param node: node to add pie chart to
    :param specific: True if by specific location or False if broader classification
    '''
    biogeo = node.biogeo
    locations = []
    colors = []
    diameter = 30

    if specific:
        # Go through each location and add it to piechart list and its probability
        for location in biogeo:
            locations.append(100 * biogeo[location])
            colors.append(NODE_COLORS[location])
    else:
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

        # add locations and probabilities
        locations.append(100 * asian_pacific)
        colors.append(NODE_COLORS['asian-pacific'])
        locations.append(100 * americas)
        colors.append(NODE_COLORS['americas'])
        locations.append(100 * africa)
        colors.append(NODE_COLORS['africa'])

    # put piechart as face and get rid of normal node circle
    node.img_style["size"] = 0
    node.add_face(faces.PieChartFace(locations, diameter, diameter, colors=colors), 0)

def add_pie_chart_leaf(node, specific):
    '''
    Adds a 'pie chart' to a leaf, but since it is only form one lcoation it just changes color of node
    :param node: leaf to color
    :param specific: True if by specific location or False if broader classification
    '''
    nstyle = NodeStyle()
    if specific:
        nstyle['fgcolor'] = get_color(node.biogeo)
    else:
        nstyle['fgcolor'] = get_color_broad(node.biogeo)
    nstyle['size'] = 20
    node.set_style(nstyle)

def add_pie_chart_all(tree, specific):
    '''
    Adds pie chart to whole of tree
    :param tree: Tree to color and add pie charts to
    :param specific: True if by specific location or False if broader classification
    '''
    for node in tree.traverse():
        if node.is_leaf():
            add_pie_chart_leaf(node, specific)
        else:
            add_pie_chart_node(node, specific)

if __name__ == '__main__':
    # load tree from file containing newick
    tree = Tree('tree.nwk')
    tree.set_outgroup('Falco_peregrinus_U83307.1')

    # Calculate biogeography chance and color nodes
    tree.ladderize(direction=0)
    parrot_node = tree.children[1]
    get_biogeography(parrot_node)
    add_pie_chart_all(parrot_node, False)

    # extra tree styling
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.rotation = -90

    # ts.mode = "c"
    # ts.arc_start = -180
    # ts.arc_span = 180

    # tree.convert_to_ultrametric()

    # show tree
    tree.show(tree_style=ts)
