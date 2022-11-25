#name: ScaffoldTreeGeneration
#description: generation scaffold tree from dataset
#language: python
#input: dataframe data [Input data table]
#input: string smiles 
#input: string names
#output: string result

!pip install scaffoldgraph
import scaffoldgraph as sg
import networkx as nx

#function that recursively adds child_nodes to hierarchies depending on the prev_scaffold value
def recurs_append_nodes(key, value, node, obj):
    if key in obj: 
        if obj[key] == value:
            obj['child_nodes'].append(node)
            return obj
    print(obj)
    for k, v in obj.items():
        if type(v) == list:
            result = recurs_append_nodes(key, value, node, v[0])

# function that returns nodes of each scaffold
def find_nodes(tree, nodes, scaffold):
    scaffold_nodes = []
    for i in range(len(nodes)):
        predecessors = list(nx.bfs_tree(tree, nodes[i], reverse=True))
        if predecessors[1] == scaffold:
            scaffold_nodes.append(nodes[i])
    return scaffold_nodes
    
#function that returns the hierarchy values
def get_hierarchies(tree):
    hierarchies = []
    counts = tree.get_hierarchy_sizes() 
    hierarchy_counts = sorted(counts.values())
    hierarchy_keys = sorted(counts.keys())
    for i in range(len(hierarchy_keys)):
        for j in range(hierarchy_counts[i]):
            hierarchies.append(hierarchy_keys[i])
    return hierarchies

#function that returns the previous scaffold
def get_prev_scaffold(tree, sorted_scaffolds, initial):
    val_dict = dict()
    for scaffold in sorted_scaffolds:
        val_dict[scaffold] = []
        for el in list(tree.successors(scaffold)):
            if el.isnumeric() is False:
                val_dict[scaffold].append(el)
    for k, v in val_dict.items():
        if initial in v:
            return k

#function that returns dict for each hierarchy depending on the input data
def get_hierarchy_dict(hierarchy_num, scaffold_str, child_nodes_list):
    hierarchy_dict = {
        'hierarchy': hierarchy_num, 
        'scaffold': scaffold_str, 
        'child_nodes': child_nodes_list
    }
    return hierarchy_dict

#function to get the json representation
def get_json_representation(tree):
    json_list = []
    scaffolds = list(tree.get_scaffold_nodes())
    sorted_scaffolds = sorted(scaffolds, key=len)
    hierarchies = get_hierarchies(tree)
    nodes = list(tree.get_molecule_nodes())
    json_list.append(get_hierarchy_dict(1, sorted_scaffolds[0], find_nodes(tree, nodes, sorted_scaffolds[0])))
    for i in range(1, len(sorted_scaffolds)):
        hierarchy_dict = get_hierarchy_dict(hierarchies[i], sorted_scaffolds[i], find_nodes(tree, nodes, sorted_scaffolds[i]))
        recurs_append_nodes('scaffold', get_prev_scaffold(tree, sorted_scaffolds, sorted_scaffolds[i]), hierarchy_dict, json_list[0])
    return json_list[0]

tree = sg.ScaffoldTree.from_dataframe(
    data, smiles_column=smiles, name_column=names, progress=True,
)
result = str(get_json_representation(tree))
