#name: ScaffoldTreeGeneration
#description: generation scaffold tree from dataset
#language: python
#input: dataframe data [Input data table]
#input: string smiles 
#input: string names
#output: string result

import scaffoldgraph as sg
import networkx as nx
import json

#function that recursively adds child_nodes to hierarchies depending on the prev_scaffold value
def recurs_append_nodes(key, value, node, obj):
    if key in obj: 
        if obj[key] == value:
            obj['child_nodes'].insert(0, node)
            return obj
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

#function that returns the scaffold parent
def get_parent(tree, scaffold):
    return ''.join(tree.get_parent_scaffolds(scaffold, max_levels=1))

#function that returns the right order of scaffolds
def get_sorted_scaffolds(tree, scaffolds):
    child_scaffolds = [scaffold for scaffold in scaffolds if len(tree.get_child_scaffolds(scaffold, max_levels=1)) == 0]
    sorted_scaffolds = []
    for i in range(0, len(child_scaffolds)):
        result = []
        result.insert(0, child_scaffolds[i])
        child = child_scaffolds[i]
        while len(get_parent(tree, child)) != 0:
            child = get_parent(tree, child)
            if child not in sorted_scaffolds:
                result.insert(0, child)
        sorted_scaffolds.extend(result)
    return sorted_scaffolds

#function that returns scaffold hierarchies
def get_hierarchies_list(tree, sorted_scaffolds):
    result = []
    for scaffold in sorted_scaffolds:
        result.append(tree.nodes[scaffold]['hierarchy'])
    return result

#function that returns dict for each hierarchy depending on the input data
def get_hierarchy_dict(scaffold_str, child_nodes_list):
    hierarchy_dict = {
        'scaffold': scaffold_str, 
        'child_nodes': child_nodes_list
    }
    return hierarchy_dict

#function to get the json representation
def get_json_representation(tree):
    json_list = []
    scaffolds = list(tree.get_scaffold_nodes())
    sorted_scaffolds = get_sorted_scaffolds(tree, scaffolds)
    hierarchies = get_hierarchies_list(tree, sorted_scaffolds)
    nodes = list(tree.get_molecule_nodes())
    json_list.append(get_hierarchy_dict(sorted_scaffolds[0], find_nodes(tree, nodes, sorted_scaffolds[0])))
    for i in range(1, len(sorted_scaffolds)):
        hierarchy_dict = get_hierarchy_dict(sorted_scaffolds[i], [])
        molecule_nodes = find_nodes(tree, nodes, sorted_scaffolds[i])
        for j in range(0, len(molecule_nodes)):
            hierarchy_dict['child_nodes'].append(get_hierarchy_dict(molecule_nodes[j], []))
        recurs_append_nodes('scaffold', ''.join(tree.get_parent_scaffolds(sorted_scaffolds[i], max_levels=1)), hierarchy_dict, json_list[0])
    return json_list[0]

tree = sg.ScaffoldTree.from_dataframe(
    data, smiles_column=smiles, name_column=names, progress=True,
)

res = get_json_representation(tree)
result = json.dumps(res)
