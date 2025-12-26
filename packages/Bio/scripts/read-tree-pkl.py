# coding: utf-8

import ete3
import numpy as np
import pandas as pd
import pickle
import six
from sklearn.cluster import AgglomerativeClustering
from Bio import Phylo as ph
from Bio.Phylo.BaseTree import TreeElement
from Bio.Phylo import parse

import click
from click_default_group import DefaultGroup


class TreeNode:

    def __init__(self, node_id: int, distance: float, children: list['TreeNode']):
        """
        :param id:
        :param distance:
        :param children: list of child nodes, set None for leaf node
        """
        self._node_id: int = node_id
        self._distance: float = distance
        self._children: [] = children

    def set_height(self, height):
        """ Sets distance based on overal height of children"""
        self._distance = height - max((ch.get_height() for ch in self._children))

    def get_height(self):
        if self._children is None or len(self._children) == 0:
            return self._distance
        else:
            return self._distance + max(ch.get_height() for ch in self._children)

    def __str__(self):
        if self._children is None or len(self._children) == 0:
            return "{id}:{distance}".format(id=self._node_id, distance=self._distance)
        else:
            return "({children}){id}:{distance}".format(
                children=','.join((str(ch) for ch in self._children)),
                id=self._node_id,
                distance=self._distance)


@click.group(cls=DefaultGroup, default='main')
def cli():
    pass


@cli.command()
@click.option('--pkl', 'pkl_f',
              help="Pickle file with tree data",
              type=click.File('rb'))
@click.option('--nwk', 'nwk_f',
              help="Output file for tree in newick format",
              type=click.File('w'))
@click.pass_context
def main(ctx, pkl_f, nwk_f):
    data = pickle.load(pkl_f)

    all_nodes = {}

    # https://stackoverflow.com/questions/27386641/how-to-traverse-a-tree-from-sklearn-agglomerativeclustering
    # The hint from user76284
    v = dict(enumerate(data.children_, data.n_leaves_))

    root: TreeNode = None  # The last created node is a root
    for row_i in range(data.children_.shape[0]):
        height = data.distances_[row_i]  # height of connection - distance between nodes (groups) -

        node_id = data.n_leaves_ + row_i
        child_id_list = data.children_[row_i, :].tolist()

        child_node_list = []
        for child_id in child_id_list:
            child_node: TreeNode = None
            if child_id < data.n_leaves_:
                # leaf
                # create child node on height
                child_node = TreeNode(child_id, height, None)
                pass
            else:
                # internal node
                child_node = all_nodes[child_id]
                child_node.set_height(height);
                pass

            child_node_list.append(child_node)

        # create node with distance 0, set distance to parent node later based on connection height
        root = TreeNode(node_id, 0, child_node_list)
        all_nodes[node_id] = root

        k = 11

    nwk_f.write(str(root))
    nwk_f.write(';')


def build_Newick_tree(children, n_leaves, X, leaf_labels, spanner):
    """
    build_Newick_tree(children,n_leaves,X,leaf_labels,spanner)

    Get a string representation (Newick tree) from the sklearn
    AgglomerativeClustering.fit output.

    Input:
        children: AgglomerativeClustering.children_
        n_leaves: AgglomerativeClustering.n_leaves_
        X: parameters supplied to AgglomerativeClustering.fit
        leaf_labels: The label of each parameter array in X
        spanner: Callable that computes the dendrite's span

    Output:
        ntree: A str with the Newick tree representation

    """
    return go_down_tree(children, n_leaves, X, leaf_labels, len(children) + n_leaves - 1, spanner)[0] + ';'


def go_down_tree(children, n_leaves, X, leaf_labels, nodename, spanner):
    """
    go_down_tree(children,n_leaves,X,leaf_labels,nodename,spanner)

    Iterative function that traverses the subtree that descends from
    nodename and returns the Newick representation of the subtree.

    Input:
        children: AgglomerativeClustering.children_
        n_leaves: AgglomerativeClustering.n_leaves_
        X: parameters supplied to AgglomerativeClustering.fit
        leaf_labels: The label of each parameter array in X
        nodename: An int that is the intermediate node name whos
            children are located in children[nodename-n_leaves].
        spanner: Callable that computes the dendrite's span

    Output:
        ntree: A str with the Newick tree representation

    """
    nodeindex = nodename - n_leaves
    if nodename < n_leaves:
        return leaf_labels[nodeindex], np.array([X[nodeindex]])
    else:
        node_children = children[nodeindex]
        branch0, branch0samples = go_down_tree(children, n_leaves, X, leaf_labels, node_children[0])
        branch1, branch1samples = go_down_tree(children, n_leaves, X, leaf_labels, node_children[1])
        node = np.vstack((branch0samples, branch1samples))
        branch0span = spanner(branch0samples)
        branch1span = spanner(branch1samples)
        nodespan = spanner(node)
        branch0distance = nodespan - branch0span
        branch1distance = nodespan - branch1span
        nodename = '({branch0}:{branch0distance},{branch1}:{branch1distance})'.format(branch0=branch0,
                                                                                      branch0distance=branch0distance,
                                                                                      branch1=branch1,
                                                                                      branch1distance=branch1distance)
        return nodename, node


def get_cluster_spanner(aggClusterer):
    """
    spanner = get_cluster_spanner(aggClusterer)

    Input:
        aggClusterer: sklearn.cluster.AgglomerativeClustering instance

    Get a callable that computes a given cluster's span. To compute
    a cluster's span, call spanner(cluster)

    The cluster must be a 2D numpy array, where the axis=0 holds
    separate cluster members and the axis=1 holds the different
    variables.

    """
    if aggClusterer.linkage == 'ward':
        if aggClusterer.affinity == 'euclidean':
            spanner = lambda x: np.sum((x - aggClusterer.pooling_func(x, axis=0)) ** 2)
    elif aggClusterer.linkage == 'complete':
        if aggClusterer.affinity == 'euclidean':
            spanner = lambda x: np.max(np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=2))
        elif aggClusterer.affinity == 'l1' or aggClusterer.affinity == 'manhattan':
            spanner = lambda x: np.max(np.sum(np.abs(x[:, None, :] - x[None, :, :]), axis=2))
        elif aggClusterer.affinity == 'l2':
            spanner = lambda x: np.max(np.sqrt(np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=2)))
        elif aggClusterer.affinity == 'cosine':
            spanner = lambda x: np.max(np.sum((x[:, None, :] * x[None, :, :])) / (
                    np.sqrt(np.sum(x[:, None, :] * x[:, None, :], axis=2, keepdims=True)) * np.sqrt(
                np.sum(x[None, :, :] * x[None, :, :], axis=2, keepdims=True))))
        else:
            raise AttributeError('Unknown affinity attribute value {0}.'.format(aggClusterer.affinity))
    elif aggClusterer.linkage == 'average':
        if aggClusterer.affinity == 'euclidean':
            spanner = lambda x: np.mean(np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=2))
        elif aggClusterer.affinity == 'l1' or aggClusterer.affinity == 'manhattan':
            spanner = lambda x: np.mean(np.sum(np.abs(x[:, None, :] - x[None, :, :]), axis=2))
        elif aggClusterer.affinity == 'l2':
            spanner = lambda x: np.mean(np.sqrt(np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=2)))
        elif aggClusterer.affinity == 'cosine':
            spanner = lambda x: np.mean(np.sum((x[:, None, :] * x[None, :, :])) / (
                    np.sqrt(np.sum(x[:, None, :] * x[:, None, :], axis=2, keepdims=True)) * np.sqrt(
                np.sum(x[None, :, :] * x[None, :, :], axis=2, keepdims=True))))
        else:
            raise AttributeError('Unknown affinity attribute value {0}.'.format(aggClusterer.affinity))
    else:
        raise AttributeError('Unknown linkage attribute value {0}.'.format(aggClusterer.linkage))
    return spanner


if __name__ == '__main__':
    cli()
