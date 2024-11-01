import {indexFromEnd} from '../utils';
import {buildTraverseD} from './graph-traverse-utils';

export type NodePathSegment = {
  id: string,
  idx: number,
};

export type NodePath = NodePathSegment[];

export type NodeAddressSegment = {
  idx: number,
}

export type NodeAddress = NodeAddressSegment[];

export class BaseTree<T> {
  public root: TreeNode<T>;

  public traverse = buildTraverseD([] as NodePath, (item: TreeNode<T>, path: NodePath) => item.getChildren().map(({id, item}, idx) => [item, [...path, {id, idx}] as NodePath] as const));

  public static isNodeChildOrEq(path: Readonly<NodeAddress>, nodeAddress: Readonly<NodeAddress>): boolean {
    for (const [level, {idx}] of path.entries()) {
      const idx2 = nodeAddress[level]?.idx;
      if (idx !== idx2)
        return false;
    }
    return true;
  }

  public static isNodeChildOffseted(path: Readonly<NodeAddress>, nodeAddress: Readonly<NodeAddress>, childOffset?: number): boolean {
    if (childOffset == null)
      return this.isNodeChildOrEq(path, nodeAddress);
    if (!this.isNodeChildOrEq(path, nodeAddress))
      return false;
    const nextSegment = nodeAddress[path.length];
    if (!nextSegment || nextSegment.idx < childOffset)
      return false;
    return true;
  }

  public static isNodeAddressEq(a1: Readonly<NodeAddress>, a2: Readonly<NodeAddress>): boolean {
    for (const [level, {idx}] of a1.entries()) {
      const idx2 = a2[level]?.idx;
      if (idx !== idx2)
        return false;
    }
    if (a1.length !== a2.length)
      return false;
    return true;
  }

  constructor(item: T) {
    this.root = new TreeNode(item);
  }

  addItem(paddress: Readonly<NodeAddress>, item: T, id: string, idx: number) {
    const nodeSeq = this.getNodesFromAddress(paddress);
    const parent = indexFromEnd(nodeSeq)!;
    parent.addChild(item, id, idx);
    return item;
  }

  getItem(address: Readonly<NodeAddress>) {
    const node = this.getNode(address);
    return node.getItem();
  }

  getNode(address: Readonly<NodeAddress>) {
    const nodeSeq = this.getNodesFromAddress(address);
    const node = indexFromEnd(nodeSeq)!;
    return node;
  }

  attachBrunch(paddress: Readonly<NodeAddress>, node: TreeNode<T>, id: string, idx: number) {
    const nodeSeq = this.getNodesFromAddress(paddress);
    const parent = indexFromEnd(nodeSeq)!;
    parent.attachNode(node, id, idx);
    return node;
  }

  replaceBrunch(paddress: Readonly<NodeAddress>, node: TreeNode<T>, id: string, idx: number) {
    const nodeSeq = this.getNodesFromAddress(paddress);
    const parent = indexFromEnd(nodeSeq)!;
    parent.removeChild({idx});
    this.attachBrunch(paddress, node, id, idx);
  }

  removeBrunch(address: Readonly<NodeAddress>) {
    const nodeSeq = this.getNodesFromAddress(address);
    const parent = indexFromEnd(nodeSeq, 1)!;
    const segment = indexFromEnd(address)!;
    return parent.removeChild(segment);
  }

  find(pred: (item: T, path: NodePath) => boolean, root = this.root) {
    return this.traverse(root, ((acc, item, path, stop) => {
      if (pred(item.getItem(), path)) {
        stop();
        return [item, path] as const;
      }
      return acc;
    }), undefined as Readonly<[TreeNode<T>, NodePath]> | undefined);
  }

  private getNodesFromAddress(address: Readonly<NodeAddress>) {
    let current = this.root;
    const nodes = [current];
    for (const segment of address) {
      const node = current.getChild(segment);
      if (!node)
        throw new Error(`NodeTree: Failed to get all nodes to: ${JSON.stringify(address)}`);
      current = node;
      nodes.push(node);
    }
    return nodes;
  }
}

export class TreeNode<T> {
  private children = new PositionedMap<TreeNode<T>>();

  constructor(private item: T) {}

  public getChild(segment: NodeAddressSegment) {
    return this.children.getItemByIndex(segment.idx);
  }

  public getChildren() {
    return this.children.getAllItems();
  }

  public addChild(item: T, id: string, index: number) {
    const cnode = new TreeNode<T>(item);
    this.children.insertItem(cnode, id, index);
    return cnode;
  }

  public attachNode(node: TreeNode<T>, id: string, index: number) {
    this.children.insertItem(node, id, index);
    return node;
  }

  public removeChild(segment: NodeAddressSegment) {
    return this.children.removeItemByIndex(segment.idx);
  }

  public getItem() {
    return this.item;
  }
}

class PositionedMap<T> {
  public data: {item: T, id: string}[] = [];

  public getItemByIndex(idx: number) {
    return this.data[idx]?.item;
  }

  public getAllItems() {
    return [...this.data];
  }

  public insertItem(item: T, id: string, index: number) {
    const dataElem = {id, item};
    this.data.splice(index, 0, dataElem);
    return item;
  }

  public removeItemByIndex(idx: number) {
    const [removedItem] = this.data.splice(idx, 1);
    return removedItem?.item;
  }
}
