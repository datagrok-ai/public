import {indexFromEnd} from '../utils';
import {buildTraverseB} from './traversable';

export type NodeAddressSegment = {
  id: string,
  idx?: undefined,
} | {
  id?: undefined,
  idx: number
};

export type NodeAddress = NodeAddressSegment[];

export type NodePath = {
  id: string,
  idx: number,
}[];

export class NodeTree<T> {
  private root = new TreeNode(this.item);

  public traverse = buildTraverseB([] as NodePath, (item: TreeNode<T>, path: NodePath) => item.getChildren().map(({id, item}, idx) => [item, [...path, {id, idx}] as NodePath] as const));

  constructor(private item: T) {}

  addItem(paddress: Readonly<NodeAddress>, item: T, id: string, idx?: number) {
    const nodeSeq = this.getNodesFromAddress(paddress);
    const parent = indexFromEnd(nodeSeq)!;
    parent.addChild(item, id, idx);
    return item;
  }

  getItem(address: Readonly<NodeAddress>) {
    const nodeSeq = this.getNodesFromAddress(address);
    const node = indexFromEnd(nodeSeq)!;
    return node.getItem();
  }

  getRoot() {
    return this.root;
  }

  attachBrunch(paddress: Readonly<NodeAddress>, node: TreeNode<T>, id: string, idx?: number) {
    const nodeSeq = this.getNodesFromAddress(paddress);
    const parent = indexFromEnd(nodeSeq)!;
    parent.attachNode(node, id, idx);
    return node;
  }

  removeBrunch(address: Readonly<NodeAddress>) {
    const nodeSeq = this.getNodesFromAddress(address);
    const parent = indexFromEnd(nodeSeq, 1)!;
    const segment = indexFromEnd(address)!;
    return parent.removeChild(segment);
  }

  private getNodesFromAddress(address: Readonly<NodeAddress>) {
    let current = this.root;
    const nodes = [current];
    for (const segment of address) {
      const node = current.getChild(segment);
      if (!node)
        throw new Error(`NodeTree: Failed to get all nodes to: ${address}`); ;
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
    if (segment.idx != null)
      return this.children.getItemByIndex(segment.idx);
    else
      return this.children.getLastItemById(segment.id);
  }

  public getChildren() {
    return this.children.getAllItems();
  }

  public addChild(item: T, id: string = '', index: number = this.children.data.length) {
    const cnode = new TreeNode<T>(item);
    this.children.insertItem(cnode, id, index);
    return cnode;
  }

  public attachNode(node: TreeNode<T>, id: string, index: number = this.children.data.length) {
    this.children.insertItem(node, id, index);
    return node;
  }

  public removeChild(segment: NodeAddressSegment) {
    if (segment.idx != null)
      return this.children.removeItemByIndex(segment.idx);
    else
      return this.children.removeLastItemById(segment.id);
  }

  public getItem() {
    return this.item;
  }
}

class PositionedMap<T> {
  public data: {item: T, id: string}[] = [];

  public getLastItemById(id: string) {
    const idx = this.findLastIndex(id);
    if (idx > 0)
      return this.data[idx]?.item;
  }

  public getFirstItemById(id: string) {
    const idx = this.findIndex(id);
    if (idx > 0)
      return this.data[idx]?.item;
  }

  public getItemByIndex(idx: number) {
    return this.data[idx]?.item;
  }

  public getAllItems() {
    return [...this.data];
  }

  public insertItem(item: T, id: string, index = this.data.length - 1) {
    const dataElem = {id, item};
    this.data.splice(index, 0, dataElem);
    return item;
  }

  public removeLastItemById(id: string) {
    const currentIndex = this.findLastIndex(id);
    if (currentIndex < 0)
      throw new Error(`PositionedMap: Removing non-existent id: ${id}`);
    const [removedItem] = this.data.splice(currentIndex, 1);
    return removedItem.item!;
  }

  public removeFirstItemById(id: string) {
    const currentIndex = this.findIndex(id);
    if (currentIndex < 0)
      throw new Error(`PositionedMap: Removing non-existent id: ${id}`);
    const [removedItem] = this.data.splice(currentIndex, 1);
    return removedItem.item!;
  }

  public removeItemByIndex(idx: number) {
    const [removedItem] = this.data.splice(idx, 1);
    return removedItem?.item;
  }

  private findIndex(targetId: string) {
    return this.data.findIndex(({id}) => id === targetId);
  }

  private findLastIndex(targetId: string) {
    return this.data.findLastIndex(({id}) => id === targetId);
  }
}
