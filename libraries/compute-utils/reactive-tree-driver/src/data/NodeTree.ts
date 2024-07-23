import {ItemPathArray} from '../config/config-processing-utils';
import {indexFromEnd, pathJoin} from '../utils';


function traverseNodeTree<T, R>(
  startNode: TreeNode<T>,
  fn: (acc: R, item: T, pathArray: ItemPathArray, stop: () => void) => R,
  acc: R,
) {
  let stop = false;
  const signal = () => stop = true;
  const q: (readonly [ItemPathArray, TreeNode<T>])[] = [[[], startNode]];
  while (q.length) {
    const [path, item] = q.shift()!;
    acc = fn(acc, item.getItem(), path, signal);
    if (stop)
      return acc;
    const next = item.getChildren().map(({id, item}) => [pathJoin(path, [id]), item] as const);
    q.push(...next);
  }
}


export class NodeTree<T> {
  private root = new TreeNode<T>(this.item);

  constructor(public item: T) {}

  addItem(pathArray: ItemPathArray, item: T, index?: number) {
    const nodeSeq = this.getNodesFromPath(pathArray, true);
    const nodeId = indexFromEnd(pathArray)!;
    const parent = indexFromEnd(nodeSeq, 1)!;
    const node = indexFromEnd(nodeSeq);
    if (node)
      throw new Error(`NodeTree: Trying to overwrite an existing node: ${node}, path: ${pathArray}`);
    parent.addChild(nodeId, item, index);
    return item;
  }

  getItem(pathArray: ItemPathArray) {
    const nodeSeq = this.getNodesFromPath(pathArray);
    const node = indexFromEnd(nodeSeq)!;
    return node.getItem();
  }

  getItemPosition(pathArray: ItemPathArray) {
    const nodeSeq = this.getNodesFromPath(pathArray);
    const parent = indexFromEnd(nodeSeq, 1)!;
    const nodeId = indexFromEnd(pathArray)!;
    return parent.getChildPosition(nodeId);
  }

  getRoot() {
    return this.root;
  }

  attachBrunch(pathArray: ItemPathArray, index: number, brunch: NodeTree<T>) {
    const nodeSeq = this.getNodesFromPath(pathArray, true);
    const parent = indexFromEnd(nodeSeq, 1)!;
    const node = indexFromEnd(nodeSeq);
    const nodeId = indexFromEnd(pathArray)!;
    if (node)
      throw new Error(`NodeTree: Trying to overwrite an existing node: ${node}, path: ${pathArray}`);
    parent.attachNode(nodeId, brunch.root, index);
    return node;
  }

  removeBrunch(pathArray: ItemPathArray) {
    const nodeSeq = this.getNodesFromPath(pathArray);
    const parent = indexFromEnd(nodeSeq, 1)!;
    const segment = indexFromEnd(pathArray, 1)!;
    return parent.removeChild(segment);
  }

  toJSON() {

  }

  traverse<R>(fn: (acc: R, item: T | undefined, pathArray: ItemPathArray) => R, acc: R) {
    return traverseNodeTree(this.root, fn, acc!);
  }

  private getNodesFromPath(pathArray: ItemPathArray, allowEmptyTarget = false) {
    let current = this.root;
    const nodes = [current];
    for (const segment of pathArray) {
      const node = current.getChild(segment);
      if (!node)
        break;

      current = node;
      nodes.push(node);
    }
    const additionalNodes = allowEmptyTarget ? 0 : 1;
    if (nodes.length < pathArray.length + additionalNodes)
      throw new Error(`NodeTree: Failed to get all nodes to: ${pathArray}`);
    return nodes;
  }

  static fromJSON() {

  }
}


class TreeNode<T> {
  private children = new PositionedMap<TreeNode<T>>();

  constructor(private item: T) {}

  public getChild(segment: string) {
    return this.children.getItemById(segment);
  }

  public getChildPosition(segment: string) {
    return this.children.getItemPositionById(segment);
  }

  public getChildren() {
    return this.children.getAll();
  }

  public addChild(segment: string, item: T, index: number = this.children.data.length) {
    const cnode = new TreeNode<T>(item);
    this.children.insertItem(segment, index, cnode);
    return cnode;
  }

  public attachNode(segment: string, node: TreeNode<T>, index: number = this.children.data.length) {
    this.children.insertItem(segment, index, node);
    return node;
  }

  public removeChild(segment: string) {
    return this.children.removeItem(segment);
  }

  public getItem() {
    return this.item;
  }
}


class PositionedMap<T> {
  // no optimization for now
  public data: {item: T, id: string}[] = [];

  public getItemById(id: string) {
    const idx = this.findIndex(id);
    if (idx > 0)
      return this.data[idx]!.item;
  }

  public getItemPositionById(id: string) {
    const idx = this.findIndex(id);
    if (idx > 0)
      return idx;
  }

  public getAll() {
    return [...this.data];
  }

  public insertItem(id: string, index: number, item: T) {
    const currentIndex = this.findIndex(id);
    if (currentIndex >= 0)
      throw new Error(`PositionedMap: Adding duplicate id: ${id}`);
    const dataElem = {id, item};
    this.data.splice(index, 0, dataElem);
    return item;
  }

  public removeItem(id: string) {
    const currentIndex = this.findIndex(id);
    if (currentIndex < 0)
      throw new Error(`PositionedMap: Removing non-existent id: ${id}`);
    const [removedItem] = this.data.splice(currentIndex, 1);
    return removedItem.item!;
  }

  private findIndex(targetId: string) {
    return this.data.findIndex(({id}) => id === targetId);
  }
}
