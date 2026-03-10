import {LGraph, LGraphNode, LLink} from 'litegraph.js';

/** Manages the LGraph instance and graph-level operations */
export class GraphManager {
  graph: LGraph;
  private _onChangeCallbacks: (() => void)[] = [];

  constructor() {
    this.graph = new LGraph();
  }

  /** Clear the graph */
  clear(): void {
    this.graph.clear();
    this.notifyChange();
  }

  /** Get all nodes using LGraph's public findNodesByClass */
  getNodes(): LGraphNode[] {
    // LGraph stores nodes in private _nodes; use serialize to get them safely
    const data = this.graph.serialize();
    if (!data || !data.nodes) return [];
    return data.nodes.map((n: any) => this.graph.getNodeById(n.id)).filter(Boolean) as LGraphNode[];
  }

  /** Get all links */
  getLinks(): Record<number, LLink> {
    return this.graph.links || {};
  }

  /** Get node by ID */
  getNodeById(id: number): LGraphNode | undefined {
    return this.graph.getNodeById(id);
  }

  /** Register a change callback */
  onChange(callback: () => void): void {
    this._onChangeCallbacks.push(callback);
  }

  /** Notify all change listeners */
  notifyChange(): void {
    for (const cb of this._onChangeCallbacks)
      cb();
  }

  /** Serialize the graph */
  serialize(): any {
    return this.graph.serialize();
  }

  /** Load a graph from serialized data */
  configure(data: any): void {
    this.graph.configure(data);
    this.notifyChange();
  }

  /** Get node count */
  getNodeCount(): number {
    const data = this.graph.serialize();
    return data?.nodes?.length ?? 0;
  }

  /** Get link count */
  getLinkCount(): number {
    return Object.keys(this.getLinks()).length;
  }
}
