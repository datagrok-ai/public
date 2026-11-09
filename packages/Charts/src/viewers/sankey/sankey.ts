/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import $ from 'cash-dom';

import {drag} from 'd3-drag';
import {ScaleOrdinal, scaleOrdinal} from 'd3-scale';
import {select} from 'd3-selection';
import * as graphlib from 'graphlib';

import {
  sankey,
  sankeyLinkHorizontal,
  sankeyCenter,
  sankeyJustify,
  sankeyLeft,
  sankeyRight,
  SankeyNode,
} from 'd3-sankey';

import {Observable} from 'rxjs';

import '../../../css/sankey-viewer.css';
import {MessageHandler, RenderSignals, subscribeToColumnRemoval, unsubscribeAll} from '../../utils/utils';
import {sankeyStatus} from './sankey-status';

interface Node {
  node: number,
  name: string,
  rowIdx: number,
}

interface Link {
  source: number,
  target: number,
  value: number,
  rowIdx: number,
}

interface Margin {
  top: number,
  right: number,
  bottom: number,
  left: number,
}

interface AlignMethod {
  center: (node: SankeyNode<{}, {}>, n: number) => number,
  justify: (node: SankeyNode<{}, {}>, n: number) => number,
  left: (node: SankeyNode<{}, {}>, n: number) => number,
  right: (node: SankeyNode<{}, {}>, n: number) => number,
}

@grok.decorators.viewer({
  name: 'Sankey',
  description: 'Creates a sankey viewer',
  icon: 'icons/sankey-viewer.svg',
})
export class SankeyViewer extends DG.JsViewer {
  sourceColumnName: string;
  targetColumnName: string;
  valueColumnName: string;
  initialized: boolean;

  graph: any;
  margin?: Margin;
  alignMethod?: AlignMethod;
  color?: ScaleOrdinal<string, number, never>;
  nodeWidth?: number;
  nodePadding?: number;

  nodeColorColumnName: string;
  linkColorColumnName: string;

  private signals = new RenderSignals();

  get onRendered(): Observable<void> {return this.signals.rendered;}

  get isRenderPending(): boolean {return this.signals.pending;}

  get renderError(): string | null {return MessageHandler._getMessage(this.root);}

  getWidgetStatus(): DG.IWidgetStatus {return sankeyStatus(this);}

  constructor() {
    super();
    // Properties
    this.sourceColumnName = this.string('sourceColumnName', null, {nullable: false, columnTypeFilter: DG.TYPE.STRING});
    this.targetColumnName = this.string('targetColumnName', null, {nullable: false, columnTypeFilter: DG.TYPE.STRING});
    this.valueColumnName = this.string('valueColumnName', null, {nullable: false, columnTypeFilter: DG.TYPE.NUMERICAL});
    this.nodeColorColumnName = this.string('nodeColorColumnName', null, {category: 'Color', columnTypeFilter: DG.TYPE.CATEGORICAL});
    this.linkColorColumnName = this.string('linkColorColumnName', null, {category: 'Color', columnTypeFilter: DG.TYPE.CATEGORICAL});
    this.addRowSourceAndFormula();

    this.initialized = false;
  }

  init() {
    // Data
    this.graph = null;

    // Chart Settings
    this.margin = {top: 10, right: 10, bottom: 10, left: 10};
    this.alignMethod = {center: sankeyCenter, justify: sankeyJustify,
      left: sankeyLeft, right: sankeyRight};
    this.color = scaleOrdinal(DG.Color.categoricalPalette);
    this.nodeWidth = 10;
    this.nodePadding = 15;
    this.initialized = true;
  }

  strColumnNames(): string[] {
    return this.dataFrame.columns.toList()
      .filter((col) => col.type === 'string' && col.categories.length <= 50).map((col) => col.name);
  }

  numColumnNames(): string[] {
    return this.dataFrame.columns.toList().filter((col) => ['double', 'int'].includes(col.type)).map((col) => col.name);
  }

  hasColumns(): boolean {
    return [this.sourceColumnName, this.targetColumnName, this.valueColumnName]
      .every((name) => this.dataFrame.col(name));
  }

  _testColumns() {
    return this.strColumnNames().length >= 2 && this.numColumnNames().length >= 1;
  }

  onTableAttached() {
    unsubscribeAll(this.subs);
    this.addSubs();
    this.init();

    if (this._testColumns()) {
      const strColumns = this.strColumnNames();
      this.sourceColumnName = strColumns[0];
      this.targetColumnName = strColumns[1];
      this.valueColumnName = this.numColumnNames()[0];
    }

    this.render();
  }

  addSubs() {
    this.subs.push(this.signals.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(this.signals.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.onMetadataChanged.subscribe((_) => this.render()));
    this.subs.push(subscribeToColumnRemoval(this, {
      sourceColumnName: () => this.strColumnNames().find((name) => name !== this.targetColumnName),
      targetColumnName: () => this.strColumnNames().find((name) => name !== this.sourceColumnName),
      valueColumnName: () => this.numColumnNames()[0],
      nodeColorColumnName: null,
      linkColorColumnName: null,
    }));
  }

  detach() {
    this.signals.reset();
    super.detach();
  }

  onSourceRowsChanged() {
    this.render();
  }

  prepareData() {
    if (!this._testColumns())
      return;

    const sourceCol = this.dataFrame.getCol(this.sourceColumnName);
    const targetCol = this.dataFrame.getCol(this.targetColumnName);
    const valueCol = this.dataFrame.getCol(this.valueColumnName);
    const linkColorCol = this.linkColorColumnName ? this.dataFrame.col(this.linkColorColumnName) : null;
    const sourceNames = sourceCol.categories;
    const targetNames = targetCol.categories;
    const sources = sourceCol.getRawData();
    const targets = targetCol.getRawData();
    const values = valueCol.getRawData();
    // a missing value is the int or float null, which would sink both nodes of the row
    const missing = valueCol.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL;
    const colors = linkColorCol?.type === DG.COLUMN_TYPE.STRING ? linkColorCol.getRawData() : null;
    const colorCount = linkColorCol === null ? 1 : colors !== null ? linkColorCol.categories.length : 2 ** 32;

    // one link per pair (and link color): d3-sankey's layout is superlinear in links, so a link per row freezes big tables
    const firstAsSource = new Int32Array(sourceNames.length).fill(-1);
    const firstAsTarget = new Int32Array(targetNames.length).fill(-1);
    const links = new Map<number, Link>();
    let negative = false;
    for (const row of this.filter.getSelectedIndexes()) {
      const s = sources[row];
      const t = targets[row];
      if (firstAsSource[s] < 0)
        firstAsSource[s] = row;
      if (firstAsTarget[t] < 0)
        firstAsTarget[t] = row;
      if (!sourceNames[s] || !targetNames[t])
        continue;

      const value = values[row] === missing ? 0 : values[row];
      negative ||= value < 0;
      const color = linkColorCol === null ? 0 : colors !== null ? colors[row] : DG.Color.getRowColor(linkColorCol, row) >>> 0;
      const key = (s * targetNames.length + t) * colorCount + color;
      const link = links.get(key);
      if (link)
        link.value += value;
      else
        links.set(key, {source: s, target: t, value: value, rowIdx: row});
    }

    // a name in both columns is one node; sources come first, each side in the order its rows come
    const nodes: Node[] = [];
    const nodeByName = new Map<string, Node>();
    for (const [names, first] of [[sourceNames, firstAsSource], [targetNames, firstAsTarget]] as const) {
      const shown = Array.from(first.keys()).filter((c) => first[c] >= 0 && names[c]).sort((a, b) => first[a] - first[b]);
      for (const c of shown) {
        const node = nodeByName.get(names[c]);
        if (node)
          node.rowIdx = Math.min(node.rowIdx, first[c]);
        else {
          nodes.push({node: nodes.length, name: names[c], rowIdx: first[c]});
          nodeByName.set(names[c], nodes[nodes.length - 1]);
        }
      }
    }

    for (const link of links.values()) {
      link.source = nodeByName.get(sourceNames[link.source])!.node;
      link.target = nodeByName.get(targetNames[link.target])!.node;
    }

    this.graph = {
      nodes: nodes,
      links: Array.from(links.values()),
      negative: negative,
    };
  }

  detectGraphCycles(graph: { nodes: Node[], links: Link[] }): string[][] {
    const g = new graphlib.Graph();

    for (const link of graph.links) {
      const source = graph.nodes[link.source].name;
      const target = graph.nodes[link.target].name;
      g.setEdge(source, target);
    }

    return graphlib.alg.findCycles(g);
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    if (this.initialized)
      this.render();
  }

  rowMatchesColumnValues(sourceCol: DG.Column | null, targetCol: DG.Column | null,
    rowIndex: number, sourceName: string, targetName: string, operator: 'AND' | 'OR' = 'AND',
  ): boolean {
    if (!sourceCol || !targetCol || !this.filter.get(rowIndex))
      return false;

    const sourceMatch = this.dataFrame.get(sourceCol.name, rowIndex) === sourceName;
    const targetMatch = this.dataFrame.get(targetCol.name, rowIndex) === targetName;

    switch (operator) {
    case 'AND': return sourceMatch && targetMatch;
    case 'OR': return sourceMatch || targetMatch;
    }
  }

  render() {
    this.signals.render(() => this.draw());
  }

  private draw() {
    $(this.root).empty();
    if (!this._testColumns()) {
      MessageHandler._showMessage(this.root, 'The Sankey viewer requires a minimum of 2 categorical (less than 50 unique categories) and 1 numerical columns.', 'd4-viewer-error');
      return;
    }
    if (!this.hasColumns()) {
      MessageHandler._showMessage(this.root, 'Select the Source, Target and Value columns.', 'd4-viewer-error');
      return;
    }

    this.prepareData();

    // d3-sankey sizes the nodes by their links: with none it throws (no rows) or yields NaN heights
    if (this.graph.links.length === 0)
      return;

    if (this.graph.negative) {
      MessageHandler._showMessage(this.root, `${this.valueColumnName} contains negative values`, 'd4-viewer-error');
      return;
    }
    if (this.graph.links.every((link: Link) => link.value === 0)) {
      MessageHandler._showMessage(this.root, `${this.valueColumnName}: all values are 0`, 'd4-viewer-error');
      return;
    }

    const cycles = this.detectGraphCycles(this.graph);
    if (cycles.length > 0) {
      MessageHandler._showMessage(this.root, 'The graph contains cycles. Please remove circular dependencies.', 'd4-viewer-error');
      return;
    }

    const width = this.root.parentElement!.clientWidth - this.margin!.left - this.margin!.right;
    const height = this.root.parentElement!.clientHeight - this.margin!.top - this.margin!.bottom;

    const generator = sankey().nodeWidth(this.nodeWidth!)
      .nodePadding(this.nodePadding!)
      .nodeAlign(this.alignMethod!.justify)
      .extent([[0, 0], [width, height]]);
    const graph = generator(this.graph!);

    const svg = select(this.root).append('svg')
      .attr('width', width + this.margin!.left + this.margin!.right)
      .attr('height', Math.abs(height + this.margin!.top + this.margin!.bottom))
      .append('g')
      .attr('transform', `translate(${this.margin!.left}, ${this.margin!.top})`);

    const nodeGroup = svg.append('g').attr('class', 'node');

    const dataFrameSourceColumn = this.dataFrame.getCol(this.sourceColumnName);
    const dataFrameTargetColumn = this.dataFrame.getCol(this.targetColumnName);
    const nodeColorCol = this.nodeColorColumnName ? this.dataFrame.col(this.nodeColorColumnName) : null;
    const linkColorCol = this.linkColorColumnName ? this.dataFrame.col(this.linkColorColumnName) : null;

    const nodes = nodeGroup
      .selectAll('rect')
      .data(graph.nodes)
      .join('rect')
      .attr('x', (d) => d.x0!)
      .attr('y', (d) => d.y0!)
      .attr('height', (d) => Math.abs(d.y1! - d.y0!))
      .attr('width', (d) => d.x1! - d.x0!)
      .attr('fill', (d: any) => {
        const color = nodeColorCol ? DG.Color.getRowColor(nodeColorCol, d.rowIdx) : this.color!(d.name);
        return DG.Color.toRgb(color);
      })
      .on('mouseover', (event, d: any) => {
        ui.tooltip.showRowGroup(this.dataFrame, (i) => {
          return this.rowMatchesColumnValues(dataFrameSourceColumn, dataFrameTargetColumn, i, d.name, d.name, 'OR');
        }, event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      //@ts-ignore
      .call(drag().subject((d) => d).on('drag', dragmove))
      .on('click', (event, d: any) => {
        if (event.defaultPrevented) return; // dragging
        this.dataFrame.selection.handleClick((i) => {
          return this.rowMatchesColumnValues(dataFrameSourceColumn, dataFrameTargetColumn, i, d.name, d.name, 'OR');
        }, event);
      });

    const links = svg.append('g')
      .selectAll('path')
      .data(graph.links)
      .join('path')
      .attr('class', 'link')
      .attr('d', sankeyLinkHorizontal())
      .style('stroke', (d: any) => {
        const color = linkColorCol ? DG.Color.getRowColor(linkColorCol, d.rowIdx) : DG.Color.darkGray;
        return DG.Color.toRgb(color);
      })
      .style('fill', 'none')
      .attr('stroke-width', (d) => Math.max(1, d.width!))
      .on('mouseover', (event, d: any) => {
        ui.tooltip.showRowGroup(this.dataFrame, (i) => {
          return this.rowMatchesColumnValues(dataFrameSourceColumn, dataFrameTargetColumn, i, d.source.name, d.target.name);
        }, event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      .on('click', (event, d: any) => {
        this.dataFrame.selection.handleClick((i) => {
          return this.rowMatchesColumnValues(dataFrameSourceColumn, dataFrameTargetColumn, i, d.source.name, d.target.name);
        }, event);
      });

    const titles = nodeGroup
      .selectAll('text')
      .data(graph.nodes)
      .join('text')
      .attr('x', (d) => d.x0! < width / 2 ? d.x1! + 6 : d.x0! - 6)
      .attr('y', (d) => (d.y1! + d.y0!) / 2)
      .attr('dy', '0.35em')
      .attr('text-anchor', (d) => d.x0! < width / 2 ? 'start' : 'end')
      .text((d: any) => d.name);

    function dragmove(this: any, event: any, d: any) {
      // this: Element???
      const rect = select(this);
      const rectX = +rect.attr('x');
      const rectY = +rect.attr('y');
      d.x0 += event.dx;
      d.x1 += event.dx;
      d.y0 += event.dy;
      d.y1 += event.dy;

      rect.attr('transform', `translate(${d.x0 - rectX}, ${d.y0 - rectY})`);
      titles.filter((t) => t === d)
        .attr('transform', `translate(${d.x0 - rectX}, ${d.y0 - rectY})`);

      generator.update(graph);
      links.attr('d', sankeyLinkHorizontal());
    };
  }
}
