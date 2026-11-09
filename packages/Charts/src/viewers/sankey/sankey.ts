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
  color: number,
}

interface CategoryCodes {
  codes: Int32Array,
  names: string[],
}

interface SankeyData {
  nodes: Node[],
  links: Link[],
  negative: boolean,
  source: CategoryCodes,
  target: CategoryCodes,
  sourceNodes: Int32Array,
  targetNodes: Int32Array,
  linkColors: CategoryCodes | null,
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

  static categoryCodes(col: DG.Column, rows: Int32Array): CategoryCodes {
    if (col.type === DG.COLUMN_TYPE.STRING && !col.isVirtual)
      return {codes: col.getRawData() as Int32Array, names: col.categories};

    const codes = new Int32Array(col.length).fill(-1);
    const codeByName = new Map<string, number>();
    for (const row of rows) {
      const name = col.getString(row);
      let code = codeByName.get(name);
      if (code === undefined)
        codeByName.set(name, code = codeByName.size);
      codes[row] = code;
    }
    return {codes: codes, names: Array.from(codeByName.keys())};
  }

  static nodesInRowOrder(sides: [names: string[], firstRows: Int32Array][]): {nodes: Node[], nodeOf: Int32Array[]} {
    const nodes: Node[] = [];
    const byName = new Map<string, Node>();
    const nodeOf = sides.map(([names]) => new Int32Array(names.length).fill(-1));
    for (let side = 0; side < sides.length; side++) {
      const [names, firstRows] = sides[side];
      const shown = Array.from(firstRows.keys())
        .filter((c) => firstRows[c] >= 0 && names[c])
        .sort((a, b) => firstRows[a] - firstRows[b]);
      for (const c of shown) {
        let node = byName.get(names[c]);
        if (node)
          node.rowIdx = Math.min(node.rowIdx, firstRows[c]);
        else {
          node = {node: nodes.length, name: names[c], rowIdx: firstRows[c]};
          nodes.push(node);
          byName.set(node.name, node);
        }
        nodeOf[side][c] = node.node;
      }
    }
    return {nodes: nodes, nodeOf: nodeOf};
  }

  private prepareData(): SankeyData {
    const rows = this.filter.getSelectedIndexes();
    const source = SankeyViewer.categoryCodes(this.dataFrame.getCol(this.sourceColumnName), rows);
    const target = SankeyViewer.categoryCodes(this.dataFrame.getCol(this.targetColumnName), rows);
    const linkColorCol = this.linkColorColumnName ? this.dataFrame.col(this.linkColorColumnName) : null;
    const linkColors = linkColorCol ? SankeyViewer.categoryCodes(linkColorCol, rows) : null;
    const valueCol = this.dataFrame.getCol(this.valueColumnName);
    const raw = !valueCol.isVirtual && (valueCol.type === DG.COLUMN_TYPE.INT || valueCol.type === DG.COLUMN_TYPE.FLOAT) ?
      valueCol.getRawData() : null;
    // a missing value is the int or float null, which would sink both nodes of the row
    const missing = valueCol.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL;
    const valueAt = raw !== null ?
      (row: number) => raw[row] === missing ? 0 : raw[row] :
      (row: number) => valueCol.isNone(row) ? 0 : Number(valueCol.get(row));

    const firstAsSource = new Int32Array(source.names.length).fill(-1);
    const firstAsTarget = new Int32Array(target.names.length).fill(-1);
    const groups = new Map<number, {sourceCode: number, targetCode: number, value: number, rowIdx: number, color: number}>();
    let negative = false;
    for (const row of rows) {
      const s = source.codes[row];
      const t = target.codes[row];
      if (firstAsSource[s] < 0)
        firstAsSource[s] = row;
      if (firstAsTarget[t] < 0)
        firstAsTarget[t] = row;
      if (!source.names[s] || !target.names[t])
        continue;

      const value = valueAt(row);
      negative ||= value < 0;
      const color = linkColors?.codes[row] ?? 0;
      const key = (s * target.names.length + t) * (linkColors?.names.length ?? 1) + color;
      const group = groups.get(key);
      if (group)
        group.value += value;
      else
        groups.set(key, {sourceCode: s, targetCode: t, value: value, rowIdx: row, color: color});
    }

    const {nodes, nodeOf: [sourceNodes, targetNodes]} =
      SankeyViewer.nodesInRowOrder([[source.names, firstAsSource], [target.names, firstAsTarget]]);
    const links: Link[] = Array.from(groups.values(), (group) => ({
      source: sourceNodes[group.sourceCode],
      target: targetNodes[group.targetCode],
      value: group.value,
      rowIdx: group.rowIdx,
      color: group.color,
    }));

    return {nodes: nodes, links: links, negative: negative, source: source, target: target,
      sourceNodes: sourceNodes, targetNodes: targetNodes, linkColors: linkColors};
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

    const data = this.prepareData();

    // d3-sankey sizes the nodes by their links: with none it throws (no rows) or yields NaN heights
    if (data.links.length === 0)
      return;

    if (data.negative) {
      MessageHandler._showMessage(this.root, `${this.valueColumnName} contains negative values`, 'd4-viewer-error');
      return;
    }
    if (data.links.every((link) => link.value === 0)) {
      MessageHandler._showMessage(this.root, `${this.valueColumnName}: all values are 0`, 'd4-viewer-error');
      return;
    }

    const cycles = this.detectGraphCycles(data);
    if (cycles.length > 0) {
      MessageHandler._showMessage(this.root, 'The graph contains cycles. Please remove circular dependencies.', 'd4-viewer-error');
      return;
    }

    const width = this.root.parentElement!.clientWidth - this.margin!.left - this.margin!.right;
    const height = this.root.parentElement!.clientHeight - this.margin!.top - this.margin!.bottom;

    const generator = sankey<Node, Link>().nodeWidth(this.nodeWidth!)
      .nodePadding(this.nodePadding!)
      .nodeAlign(this.alignMethod!.justify)
      .extent([[0, 0], [width, height]]);
    const graph = generator(data);

    const svg = select(this.root).append('svg')
      .attr('width', width + this.margin!.left + this.margin!.right)
      .attr('height', Math.abs(height + this.margin!.top + this.margin!.bottom))
      .append('g')
      .attr('transform', `translate(${this.margin!.left}, ${this.margin!.top})`);

    const nodeGroup = svg.append('g').attr('class', 'node');

    const {source, target, sourceNodes, targetNodes, linkColors} = data;
    const filter = this.filter.getBuffer();
    const passes = (i: number) => (filter[i >>> 5] & (1 << (i & 31))) !== 0;
    const nodeRows = (node: number) => (i: number) => passes(i) &&
      (sourceNodes[source.codes[i]] === node || targetNodes[target.codes[i]] === node);
    const linkRows = (link: any) => (i: number) => passes(i) &&
      sourceNodes[source.codes[i]] === link.source.node && targetNodes[target.codes[i]] === link.target.node &&
      (linkColors === null || linkColors.codes[i] === link.color);
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
        ui.tooltip.showRowGroup(this.dataFrame, nodeRows(d.node), event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      //@ts-ignore
      .call(drag().subject((d) => d).on('drag', dragmove))
      .on('click', (event, d: any) => {
        if (event.defaultPrevented) return; // dragging
        this.dataFrame.selection.handleClick(nodeRows(d.node), event);
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
        ui.tooltip.showRowGroup(this.dataFrame, linkRows(d), event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      .on('click', (event, d: any) => {
        this.dataFrame.selection.handleClick(linkRows(d), event);
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
