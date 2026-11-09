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
  sourceCode: number,
  targetCode: number,
}

interface Link {
  source: number,
  target: number,
  value: number,
  rowIdx: number,
  sourceCode: number,
  targetCode: number,
  color: number,
}

interface Side {
  names: string[],
  codes: Int32Array,
  firstRows: Int32Array,
  nodeOf: Int32Array,
}

interface SankeyData {
  nodes: Node[],
  links: Link[],
  negative: boolean,
  sources: Int32Array,
  targets: Int32Array,
  colors: Int32Array | null,
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
      .filter((col) => col.type === 'string' && !col.isVirtual && col.categories.length <= 50).map((col) => col.name);
  }

  numColumnNames(): string[] {
    return this.dataFrame.columns.toList().filter((col) => ['double', 'int'].includes(col.type) && !col.isVirtual).map((col) => col.name);
  }

  hasColumns(): boolean {
    const plain = (name: string, types: string[]) => {
      const col = this.dataFrame.col(name);
      return col !== null && !col.isVirtual && types.includes(col.type);
    };
    return plain(this.sourceColumnName, [DG.COLUMN_TYPE.STRING]) && plain(this.targetColumnName, [DG.COLUMN_TYPE.STRING]) &&
      plain(this.valueColumnName, [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT]);
  }

  private side(name: string): Side {
    const col = this.dataFrame.getCol(name);
    const names = col.categories;
    return {names: names, codes: col.getRawData() as Int32Array,
      firstRows: new Int32Array(names.length).fill(-1), nodeOf: new Int32Array(names.length)};
  }

  private stringLinkColorCol(): DG.Column | null {
    const col = this.linkColorColumnName ? this.dataFrame.col(this.linkColorColumnName) : null;
    return col !== null && col.type === DG.COLUMN_TYPE.STRING && !col.isVirtual ? col : null;
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

  private prepareData(): SankeyData {
    const source = this.side(this.sourceColumnName);
    const target = this.side(this.targetColumnName);
    const valueCol = this.dataFrame.getCol(this.valueColumnName);
    const values = valueCol.getRawData() as Int32Array | Float32Array | Float64Array;
    // a missing value is the int or float null, which would sink both nodes of the row
    const missing = valueCol.type === DG.COLUMN_TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL;
    const colorCol = this.stringLinkColorCol();
    const colors = (colorCol?.getRawData() ?? null) as Int32Array | null;
    const colorCount = colorCol?.categories.length ?? 1;
    const links = new Map<number, Link>();
    let negative = false;
    for (const row of this.filter.getSelectedIndexes()) {
      const s = source.codes[row];
      const t = target.codes[row];
      if (source.firstRows[s] < 0)
        source.firstRows[s] = row;
      if (target.firstRows[t] < 0)
        target.firstRows[t] = row;
      if (!source.names[s] || !target.names[t])
        continue;

      const value = values[row] === missing ? 0 : values[row];
      negative ||= value < 0;
      const color = colors?.[row] ?? 0;
      const key = (s * target.names.length + t) * colorCount + color;
      const link = links.get(key);
      if (link)
        link.value += value;
      else
        links.set(key, {source: -1, target: -1, value: value, rowIdx: row, sourceCode: s, targetCode: t, color: color});
    }

    const nodes: Node[] = [];
    const byName = new Map<string, Node>();
    const addNodes = ({names, firstRows, nodeOf}: Side, code: 'sourceCode' | 'targetCode') => {
      const shown = Array.from(firstRows.keys())
        .filter((c) => firstRows[c] >= 0 && names[c])
        .sort((a, b) => firstRows[a] - firstRows[b]);
      for (const c of shown) {
        let node = byName.get(names[c]);
        if (node === undefined) {
          node = {node: nodes.length, name: names[c], rowIdx: firstRows[c], sourceCode: -1, targetCode: -1};
          nodes.push(node);
          byName.set(node.name, node);
        }
        node.rowIdx = Math.min(node.rowIdx, firstRows[c]);
        node[code] = c;
        nodeOf[c] = node.node;
      }
    };
    addNodes(source, 'sourceCode');
    addNodes(target, 'targetCode');
    for (const link of links.values()) {
      link.source = source.nodeOf[link.sourceCode];
      link.target = target.nodeOf[link.targetCode];
    }

    return {nodes: nodes, links: Array.from(links.values()), negative: negative,
      sources: source.codes, targets: target.codes, colors: colors};
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

    const {sources, targets, colors} = data;
    const filter = this.filter.getBuffer();
    const passes = (i: number) => (filter[i >>> 5] & (1 << (i & 31))) !== 0;
    const nodeRows = (node: Node) => (i: number) => passes(i) &&
      (sources[i] === node.sourceCode || targets[i] === node.targetCode);
    const linkRows = (link: Link) => (i: number) => passes(i) &&
      sources[i] === link.sourceCode && targets[i] === link.targetCode && (colors === null || colors[i] === link.color);
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
      .attr('fill', (d) => {
        const color = nodeColorCol ? DG.Color.getRowColor(nodeColorCol, d.rowIdx) : this.color!(d.name);
        return DG.Color.toRgb(color);
      })
      .on('mouseover', (event, d) => {
        ui.tooltip.showRowGroup(this.dataFrame, nodeRows(d), event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      //@ts-ignore
      .call(drag().subject((d) => d).on('drag', dragmove))
      .on('click', (event, d) => {
        if (event.defaultPrevented) return; // dragging
        this.dataFrame.selection.handleClick(nodeRows(d), event);
      });

    const links = svg.append('g')
      .selectAll('path')
      .data(graph.links)
      .join('path')
      .attr('class', 'link')
      .attr('d', sankeyLinkHorizontal())
      .style('stroke', (d) => {
        const color = linkColorCol ? DG.Color.getRowColor(linkColorCol, d.rowIdx) : DG.Color.darkGray;
        return DG.Color.toRgb(color);
      })
      .style('fill', 'none')
      .attr('stroke-width', (d) => Math.max(1, d.width!))
      .on('mouseover', (event, d) => {
        ui.tooltip.showRowGroup(this.dataFrame, linkRows(d), event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      .on('click', (event, d) => {
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
      .text((d) => d.name);

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
