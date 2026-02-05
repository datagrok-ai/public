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

import '../../../css/sankey-viewer.css';
import { MessageHandler } from '../../utils/utils';

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

  sourceCol: DG.Column | null = null;
  targetCol: DG.Column | null = null;

  nodeColorColumnName: string;
  linkColorColumnName: string;

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

  _testColumns() {
    const columns = this.dataFrame.columns.toList();
    const strColumns = columns.filter((col) => col.type === 'string' && col.categories.length <= 50);
    const numColumns = columns.filter((col) => ['double', 'int'].includes(col.type));
    return (strColumns!.length >= 2 && numColumns!.length >= 1);
  }

  onTableAttached() {
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.onMetadataChanged.subscribe((_) => this.render()));

    this.init();

    if (this._testColumns()) {
      const columns = this.dataFrame.columns.toList();
      const strColumns = columns.filter((col) => col.type === 'string' && col.categories.length <= 50);
      const numColumns = columns.filter((col) => ['double', 'int'].includes(col.type));

      this.sourceColumnName = strColumns[0].name;
      this.targetColumnName = strColumns[1].name;
      this.valueColumnName = numColumns[0].name;
    }

    this.prepareData();
    this.render();
  }

  onSourceRowsChanged() {
    this.render();
  }

  prepareData() {
    if (!this._testColumns())
      return;

    const dataFrameSourceColumn = this.dataFrame.getCol(this.sourceColumnName);
    const dataFrameTargetColumn = this.dataFrame.getCol(this.targetColumnName);
    const dataFrameValueColumn = this.dataFrame.getCol(this.valueColumnName);
    const selectedIndexes = this.filter.getSelectedIndexes();
    const { rowCount } = this.dataFrame;
    const filteredIndexList = selectedIndexes.length > 0 ?
      selectedIndexes :
      Array.from({ length: rowCount }, (_, index) => index);

    const sourceList = new Array<string>(filteredIndexList.length);
    const targetList = new Array<string>(filteredIndexList.length);
    const valueList = new Array<number>(filteredIndexList.length);

    for (let i = 0; i < filteredIndexList.length; i++) {
      sourceList[i] = dataFrameSourceColumn.get(filteredIndexList[i]);
      targetList[i] = dataFrameTargetColumn.get(filteredIndexList[i]);
      valueList[i] = dataFrameValueColumn.get(filteredIndexList[i]);
    }

    this.sourceCol = DG.Column.fromList('string', this.sourceColumnName, sourceList);
    this.targetCol = DG.Column.fromList('string', this.targetColumnName, targetList);

    const nodeRowIdx = new Map<string, number>();
    const filteredRowCount = filteredIndexList.length;
    for (let i = 0; i < filteredRowCount; i++) {
      const rowIdx = filteredIndexList[i];
      if (!nodeRowIdx.has(sourceList[i]))
        nodeRowIdx.set(sourceList[i], rowIdx);
      if (!nodeRowIdx.has(targetList[i]))
        nodeRowIdx.set(targetList[i], rowIdx);
    }

    const nodes: Node[] = Array.from(new Set(sourceList.concat(targetList)))
      .filter((element) => element)
      .map((name: string, index: number) => ({
        node: index,
        name: name,
        rowIdx: nodeRowIdx.get(name)!,
      }));

    const links: Link[] = [];
    for (let i = 0; i < filteredRowCount; i++) {
      const source = nodes.findIndex((node) => node.name === sourceList[i]);
      const target = nodes.findIndex((node) => node.name === targetList[i]);
      if (source === -1 || target === -1)
        continue;

      links.push({
        source: source,
        target: target,
        value: valueList[i],
        rowIdx: filteredIndexList[i],
      });
    }

    this.graph = {
      nodes: nodes,
      links: links,
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
    if (this.initialized) {
      this.prepareData();
      this.render();
    }
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  rowMatchesColumnValues(sourceCol: DG.Column | null, targetCol: DG.Column | null,
    rowIndex: number, sourceName: string, targetName: string, operator: 'AND' | 'OR' = 'AND',
  ): boolean {
    if (!sourceCol || !targetCol)
      return false;

    const sourceMatch = this.dataFrame.get(sourceCol.name, rowIndex) === sourceName;
    const targetMatch = this.dataFrame.get(targetCol.name, rowIndex) === targetName;

    switch (operator) {
    case 'AND': return sourceMatch && targetMatch;
    case 'OR': return sourceMatch || targetMatch;
    }
  }

  render() {
    $(this.root).empty();
    if (!this._testColumns()) {
      MessageHandler._showMessage(this.root, 'The Sankey viewer requires a minimum of 2 categorical (less than 50 unique categories) and 1 numerical columns.', 'd4-viewer-error');
      return;
    }

    this.prepareData();

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
          return this.rowMatchesColumnValues(this.sourceCol, this.targetCol, i, d.name, d.name, 'OR');
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
          return this.rowMatchesColumnValues(this.sourceCol, this.targetCol, i, d.source.name, d.target.name);
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
      .attr('class', (d: any) => (d.name).split(' ').join('_'))
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
      nodeGroup.select(`text.${d.name.split(' ').join('_')}`)
        .attr('transform', `translate(${d.x0 - rectX}, ${d.y0 - rectY})`);

      generator.update(graph);
      links.attr('d', sankeyLinkHorizontal());
    };
  }
}
