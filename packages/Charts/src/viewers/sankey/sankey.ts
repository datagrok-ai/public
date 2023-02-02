import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';

import {drag} from 'd3-drag';
import {ScaleOrdinal, scaleOrdinal} from 'd3-scale';
import {select} from 'd3-selection';

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

interface Node {
  node: number,
  name: string,
}

interface Link {
  source: number,
  target: number,
  value: number,
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

  sourceCol?: DG.Column;
  targetCol?: DG.Column;

  constructor() {
    super();
    // Properties
    this.sourceColumnName = this.string('sourceColumnName');
    this.targetColumnName = this.string('targetColumnName');
    this.valueColumnName = this.string('valueColumnName');

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

  onTableAttached() {
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));

    this.init();

    const columns = this.dataFrame.columns.toList();
    const strColumns = columns.filter((col) => col.type === 'string');
    const numColumns = columns.filter((col) => ['double', 'int'].includes(col.type));
    this.sourceColumnName = strColumns[0].name;
    this.targetColumnName = strColumns[1].name;
    this.valueColumnName = numColumns[0].name;

    this.prepareData();
    this.render();
  }

  prepareData() {
    const dataFrameSourceColumn = this.dataFrame.getCol(this.sourceColumnName);
    const dataFrameTargetColumn = this.dataFrame.getCol(this.targetColumnName);
    const dataFrameValueColumn = this.dataFrame.getCol(this.valueColumnName);
    const filteredIndexList = this.dataFrame.filter.getSelectedIndexes();

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

    const nodes: Node[] = Array.from(new Set(sourceList.concat(targetList)))
      .map((node: string, index: number) => ({node: index, name: node}));

    const links: Link[] = [];
    const rowCount = filteredIndexList.length;
    for (let i = 0; i < rowCount; i++) {
      links.push({
        source: nodes.findIndex((node) => node.name === sourceList[i]),
        target: nodes.findIndex((node) => node.name === targetList[i]),
        value: valueList[i],
      });
    }

    this.graph = {
      nodes: nodes,
      links: links,
    };
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

  render() {
    this.prepareData();

    $(this.root).empty();
    const width = this.root.parentElement!.clientWidth - this.margin!.left - this.margin!.right;
    const height = this.root.parentElement!.clientHeight - this.margin!.top - this.margin!.bottom;

    const generator = sankey().nodeWidth(this.nodeWidth!)
      .nodePadding(this.nodePadding!)
      .nodeAlign(this.alignMethod!.justify)
      .extent([[0, 0], [width, height]]);
    const graph = generator(this.graph!);

    const svg = select(this.root).append('svg')
      .attr('width', width + this.margin!.left + this.margin!.right)
      .attr('height', height + this.margin!.top + this.margin!.bottom)
      .append('g')
      .attr('transform', `translate(${this.margin!.left}, ${this.margin!.top})`);

    const nodeGroup = svg.append('g').attr('class', 'node');

    const dataFrameSourceColumn = this.dataFrame.getCol('source');
    const dataFrameTargetColumn = this.dataFrame.getCol('target');

    const nodes = nodeGroup
      .selectAll('rect')
      .data(graph.nodes)
      .join('rect')
      .attr('x', (d) => d.x0!)
      .attr('y', (d) => d.y0!)
      .attr('height', (d) => d.y1! - d.y0!)
      .attr('width', (d) => d.x1! - d.x0!)
      .attr('fill', (d: any) => DG.Color.toRgb(this.color!(d.name)))
      .on('mouseover', (event, d: any) => {
        ui.tooltip.showRowGroup(this.dataFrame, (i) => {
          return this.sourceCol!.get(i) === d.name ||
            this.targetCol!.get(i) === d.name;
        }, event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      //@ts-ignore
      .call(drag().subject((d) => d).on('drag', dragmove))
      .on('click', (event, d: any) => {
        if (event.defaultPrevented) return; // dragging
        this.dataFrame.selection.handleClick((i) => {
          return dataFrameSourceColumn.get(i) === d.name ||
            dataFrameTargetColumn.get(i) === d.name;
        }, event);
      });

    const links = svg.append('g')
      .selectAll('path')
      .data(graph.links)
      .join('path')
      .attr('class', 'link')
      .attr('d', sankeyLinkHorizontal())
      .attr('stroke-width', (d) => Math.max(1, d.width!))
      .on('mouseover', (event, d: any) => {
        ui.tooltip.showRowGroup(this.dataFrame, (i) => {
          return this.sourceCol!.get(i) === d.source.name &&
            this.targetCol!.get(i) === d.target.name;
        }, event.x, event.y);
      })
      .on('mouseout', () => ui.tooltip.hide())
      .on('click', (event, d: any) => {
        this.dataFrame.selection.handleClick((i) => {
          return dataFrameSourceColumn.get(i) === d.source.name &&
            dataFrameTargetColumn.get(i) === d.target.name;
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
