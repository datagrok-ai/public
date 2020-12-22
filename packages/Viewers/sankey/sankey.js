import * as d3 from 'd3';
import {sankey, sankeyLinkHorizontal} from 'd3-sankey';

export class SankeyViewer extends DG.JsViewer {
  constructor() {
    super();

    // Properties
    this.sourceColumnName = this.string('sourceColumnName');
    this.targetColumnName = this.string('targetColumnName');
    this.valueColumnName = this.float('valueColumnName');
  }

  init() {
    // Data
    this.graph = {};
    // Chart Settings
    this.margin = {top: 10, right: 10, bottom: 10, left: 10};
    this.color = d3.scaleOrdinal(d3.schemeCategory10); // TODO: use DG.Color.categoricalPalette
    this.initialized = true;
  }

  onTableAttached() {
    this.init();

    let columns = this.dataFrame.columns.toList();
    this.strColumns = columns.filter(col => col.type === 'string');
    this.numColumns = columns.filter(col => ['double', 'int'].includes(col.type));

    this.sourceColumnName = this.strColumns[0].name;
    this.targetColumnName = this.strColumns[1].name;
    this.valueColumnName = this.numColumns[0].name;
    this.prepareData();

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.onSizeChanged, 50).subscribe((_) => this.render()));

    this.render();
  }

  prepareData() {

    let sourceCol = this.dataFrame.getCol(this.sourceColumnName);
    let targetCol = this.dataFrame.getCol(this.targetColumnName);
    let valueCol = this.dataFrame.getCol(this.valueColumnName);
    let sourceCats = sourceCol.categories;
    let targetCats = targetCol.categories;
    let nodes = Array.from(new Set(sourceCats.concat(targetCats)))
      .map((node, index) => ({node: index, name: node}));

    let links = [];
    let rowCount = this.dataFrame.rowCount;
    let source = sourceCol.getRawData();
    let target = targetCol.getRawData();
    let value = valueCol.getRawData();
    for (let i = 0; i < rowCount; i++) {
      links.push({
        source: nodes.findIndex(node => node.name === sourceCats[source[i]]),
        target: nodes.findIndex(node => node.name === targetCats[target[i]]),
        value: value[i]
      });
    }

    this.graph = {
      nodes: nodes,
      links: links
    };

  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  render() {
    $(this.root).empty();
    let width = this.root.parentElement.clientWidth - this.margin.left - this.margin.right;
    let height = this.root.parentElement.clientHeight - this.margin.top - this.margin.bottom;

    let {nodes, links} = sankey().extent([[0, 0], [width, height]])(this.graph);

    let svg = d3.select(this.root).append("svg")
        .attr("width", width + this.margin.left + this.margin.right)
        .attr("height", height + this.margin.top + this.margin.bottom)
      .append("g")
        .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);

    svg.append('g').attr("stroke", "#000")
      .selectAll("rect")
      .data(nodes)
      .join("rect")
        .attr("x", d => d.x0)
        .attr("y", d => d.y0)
        .attr("height", d => d.y1 - d.y0)
        .attr("width", d => d.x1 - d.x0)
        .attr("fill", d => this.color(d.name))
      .on('mouseover', d => ui.tooltip.showRowGroup(this.dataFrame, i => true, d.x, d.y))
      .on('mouseout', () => ui.tooltip.hide())
      .on('mousedown', d => this.dataFrame.selection.handleClick(i => true, d));

    svg.append("g")
        .attr("fill", "none")
        .attr("stroke", "#000")
        .attr("stroke-opacity", 0.2)
      .selectAll("path")
      .data(links)
      .join("path")
        .attr("d", sankeyLinkHorizontal())
        .attr("stroke-width", d => Math.max(1, d.width));

    svg.append("g")
        .attr("font-family", "sans-serif")
        .attr("font-size", 10)
      .selectAll("text")
      .data(nodes)
      .join("text")
        .attr("x", d => d.x0 < width / 2 ? d.x1 + 6 : d.x0 - 6)
        .attr("y", d => (d.y1 + d.y0) / 2)
        .attr("dy", "0.35em")
        .attr("text-anchor", d => d.x0 < width / 2 ? "start" : "end")
        .text(d => d.name);
  }
}
