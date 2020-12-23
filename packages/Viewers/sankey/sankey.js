import {drag} from 'd3-drag';
import {scaleOrdinal} from 'd3-scale';
import {select} from 'd3-selection';
import {sankey, sankeyLinkHorizontal} from 'd3-sankey';

export class SankeyViewer extends DG.JsViewer {
  constructor() {
    super();

    // Properties
    this.sourceColumnName = this.string('sourceColumnName');
    this.targetColumnName = this.string('targetColumnName');
    this.valueColumnName = this.float('valueColumnName');

    this.initialized = false;
  }

  init() {
    // Data
    this.graph = {};
    // Chart Settings
    this.margin = {top: 10, right: 10, bottom: 10, left: 10};
    this.color = scaleOrdinal(DG.Color.categoricalPalette);
    this.nodeWidth = 10;
    this.nodePadding = 15;
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

    this.sourceCol = this.dataFrame.getCol(this.sourceColumnName);
    this.targetCol = this.dataFrame.getCol(this.targetColumnName);
    this.valueCol = this.dataFrame.getCol(this.valueColumnName);
    let sourceCats = this.sourceCol.categories;
    let targetCats = this.targetCol.categories;
    let nodes = Array.from(new Set(sourceCats.concat(targetCats)))
      .map((node, index) => ({node: index, name: node}));

    let links = [];
    let rowCount = this.dataFrame.rowCount;
    let source = this.sourceCol.getRawData();
    let target = this.targetCol.getRawData();
    let value = this.valueCol.getRawData();
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

  onPropertyChanged(property) {
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
    $(this.root).empty();
    let width = this.root.parentElement.clientWidth - this.margin.left - this.margin.right;
    let height = this.root.parentElement.clientHeight - this.margin.top - this.margin.bottom;

    let generator = sankey().nodeWidth(this.nodeWidth)
      .nodePadding(this.nodePadding)
      .extent([[0, 0], [width, height]]);
    let graph = generator(this.graph);

    let svg = select(this.root).append("svg")
        .attr("width", width + this.margin.left + this.margin.right)
        .attr("height", height + this.margin.top + this.margin.bottom)
      .append("g")
        .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);

    let nodes = svg.append('g')
      .selectAll("rect")
      .data(graph.nodes)
      .join("rect")
        .attr("x", d => d.x0)
        .attr("y", d => d.y0)
        .attr("height", d => d.y1 - d.y0)
        .attr("width", d => d.x1 - d.x0)
        .attr("fill", d => DG.Color.toRgb(this.color(d.name)))
      .on("mouseover", (event, d) => {
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.sourceCol.get(i) === d.name ||
            this.targetCol.get(i) === d.name;
        }, event.x, event.y);
      })
      .on("mouseout", () => ui.tooltip.hide())
      // .call(drag().subject(d => d)
      //   .on("start", () => select(this).attr("stroke", "#fff"))
      //   .on("drag", dragmove)
      //   .on("end", () => select(this).attr("stroke", null)))
      .on("click", (event, d) => {
        if (event.defaultPrevented) return;  // dragging
        this.dataFrame.selection.handleClick(i => {
          return this.sourceCol.get(i) === d.name ||
            this.targetCol.get(i) === d.name;
        }, event);
      });

    let links = svg.append("g")
      .selectAll("path")
      .data(graph.links)
      .join("path")
        .attr("class", "link")
        .attr("d", sankeyLinkHorizontal())
        .attr("stroke-width", d => Math.max(1, d.width));

    let titles = svg.append("g")
        .attr("font-family", "'Roboto', 'Roboto Local', sans-serif")
        .attr("font-size", 13)
      .selectAll("text")
      .data(graph.nodes)
      .join("text")
        .attr("x", d => d.x0 < width / 2 ? d.x1 + 6 : d.x0 - 6)
        .attr("y", d => (d.y1 + d.y0) / 2)
        .attr("dy", "0.35em")
        .attr("text-anchor", d => d.x0 < width / 2 ? "start" : "end")
        .text(d => d.name);

    function dragmove(event, d) {
      d.y0 = Math.max(0, Math.min(height - d.value, event.y));
      select(this).attr("transform", `translate(${d.x0}, ${d.y0})`);
      generator.update(graph);
      links.attr("d", sankeyLinkHorizontal());
    };
  }
}
