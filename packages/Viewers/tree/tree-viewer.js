import * as d3 from 'd3';

export class TreeViewer extends DG.JsViewer {
  constructor() {
    super();
    this.legendDepth = this.int('legendDepth', 3);
    this.innerRadiusMargin = this.int('innerRadiusMargin', 100);
    this.labelFontSize = this.string('labelFontSize', '6px');
    this.legendFontSize = this.string('legendFontSize', '18px');
    this.showLabels = this.bool('showLabels', false);
    this.fixedDistance = '0.5';
    this.defaultSize = 400;
  }

  getLegendDomain(data, depth) {
    if (!Array.isArray(data)) data = [data];
    let domain = [];
  
    while (depth) {
      let subtrees = [];
      for (let child of data) {
        if (child.name && child.name !== 'root') domain.push(child.name);
        if (child.branchset) subtrees = [...subtrees, ...child.branchset];
      }
      data = subtrees;
      depth -= 1;
    }

    return domain;
  }

  onTableAttached() {
    this.newick = this.dataFrame.getTag('.newick');
    this.parsedNewick = JSON.parse(this.dataFrame.getTag('.newickJson')
      .replace(/children/g, 'branchset').replace(/attribute/g, 'length'), (key, value) => {
        return key === 'length' ? parseFloat(value || this.fixedDistance) : value;
      }).json;
    this.render();
  }

  onPropertyChanged() { this.render() }

  render() {
    $(this.root).empty();

    if (this.newick == null) {
      this.root.appendChild(ui.divText('Newick tag not found.', 'd4-viewer-error'));
      return;
    }

    const domain = this.getLegendDomain(this.parsedNewick, this.legendDepth);

    /** https://observablehq.com/@mbostock/tree-of-life */

    let width = this.root.clientWidth;
    if (!width) {
      width = this.defaultSize;
      this.showLabels = false;
    }
    let outerRadius = width / 2;
    let innerRadius = outerRadius - (this.showLabels ? this.innerRadiusMargin : 0);

    let cluster = d3.cluster()
      .size([360, innerRadius])
      .separation((a, b) => 1);

    let color = d3.scaleOrdinal()
      .domain(domain)
      .range(d3.schemeCategory10);

    function maxLength(d) {
      return d.data.length + (d.children ? d3.max(d.children, maxLength) : 0);
    }

    function setRadius(d, y0, k) {
      d.radius = (y0 += d.data.length) * k;
      if (d.children) d.children.forEach(d => setRadius(d, y0, k));
    }

    function setColor(d) {
      var name = d.data.name;
      d.color = color.domain().indexOf(name) >= 0 ? color(name) : d.parent ? d.parent.color : null;
      if (d.children) d.children.forEach(setColor);
    }

    function linkVariable(d) {
      return linkStep(d.source.x, d.source.radius, d.target.x, d.target.radius);
    }

    function linkConstant(d) {
      return linkStep(d.source.x, d.source.y, d.target.x, d.target.y);
    }

    function linkExtensionVariable(d) {
      return linkStep(d.target.x, d.target.radius, d.target.x, innerRadius);
    }

    function linkExtensionConstant(d) {
      return linkStep(d.target.x, d.target.y, d.target.x, innerRadius);
    }

    function linkStep(startAngle, startRadius, endAngle, endRadius) {
      const c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI);
      const s0 = Math.sin(startAngle);
      const c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI);
      const s1 = Math.sin(endAngle);
      return "M" + startRadius * c0 + "," + startRadius * s0
          + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
          + "L" + endRadius * c1 + "," + endRadius * s1;
    }

    let legend = svg => {
      const g = svg
        .selectAll("g")
        .data(color.domain())
        .join("g")
          .attr("transform", (d, i) => `translate(${this.showLabels ?  (-outerRadius - 50) : -width},${-outerRadius + i * 20})`);

      g.append("rect")
          .attr("width", 18)
          .attr("height", 18)
          .attr("fill", color);

      g.append("text")
          .attr("x", 24)
          .attr("y", 9)
          .attr("dy", "0.35em")
          .attr("font-size", Number.parseInt(this.legendFontSize))
          .text(d => d);
    }

    const root = d3.hierarchy(this.parsedNewick, d => d.branchset)
        .sum(d => d.branchset ? 0 : 1)
        .sort((a, b) => (a.value - b.value) || d3.ascending(a.data.length, b.data.length));


    cluster(root);
    setRadius(root, root.data.length = 0, innerRadius / maxLength(root));
    setColor(root);

    const svg = d3.select(this.root).append("svg")
        .attr("viewBox", [-outerRadius, -outerRadius, width, width])
        .attr("font-family", "sans-serif")
        //.attr("font-size", Number.parseInt(this.fontSize));

    svg.append("g")
        .call(legend);

    svg.append("style").text(`

    .link--active {
      stroke: #000 !important;
      stroke-width: 1.5px;
    }

    .link-extension--active {
      stroke-opacity: .6;
    }

    .label--active {
      font-weight: bold;
    }

    `);

    const linkExtension = svg.append("g")
        .attr("fill", "none")
        .attr("stroke", "#000")
        .attr("stroke-opacity", 0.25)
      .selectAll("path")
      .data(root.links().filter(d => !d.target.children))
      .join("path")
        .each(function(d) { d.target.linkExtensionNode = this; })
        .attr("d", linkExtensionVariable);

    const link = svg.append("g")
        .attr("fill", "none")
        .attr("stroke", "#000")
      .selectAll("path")
      .data(root.links())
      .join("path")
        .each(function(d) { d.target.linkNode = this; })
        .attr("d", linkVariable)
        .attr("stroke", d => d.target.color);

    if (this.showLabels) {
      svg.append("g")
      .selectAll("text")
      .data(root.leaves())
      .join("text")
        .attr("font-size", Number.parseInt(this.labelFontSize))
        .attr("dy", ".31em")
        .attr("transform", d => `rotate(${d.x - 90}) translate(${innerRadius + 4},0)${d.x < 180 ? "" : " rotate(180)"}`)
        .attr("text-anchor", d => d.x < 180 ? "start" : "end")
        .text(d => d.data.name.replace(/_/g, " "))
        .on("mouseover", mouseovered(true))
        .on("mouseout", mouseovered(false));
    }

    function mouseovered(active) {
      return function(event, d) {
        d3.select(this).classed("label--active", active);
        d3.select(d.linkExtensionNode).classed("link-extension--active", active).raise();
        do d3.select(d.linkNode).classed("link--active", active).raise();
        while (d = d.parent);
      };
    }

    this.root.style = 'position: absolute; left: 0; right: 0; top: 0; bottom: 0;';
  }
}
