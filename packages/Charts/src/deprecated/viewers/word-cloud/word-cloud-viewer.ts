import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as d3 from 'd3';

import cloud from 'd3-cloud';
import $ from 'cash-dom';


export class WordCloudViewer extends DG.JsViewer {
  strColumnName: string;
  strColumns: any;
  initialized: boolean;

  freqMap?: any;

  constructor() {
    super();

    // Properties
    this.strColumnName = this.string('strColumnName');

    this.strColumns = [];
    this.initialized = false;
  }

  init() {
    this.initialized = true;
  }

  testColumns() {
    return (this.strColumns.length >= 1);
  }

  onTableAttached() {
    this.init();

    const columns = this.dataFrame.columns.toList();
    this.strColumns = columns.filter((col) => col.type === 'string');

    if (this.testColumns()) {
      // Find a string column with the smallest number of unique values
      this.strColumnName = this.strColumns.reduce((prev: DG.Column, curr: DG.Column) =>
        prev.categories.length < curr.categories.length ? prev : curr).name;
    }

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));

    this.render();
  }

  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    if (this.initialized && this.testColumns()) {
      if (property.name === 'strColumnName') this.strColumnName = property.get(this);
      this.render();
    }
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  render() {
    if (!this.testColumns()) {
      this.root.innerText = 'Not enough data to produce the result.';
      return;
    }

    $(this.root).empty();
    const margin = {top: 10, right: 10, bottom: 10, left: 10};
    const width = this.root.parentElement!.clientWidth - margin.left - margin.right;
    const height = this.root.parentElement!.clientHeight - margin.top - margin.bottom;

    const svg = d3.select(this.root).append('svg')
      .attr('width', width + margin.left + margin.right)
      .attr('height', height + margin.top + margin.bottom);

    svg
      .append('g')
      .attr('transform', `translate(${margin.left}, ${margin.top})`);

    const strColumn = this.dataFrame.getCol(this.strColumnName);
    const words = strColumn.categories;
    this.freqMap = {};
    strColumn.toList().forEach((w) => this.freqMap[w] = (this.freqMap[w] || 0) + 1);
    const sortedWords = Object.keys(this.freqMap).sort((a, b) => {
      return this.freqMap[b] - this.freqMap[a];
    });
    const fontSize = d3.scaleLinear()
      .domain([0, this.freqMap[sortedWords[0]]])
      .range([10, 100]);

    const fontColor = () => {
      const nColors = DG.Color.categoricalPalette.length;
      const color = DG.Color.getCategoricalColor(Math.floor(Math.random() * nColors));
      return DG.Color.toRgb(color);
    };

    const layout = cloud()
      .size([width, height])
      .words(words.map((d) => {
        return {
          text: d,
          size: fontSize(this.freqMap[d]),
          color: fontColor(),
        };
      }))
      .padding(10)
      .rotate(() => (~~(Math.random() * 6) - 3) * 30)
      .fontSize((d: any) => d.size)
      .on('end', draw);

    layout.start();
    const table = this.dataFrame;

    function draw(words: string[]) {
      svg
        .append('g')
        .attr('transform', `translate(${layout.size()[0] / 2}, ${layout.size()[1] / 2})`)
        .selectAll('text')
        .data(words)
        .enter().append('text')
        .style('font-size', (d: any) => d.size)
        .style('fill', (d: any) => d.color)
        .attr('text-anchor', 'middle')
        .attr('transform', (d: any) => `translate(${[d.x, d.y]}) rotate(${d.rotate})`)
        .text((d: any) => d.text)
        .on('mouseover', (d) => ui.tooltip.showRowGroup(table, (i) => {
          return d.srcElement.innerHTML === strColumn.get(i);
        }, d.x, d.y))
        .on('mouseout', () => ui.tooltip.hide())
        .on('mousedown', (d) => {
          table.selection.handleClick((i) => {
            return d.srcElement.innerHTML === strColumn.get(i);
          }, d);
        });
    }

    this.root.appendChild(svg.node()!);
  }
}
