import {
  select,
  scaleLinear,
  scalePoint,
  axisLeft,
  axisBottom
} from 'd3';

export class DumbbellViewer extends DG.JsViewer {
  constructor() {
    super();
    this.subjectColumnName = this.string('subjectColumnName', 'USUBJID');
    this.startColumnName = this.int('startColumnName', 'AESTDY');
    this.endColumnName = this.int('startColumnName', 'AEENDY');
    this.initialized = false;
  }

  init() {
    this.margin = {top: 10, right: 10, bottom: 30, left: 90};
    this.xScale = scaleLinear();
    this.yScale = scalePoint();
    this.data = [];
    this.initialized = true;
  }

  onTableAttached() {
    this.init();

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

    this.render();
  }

  detach() {
    this.subs.forEach(sub => sub.unsubscribe());
  }

  render(computeData = true) {

    if (computeData) {
      this.data.length = 0;
      let subjectCol = this.dataFrame.getCol(this.subjectColumnName);
      let startCol = this.dataFrame.getCol(this.startColumnName).getRawData();
      let endCol = this.dataFrame.getCol(this.endColumnName);
      //this.xMax = endCol.max;
      //this.xMin = endCol.min;
      endCol = endCol.getRawData();
      this.xMax = 0;
      for (let i of this.dataFrame.filter.getSelectedIndexes()) {
        if (this.xMax < endCol[i]) this.xMax = endCol[i];
      }
      for (let i of this.dataFrame.filter.getSelectedIndexes()) {
        let isNull = endCol[i] === -2147483648;
        let curSubj = subjectCol.get(i);
        if (this.data.some(obj => obj.usubjid === curSubj)) curSubj += ' (1)';
        this.data.push({
          usubjid: curSubj,
          aestdy: startCol[i],
          imputed: isNull,
          aeendy: isNull ? this.xMax : endCol[i]
        });
      }
    }
    console.log(this.data);

    let width = this.root.parentElement.clientWidth;
    let height = this.root.parentElement.clientHeight;
    let innerWidth = width - this.margin.left - this.margin.right;
    let innerHeight = height - this.margin.top - this.margin.bottom;

    $(this.root).empty();
    this.root.style = 'width: 100%; height: 100%;';
    this.root.id = 'chart';

    let svg = select("#chart").append("svg")
      .attr("width", width)
      .attr("height", height);
    let g = svg.append("g")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);

    this.xScale
      .domain([0, this.xMax])
      .rangeRound([0, innerWidth])
      .nice();
    this.yScale
      .domain(this.data.map(d => d.usubjid))
      .rangeRound([this.margin.top, innerHeight])
      .padding(0.5);

    let yAxis = g.append("g").call(axisLeft(this.yScale));
    let xAxis = g.append("g").call(axisBottom(this.xScale))
      .attr("transform", `translate(0, ${innerHeight})`);
    let dumbbell = g.append("g");

    dumbbell.selectAll("line")
      .data(this.data)
      .join("line")
        .attr("stroke", "#808080")
        .attr("stroke-width", 2)
        .attr("x1", d => this.xScale(d.aestdy))
        .attr("x2", d => this.xScale(d.aeendy))
        .attr("y1", d => this.yScale(d.usubjid))
        .attr("y2", d => this.yScale(d.usubjid));

    dumbbell.selectAll("circle")
      .data(this.data)
      .join("circle")
        .attr("fill", "#808080")
        .attr("cx", d => this.xScale(d.aestdy))
        .attr("cy", d => this.yScale(d.usubjid))
        .attr("r", 5)
        .clone(true)
          .attr("cx", d => this.xScale(d.aeendy))
          .attr("fill", d => d.imputed ? "#fc8d62" : "#62d1fc");
  }
}
