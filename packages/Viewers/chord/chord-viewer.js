import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import {select, scaleOrdinal, color} from 'd3';
import {layoutConf, topSort} from './configuration.js';


export class ChordViewer extends DG.JsViewer {

  constructor() {
    super();

    // Properties
    this.fromColumnName = this.string('fromColumnName');
    this.toColumnName = this.string('toColumnName');
    this.aggType = this.string('aggType', 'count', { choices: ['count', 'sum'] });
    this.chordLengthColumnName = this.float('chordLengthColumnName');
    this.colorBy = this.string('colorBy', 'source', { choices: ['source', 'target'] });
    this.sortBy = this.string('sortBy', 'topology', { choices: ['alphabet', 'frequency', 'topology'] });
    this.direction = this.string('direction', 'clockwise', { choices: ['clockwise', 'counterclockwise'] });

    this.initialized = false;
    this.numColumns = [];
    this.strColumns = [];
    this.fromCol;
    this.toCol;
    this.data = [];
    this.conf = layoutConf;
    this.chords = [];
    this.chordConf = {};
    this.segments = {};
  }

  init() {
    this.innerRadiusMargin = 80;
    this.outerRadiusMargin = 60;
    this.minSegmentWidth = 10;
    this.colorScale = scaleOrdinal(DG.Color.categoricalPalette);
    this.color = c => DG.Color.toRgb(this.colorScale(c));
    this.chordConf.color = datum => this.color(datum[this.colorBy]['label']);
    this.chordConf.opacity = 0.7;
    this.labelConf = {
      innerRadius: 1.02,
      style: { 'font-size': 12, fill: '#7f7f7f' }
    };
    this.initialized = true;
  }

  testColumns() {
    return (this.strColumns.length >= 2 && this.numColumns.length >= 1);
  }

  onTableAttached() {
    this.init();

    this.strColumns = this.dataFrame.columns.toList()
      .filter(col => col.type === 'string')
      .sort((a, b) => a.categories.length - b.categories.length);
    this.numColumns = [...this.dataFrame.columns.numerical];

    // TODO: Choose the most relevant columns
    if (this.testColumns()) {
      this.fromColumnName = this.strColumns[0].name;
      this.toColumnName = this.strColumns[1].name;
      this.chordLengthColumnName = this.numColumns[0].name;
    }

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

    this.render();
  }

  onPropertyChanged(property) {
    super.onPropertyChanged(property);
    if (this.initialized && this.testColumns()) this.render();
  }

  detach() {
    this.subs.forEach(sub => sub.unsubscribe());
  }

  generateData() {
    this.data.length = 0;

    this.fromColumn = this.dataFrame.getCol(this.fromColumnName);
    this.toColumn = this.dataFrame.getCol(this.toColumnName);
    this.conf.events = {
      mouseover: (datum, index, nodes, event) => {
        select(nodes[index]).select(`#${datum.id}`).attr('stroke', color(this.color(datum.label)).darker());
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.fromColumn.get(i) === datum.label ||
            this.toColumn.get(i) === datum.label;
        }, event.x, event.y);
      },
      mouseout: (datum, index, nodes, event) => {
        select(nodes[index]).select(`#${datum.id}`).attr('stroke', 'none');
        ui.tooltip.hide()
      },
      mousedown: (datum, index, nodes, event) => {
        this.dataFrame.selection.handleClick(i => {
          return this.fromColumn.get(i) === datum.label ||
            this.toColumn.get(i) === datum.label;
        }, event);
      }
    };

    this.freqMap = {};
    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      let from = this.fromColumn.isNone(i) ? "" : this.fromColumn.get(i);
      this.freqMap[from] = (this.freqMap[from] || 0) + 1;
      let to = this.toColumn.isNone(i) ? "" : this.toColumn.get(i);
      this.freqMap[to] = (this.freqMap[to] || 0) + 1;
    }

    if (this.aggType === 'count') this.chordLengthColumnName = null;
    if (this.fromColumnName !== this.toColumnName) {
      this.aggregatedTable = this.dataFrame
        .groupBy([this.fromColumnName, this.toColumnName])
        .add(this.aggType, this.chordLengthColumnName, 'result')
        .aggregate();

      this.fromCol = this.aggregatedTable.getCol(this.fromColumnName);
      this.toCol = this.aggregatedTable.getCol(this.toColumnName);
      this.aggVal = this.aggregatedTable.getCol('result').getRawData();
      this.rowCount = this.aggregatedTable.rowCount;
  
      this.categories = Array.from(new Set(this.fromCol.categories.concat(this.toCol.categories)));
    } else {
      this.categories = Array.from(this.fromColumn.categories);
    }

    this.data = this.categories
      .sort((this.sortBy === 'frequency') ? (a, b) => this.freqMap[b] - this.freqMap[a] : undefined)
      .map((s, ind) => {
        return {
          id: `id-${ind}`,
          label: s,
          len: this.freqMap[s],
          color: this.color(s)
        }
    });

    if (this.fromColumnName !== this.toColumnName) {
      for (const prop of Object.getOwnPropertyNames(this.segments)) {
        delete this.segments[prop];
      }
  
      this.data.forEach(s => {
        this.segments[s.label] = { datum: s, targets: [], aggTotal: null, visited: false };
        s.pos = 0;
      });
  
      for (let i = 0; i < this.rowCount; i++) {
        let from = this.fromCol.get(i);
        let to = this.toCol.get(i);
  
        this.segments[from]['targets'].push(to);
        this.segments[from]['aggTotal'] = (this.segments[from]['aggTotal'] || 0) + this.aggVal[i];
        if (from !== to) this.segments[to]['aggTotal'] = (this.segments[to]['aggTotal'] || 0) + this.aggVal[i];
      }
    }

    if (this.sortBy === 'topology') {
      if (this.fromColumnName === this.toColumnName) {
        this.props.sortBy = 'alphabet';
      } else {
        this.data = topSort(this.segments);
      }
    }
    if (this.direction === 'counterclockwise') this.data.reverse();
  }

  computeChords() {
    this.chords.length = 0;
    let source = this.fromCol.getRawData();
    let fromCatList = this.fromCol.categories;
    let target = this.toCol.getRawData();
    let toCatList = this.toCol.categories;

    for (let i = 0; i < this.rowCount; i++) {
      let sourceLabel = fromCatList[source[i]];
      let targetLabel = toCatList[target[i]];
      let sourceBlock = this.segments[sourceLabel]['datum'];
      let targetBlock = this.segments[targetLabel]['datum'];
      let sourceStep = sourceBlock.len * (this.aggVal[i] / this.segments[sourceLabel]['aggTotal']);
      let targetStep = targetBlock.len * (this.aggVal[i] / this.segments[targetLabel]['aggTotal']);

      this.chords.push({
        source: {
          id: sourceBlock.id,
          start: sourceBlock.pos,
          end: sourceBlock.pos + sourceStep,
          label: sourceLabel
        },
        target: {
          id: targetBlock.id,
          start: targetBlock.pos,
          end: targetBlock.pos + targetStep,
          label: targetLabel
        },
        value: this.aggVal[i],
      });

      sourceBlock.pos += sourceStep;
      if (sourceLabel === targetLabel) continue;
      targetBlock.pos += targetStep;
    }

    this.chordConf.events = {
      mouseover: (datum, index, nodes, event) => {
        select(nodes[index]).attr('opacity', 0.9);
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.fromColumn.get(i) === datum.source.label &&
            this.toColumn.get(i) === datum.target.label;
        }, event.x, event.y);
      },
      mouseout: (datum, index, nodes, event) => {
        select(nodes[index]).attr('opacity', 0.7);
        ui.tooltip.hide()
      },
      mousedown: (datum, index, nodes, event) => {
        this.dataFrame.selection.handleClick(i => {
          return this.fromColumn.get(i) === datum.source.label &&
            this.toColumn.get(i) === datum.target.label;
        }, event);
      }
    };

  }

  render(computeData = true) {

    if (!this.testColumns()) {
      this.root.innerText = 'Not enough data to produce the result.';
      return;
    }

    if (computeData) {
      this.generateData();
      if (this.fromColumnName !== this.toColumnName) this.computeChords();
    }

    $(this.root).empty();
    let width = this.root.parentElement.clientWidth;
    let height = this.root.parentElement.clientHeight;
    let size = Math.min(width, height);

    let circos = Circos({
      container: this.root,
      width: size,
      height: size
    });

    this.conf.innerRadius = Math.max(0, size/2 - this.innerRadiusMargin);
    this.conf.outerRadius = Math.max(0, size/2 - this.outerRadiusMargin);
    // this.chordConf.radius = d => (d.source.id === d.target.id) ? this.conf.outerRadius : null;

    circos.layout(this.data, this.conf);

    if (this.fromColumnName !== this.toColumnName) {
      circos.chords('chords-track', this.chords, this.chordConf);
    }

    circos.text('labels', this.data.map(d => { return {
      block_id: d.id,
      position: this.freqMap[d.label] / 2,
      value: d.label
    }}), this.labelConf);
    circos.render();

    let labels = select(this.root).selectAll('.block');

    // fix label rotation past 180
    labels.filter((d, i, nodes) => {
      return +(select(nodes[i]).attr('transform').match(/\d+\.?\d*/g)[0]) >= 180;
    }).selectAll('text')
        .attr('transform', (d, i, nodes) => select(nodes[i]).attr('transform') + ' rotate(180) ')
        .attr('text-anchor', 'end');

    // ellipsis
    labels.selectAll('text').each((d, i, nodes) => {
      let el = select(nodes[i]);
      let textLength = el.node().getComputedTextLength();
      let text = el.text();
      while (text.length && textLength > (this.outerRadiusMargin - 15)) {
        text = text.slice(0, -1);
        el.text(text + '\u2026');
        textLength = el.node().getComputedTextLength();
      }
    });

    this.root.firstChild.style = 'position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);';
  }
}
