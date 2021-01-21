import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import {select, scaleOrdinal} from 'd3';
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
    this.sortBy = this.string('sortBy', 'frequency', { choices: ['alphabet', 'frequency', 'topology'] });

    this.initialized = false;
    this.numColumns = [];
    this.strColumns = [];
    this.fromCol;
    this.toCol;
    this.data = [];
    this.conf = layoutConf;
    this.chords = [];
    this.chordConf = {};
    this.labels = [];
    this.segments = {};
  }

  init() {
    this.innerRadiusMargin = 60;
    this.outerRadiusMargin = 40;
    this.color = scaleOrdinal(DG.Color.categoricalPalette);
    this.chordConf.color = (datum, index) => DG.Color.toRgb(this.color(datum[this.colorBy]['id']));
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
    this.labels.length = 0;

    this.fromColumn = this.dataFrame.getCol(this.fromColumnName);
    this.toColumn = this.dataFrame.getCol(this.toColumnName);
    this.conf.events = {
      mouseover: (datum, index, nodes, event) => {
        select(nodes[index]).select(`#${datum.id}`).attr('stroke', 'black');
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.fromColumn.get(i) === datum.id ||
            this.toColumn.get(i) === datum.id;
        }, event.x, event.y);
      },
      mouseout: (datum, index, nodes, event) => {
        select(nodes[index]).select(`#${datum.id}`).attr('stroke', 'none');
        ui.tooltip.hide()
      },
      mousedown: (datum, index, nodes, event) => {
        this.dataFrame.selection.handleClick(i => {
          return this.fromColumn.get(i) === datum.id ||
            this.toColumn.get(i) === datum.id;
        }, event);
      }
    };

    if (this.aggType === 'count') this.chordLengthColumnName = null;
    this.aggregatedTable = this.dataFrame
      .groupBy([this.fromColumnName, this.toColumnName])
      .add(this.aggType, this.chordLengthColumnName, 'result')
      .aggregate();

    this.freqMap = {};
    this.fromColumn.toList().forEach(k => this.freqMap[k] = (this.freqMap[k] || 0) + 1);
    this.toColumn.toList().forEach(k => this.freqMap[k] = (this.freqMap[k] || 0) + 1);

    this.fromCol = this.aggregatedTable.getCol(this.fromColumnName);
    this.toCol = this.aggregatedTable.getCol(this.toColumnName);
    this.aggVal = this.aggregatedTable.getCol('result').getRawData();
    this.rowCount = this.aggregatedTable.rowCount;

    this.categories = Array.from(new Set(this.fromCol.categories.concat(this.toCol.categories)));
    this.data = this.categories
      .sort((this.sortBy === 'frequency') ? (a, b) => this.freqMap[b] - this.freqMap[a] : undefined)
      .map(s => {
        this.labels.push({ block_id: s, position: this.freqMap[s] / 2, value: s });
        return {
          id: s,
          len: this.freqMap[s],
          color: DG.Color.toRgb(this.color(s))
        }
    });

    for (const prop of Object.getOwnPropertyNames(this.segments)) {
      delete this.segments[prop];
    }

    this.data.forEach(s => {
      this.segments[s.id] = { datum: s, targets: [], aggTotal: null, visited: false };
      s.pos = 0;
    });

    for (let i = 0; i < this.rowCount; i++) {
      let from = this.fromCol.get(i);
      let to = this.toCol.get(i);

      this.segments[from]['targets'].push(to);
      this.segments[from]['aggTotal'] = (this.segments[from]['aggTotal'] || 0) + this.aggVal[i];
      if (from !== to) this.segments[to]['aggTotal'] = (this.segments[to]['aggTotal'] || 0) + this.aggVal[i];
    }

    if (this.sortBy === 'topology') this.data = topSort(this.segments);
  }

  computeChords() {
    this.chords.length = 0;
    let source = this.fromCol.getRawData();
    let fromCatList = this.fromCol.categories;
    let target = this.toCol.getRawData();
    let toCatList = this.toCol.categories;

    for (let i = 0; i < this.rowCount; i++) {
      let sourceId = fromCatList[source[i]];
      let targetId = toCatList[target[i]];
      let sourceBlock = this.segments[sourceId]['datum'];
      let targetBlock = this.segments[targetId]['datum'];
      let sourceStep = sourceBlock.len * (this.aggVal[i] / this.segments[sourceId]['aggTotal']);
      let targetStep = targetBlock.len * (this.aggVal[i] / this.segments[targetId]['aggTotal']);

      this.chords.push({
        source: {
          id: sourceId,
          start: sourceBlock.pos,
          end: sourceBlock.pos + sourceStep
        },
        target: {
          id: targetId,
          start: targetBlock.pos,
          end: targetBlock.pos + targetStep
        },
        value: this.aggVal[i]
      });

      sourceBlock.pos += sourceStep;
      if (sourceId === targetId) continue;
      targetBlock.pos += targetStep;
    }

    this.chordConf.events = {
      mouseover: (datum, index, nodes, event) => {
        select(nodes[index]).attr('opacity', 0.9);
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.fromColumn.get(i) === datum.source.id &&
            this.toColumn.get(i) === datum.target.id;
        }, event.x, event.y);
      },
      mouseout: (datum, index, nodes, event) => {
        select(nodes[index]).attr('opacity', 0.7);
        ui.tooltip.hide()
      },
      mousedown: (datum, index, nodes, event) => {
        this.dataFrame.selection.handleClick(i => {
          return this.fromColumn.get(i) === datum.source.id &&
            this.toColumn.get(i) === datum.target.id;
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
      this.computeChords();
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

    this.conf.innerRadius = size/2 - this.innerRadiusMargin;
    this.conf.outerRadius = size/2 - this.outerRadiusMargin;
    // this.chordConf.radius = d => (d.source.id === d.target.id) ? this.conf.outerRadius : null;

    circos.layout(this.data, this.conf);
    circos.chords('chords-track', this.chords, this.chordConf);
    circos.text('labels', this.labels, this.labelConf);
    circos.render();

    this.root.firstChild.style = 'position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);';
  }
}
