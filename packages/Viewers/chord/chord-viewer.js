import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import {select, scaleOrdinal} from 'd3';
import {layoutConf} from './configuration.js';


export class ChordViewer extends DG.JsViewer {

  constructor() {
    super();

    // Properties
    this.fromColumnName = this.string('fromColumnName');
    this.toColumnName = this.string('toColumnName');
    this.aggType = this.string('aggType', 'count');
    this.getProperty('aggType').choices = ['count', 'sum'];
    this.chordLengthColumnName = this.float('chordLengthColumnName');
    this.colorBy = this.string('colorBy', 'source');
    this.getProperty('colorBy').choices = ['source', 'target'];

    this.initialized = false;
    this.numColumns = [];
    this.strColumns = [];
    this.fromCol;
    this.toCol;
    this.data = [];
    this.conf = layoutConf;
    this.chords = [];
    this.chordConf = {};
  }

  init() {
    this.innerRadiusMargin = 60;
    this.outerRadiusMargin = 40;
    this.color = scaleOrdinal(DG.Color.categoricalPalette);
    this.chordConf.color = (datum, index) => DG.Color.toRgb(this.color(datum[this.colorBy]['id']));
    this.chordConf.opacity = 0.7;
    this.initialized = true;
  }

  testColumns() {
    return (this.strColumns.length >= 2 && this.numColumns.length >= 1);
  }

  onTableAttached() {
    this.init();

    this.strColumns = this.dataFrame.columns.toList().filter(col => col.type === 'string');
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

    this.categories = Array.from(new Set(this.fromCol.categories.concat(this.toCol.categories)));
    this.data = this.categories
      .sort((a, b) => this.freqMap[b] - this.freqMap[a])
      .map(s => {
        return {
          id: s,
          label: s,
          len: this.freqMap[s],
          color: DG.Color.toRgb(this.color(s))
        }
      });
  }

  computeChords() {
    this.chords.length = 0;
    let source = this.fromCol.getRawData();
    let fromCatList = this.fromCol.categories;
    let target = this.toCol.getRawData();
    let toCatList = this.toCol.categories;
    let aggVal = this.aggregatedTable.getCol('result').getRawData();
    let rowCount = this.aggregatedTable.rowCount;

    let aggTotal = {};
    for (let i = 0; i < rowCount; i++) {
      const from = this.fromCol.get(i);
      const to = this.toCol.get(i);
      if (from === to) {
        aggTotal[from] = (aggTotal[from] || 0) + aggVal[i];
      } else {
        aggTotal[from] = (aggTotal[from] || 0) + aggVal[i];
        aggTotal[to] = (aggTotal[to] || 0) + aggVal[i];
      }
    }

    let segments = {};
    this.data.forEach(s => {
      segments[s.id] = s;
      s.pos = 0;
    });

    for (let i = 0; i < rowCount; i++) {
      let sourceId = fromCatList[source[i]];
      let targetId = toCatList[target[i]];
      let sourceBlock = segments[sourceId];
      let targetBlock = segments[targetId];
      let sourceStep = sourceBlock.len * (aggVal[i] / aggTotal[sourceId]);
      let targetStep = targetBlock.len * (aggVal[i] / aggTotal[targetId]);

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
        value: aggVal[i]
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

    this.root.style = 'width: 100%; height: 100%;';
    this.root.id = 'chart';

    let circos = Circos({
      container: '#chart',
      width: size,
      height: size
    });

    this.conf.innerRadius = size/2 - this.innerRadiusMargin;
    this.conf.outerRadius = size/2 - this.outerRadiusMargin;
    this.chordConf.radius = d => (d.source.id === d.target.id) ? this.conf.outerRadius : null;

    circos.layout(this.data, this.conf);
    circos.chords('chords-track', this.chords, this.chordConf);
    circos.render();

    document.getElementById('chart').children[0]
      .setAttribute('style', 'position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);');
  }
}
