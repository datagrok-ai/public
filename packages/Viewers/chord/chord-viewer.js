import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import {layoutConf, chordConf} from './configuration.js';


export class ChordViewer extends DG.JsViewer {

  constructor() {
    super();

    // Properties
    this.fromColumnName = this.string('fromColumnName');
    this.toColumnName = this.string('toColumnName');
    this.aggType = this.string('aggType', 'count');
    this.getProperty('aggType').choices = ['count', 'sum'];

    this.initialized = false;
    this.numColumns = [];
    this.strColumns = [];
    this.fromCol;
    this.toCol;
    this.data = [];
    this.conf = layoutConf;
    this.chords = [];
    this.chordConf = chordConf;
  }

  init() {
    this.initialized = true;
  }

  testColumns() {
    return (this.strColumns.length >= 2 && this.numColumns.length >= 1);
  }

  onTableAttached() {
    this.init();

    let columns = this.dataFrame.columns.toList();
    this.strColumns = columns.filter(col => col.type === 'string');
    this.numColumns = columns.filter(col => ['double', 'int'].includes(col.type));

    // TODO: Choose the most relevant columns
    if (this.testColumns()) {
      this.fromColumnName = this.strColumns[0].name;
      this.toColumnName = this.strColumns[1].name;
      this.generateData();
    }

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));

    this.render();
  }

  onPropertyChanged(property) {
    super.onPropertyChanged(property);
    if (this.initialized && this.testColumns()) {
      this.generateData();
      this.render();
    }
  }

  detach() {
    this.subs.forEach(sub => sub.unsubscribe());
  }

  generateData() {

    this.fromColumn = this.dataFrame.getCol(this.fromColumnName);
    this.toColumn = this.dataFrame.getCol(this.toColumnName);
    this.conf.events = {
      mouseover: (datum, index, nodes, event) => {
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.fromColumn.get(i) === datum.id ||
            this.toColumn.get(i) === datum.id;
        }, event.x, event.y);
      },
      mouseout: () => ui.tooltip.hide(),
      mousedown: (datum, index, nodes, event) => {
        this.dataFrame.selection.handleClick(i => {
          return this.fromColumn.get(i) === datum.id ||
            this.toColumn.get(i) === datum.id;
        }, event);
      }
    };

    // For now, applies the aggregation function to the first numeric column
    this.aggregatedTable = this.dataFrame
      .groupBy([this.fromColumnName, this.toColumnName])
      .add(this.aggType, this.numColumns[0].name, 'result')
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
          color: "#80b1d3"
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

    let nMaxChords = {};
    for (let cat of this.categories) {
      for (let i = 0; i < rowCount; i++) {
        if (cat === this.fromCol.get(i) || cat === this.toCol.get(i)) {
          nMaxChords[cat] = (nMaxChords[cat] || 0) + 1;
        }
      }
    }

    this.segments = {};
    this.data.forEach(s => {
      this.segments[s.id] = s;
      s.step = s.len / nMaxChords[s.id];
      s.pos = 0;
    });

    for (let i = 0; i < rowCount; i++) {
      let sourceId = fromCatList[source[i]];
      let targetId = toCatList[target[i]];
      let sourceBlock = this.segments[sourceId];
      let targetBlock = this.segments[targetId];

      let chord = {
        source: {
          id: sourceId,
          start: sourceBlock.pos,
          end: sourceBlock.pos + sourceBlock.step
        },
        target: {
          id: targetId,
          start: targetBlock.pos,
          end: targetBlock.pos + targetBlock.step
        },
        value: aggVal[i]
      };

      if (sourceId === targetId) {
        chord.source.start = 0;
        chord.target.start = chord.source.end;
        chord.target.end = targetBlock.len;
      }

      this.chords.push(chord);

      if (sourceId === targetId) continue;
      sourceBlock.pos += sourceBlock.step;
      targetBlock.pos += targetBlock.step;
    }

    this.chordConf.events = {
      mouseover: (datum, index, nodes, event) => {
        ui.tooltip.showRowGroup(this.dataFrame, i => {
          return this.fromColumn.get(i) === datum.source.id &&
            this.toColumn.get(i) === datum.target.id;
        }, event.x, event.y);
      },
      mouseout: () => ui.tooltip.hide(),
      mousedown: (datum, index, nodes, event) => {
        this.dataFrame.selection.handleClick(i => {
          return this.fromColumn.get(i) === datum.source.id &&
            this.toColumn.get(i) === datum.target.id;
        }, event);
      }
    };

  }

  render() {

    if (!this.testColumns()) {
      this.root.innerText = "Not enough data to produce the result.";
      return;
    }

    $(this.root).empty();
    let width = this.root.parentElement.clientWidth;
    let height = this.root.parentElement.clientHeight;
    let size = Math.min(width, height);

    this.root.style = "width: 100%; height: 100%;";
    this.root.id = 'chart';

    let circos = Circos({
      container: '#chart',
      width: size,
      height: size
    });

    this.conf.innerRadius = size/2 - 60;
    this.conf.outerRadius = size/2 - 40;
    circos.layout(this.data, this.conf);

    this.computeChords();
    circos.chords('beta-track', this.chords, this.chordConf);
    circos.render();
  }
}
