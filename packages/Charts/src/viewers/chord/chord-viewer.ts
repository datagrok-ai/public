import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import * as Circos from 'circos';
import $ from 'cash-dom';

import {select, scaleLinear, scaleOrdinal, color, ScaleLinear, ScaleOrdinal} from 'd3';
import {layoutConf, topSort} from './utils';

export class ChordViewer extends DG.JsViewer {
  fromColumnName: string | null;
  toColumnName: string | null;
  aggType: string;
  chordLengthColumnName: string;
  colorBy: string;
  sortBy: string;
  direction: string;
  includeNulls: boolean;

  initialized: boolean;
  data: any;
  chords: any;
  segments: any;

  maxAvgNameLen?: number;
  maxLinkNumber?: number;
  minLengthThreshold?: number;
  innerRadiusMargin?: number;
  outerRadiusMargin?: number;

  gapScale?: ScaleLinear<number, number, never>;
  reservedColor?: string;
  colorScale?: ScaleOrdinal<string, number, never>;
  color?: (c: any) => string;

  chordOpacity?: number;
  highlightedChordOpacity?: number;

  conf: any;
  chordConf: any;
  labelConf?: {};

  strColumns?: DG.Column[];
  numColumns?: DG.Column[];

  aggregatedTable?: DG.DataFrame;

  distinctCols?: boolean;

  fromColumnAggr?: DG.Column;
  toColumnAggr?: DG.Column;
  chordWeights?: Int32Array | Float32Array | Float64Array | Uint32Array;
  rowCountAggr?: number;
  categories?: string[];

  indexes?: Int32Array;
  freqMap?: any;

  fromColumn?: DG.Column;
  toColumn?: DG.Column;

  constructor() {
    super();

    // Properties
    this.fromColumnName = this.string('fromColumnName');
    this.toColumnName = this.string('toColumnName');
    this.aggType = this.string('aggType', 'count', {choices: ['count', 'sum']});
    this.chordLengthColumnName = this.string('chordLengthColumnName');
    this.colorBy = this.string('colorBy', 'source', {choices: ['source', 'target']});
    this.sortBy = this.string('sortBy', 'topology', {choices: ['alphabet', 'frequency', 'topology']});
    this.direction = this.string('direction', 'clockwise', {choices: ['clockwise', 'counterclockwise']});
    this.includeNulls = this.bool('includeNulls', true);

    this.initialized = false;
    this.data = [];
    this.chords = [];
    this.segments = {};
  }

  init() {
    this.maxAvgNameLen = 15;
    this.maxLinkNumber = 1500;
    this.minLengthThreshold = 0.005;
    this.innerRadiusMargin = 80;
    this.outerRadiusMargin = 60;

    this.gapScale = scaleLinear([1, 100], [0.04, this.minLengthThreshold / 2]).clamp(true);
    this.reservedColor = 'rgb(127,127,127)';
    const colors = DG.Color.categoricalPalette.slice();
    colors.splice(colors.indexOf(4286545791), 1); // represents the reserved color
    this.colorScale = scaleOrdinal(colors);
    this.color = (c) => c ? DG.Color.toRgb(this.colorScale!(c)) : this.reservedColor!;

    this.chordOpacity = 0.7;
    this.highlightedChordOpacity = 0.9;

    this.conf = layoutConf;
    this.chordConf = {
      color: (datum: any) => this.color!(datum[this.colorBy]['label']),
      opacity: this.chordOpacity,
    };

    this.labelConf = {
      innerRadius: 1.02,
      style: {'font-size': 12, 'fill': this.reservedColor},
    };

    this.initialized = true;
  }

  _testColumns() {
    return (this.strColumns!.length >= 2 && this.numColumns!.length >= 1);
  }

  onTableAttached() {
    this.init();

    this.strColumns = this.dataFrame.columns.toList()
      .filter((col) => col.type === 'string')
      .sort((a, b) => a.categories.length - b.categories.length);
    this.numColumns = [...this.dataFrame.columns.numerical];

    // TODO: Choose the most relevant columns
    if (this._testColumns()) {
      this.fromColumnName = this.strColumns[0].name;
      this.toColumnName = this.strColumns[1].name;
      this.chordLengthColumnName = this.numColumns[0].name;
    }

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

    this.render();
  }

  onPropertyChanged(property: DG.Property) {
    if (this.initialized && this._testColumns()) {
      if (property.name === 'colorBy' && this.chords.length) this.render(false);
      else if (property.name === 'sortBy' && this.sortBy === 'alphabet' && !this.distinctCols) return;
      else this.render();
    }
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  _getFrequencies(sourceCol: DG.Column, targetCol: DG.Column, indexes: Int32Array) {
    const map: any = {};
    for (const i of indexes) {
      if (!this.includeNulls && (sourceCol.isNone(i) || targetCol.isNone(i))) continue;
      const from = sourceCol.isNone(i) ? '' : sourceCol.get(i);
      map[from] = (map[from] || 0) + 1;
      const to = targetCol.isNone(i) ? '' : targetCol.get(i);
      map[to] = (map[to] || 0) + 1;
    }
    return map;
  }

  _aggregate() {
    this.aggregatedTable = this.dataFrame
      .groupBy([this.fromColumnName!, (this.distinctCols) ? this.toColumnName! : ''])
      .whereRowMask(this.dataFrame.filter)
      .add(this.aggType as DG.AggregationType, this.chordLengthColumnName, 'result')
      .aggregate();

    if (!this.includeNulls) {
      this.aggregatedTable = this.aggregatedTable.rows
        .match(`${this.fromColumnName} regex .+ and ${this.toColumnName} regex .+`)
        .toDataFrame();
    }

    this.fromColumnAggr = this.aggregatedTable!.getCol(this.fromColumnName!);
    this.toColumnAggr = this.aggregatedTable!.getCol(this.toColumnName!);
    this.chordWeights = this.aggregatedTable!.getCol('result').getRawData();
    this.rowCountAggr = this.aggregatedTable!.rowCount;

    this.categories = Array.from(new Set(this.fromColumnAggr.categories.concat(this.toColumnAggr.categories)));
  }

  _updateBlockMap() {
    for (const prop of Object.getOwnPropertyNames(this.segments))
      delete this.segments[prop];

    this.data.forEach((s: any) => {
      this.segments[s.label] = {datum: s, targets: [], aggTotal: null, visited: false};
      s.pos = 0;
    });

    for (let i = 0; i < this.rowCountAggr!; i++) {
      const from = this.fromColumnAggr!.get(i);
      const to = this.toColumnAggr!.get(i);

      this.segments[from]['targets'].push(to);
      this.segments[from]['aggTotal'] = (this.segments[from]['aggTotal'] || 0) + this.chordWeights![i];
      if (from !== to) this.segments[to]['aggTotal'] = (this.segments[to]['aggTotal'] || 0) + this.chordWeights![i];
    }
  }

  _generateData() {
    this.data.length = 0;
    this.distinctCols = this.fromColumnName !== this.toColumnName;

    this.indexes = this.dataFrame.filter.getSelectedIndexes();
    this.freqMap = this._getFrequencies(this.fromColumn!, this.toColumn!, this.indexes);

    this.conf.events = {
      mouseover: (datum: any, index: number, nodes: any, event: MouseEvent) => {
        //@ts-ignore
        select(nodes[index]).select(`#${datum.id}`).attr('stroke', color(this.color(datum.label))!.darker());
        ui.tooltip.showRowGroup(this.dataFrame, (i) => {
          return this.dataFrame.filter.get(i) && (this.fromColumn!.get(i) === datum.label ||
            this.toColumn!.get(i) === datum.label);
        }, event.x, event.y);
      },
      mouseout: (datum: any, index: number, nodes: any, event: MouseEvent) => {
        select(nodes[index]).select(`#${datum.id}`).attr('stroke', 'none');
        ui.tooltip.hide();
      },
      mousedown: (datum: any, index: number, nodes: any, event: MouseEvent) => {
        this.dataFrame.selection.handleClick((i) => {
          return this.dataFrame.filter.get(i) && (this.fromColumn!.get(i) === datum.label ||
            this.toColumn!.get(i) === datum.label);
        }, event);
      },
    };

    if (this.aggType === 'sum' && this.chordLengthColumnName === null)
      this.chordLengthColumnName = this.numColumns![0].name;


    this._aggregate();

    this.conf.gap = this.gapScale!(this.categories!.length);

    if (this.sortBy !== 'topology') {
      this.categories!.sort((this.sortBy === 'frequency') ?
        (a, b) => this.freqMap![b] - this.freqMap[a] : undefined);
    }

    this.data = this.categories!.map((s, ind) => ({
      id: `id-${ind}`,
      label: s,
      len: this.freqMap[s],
      color: this.color!(s),
    }));

    if (this.distinctCols) this._updateBlockMap();

    if (this.sortBy === 'topology') {
      if (!this.distinctCols) {
        grok.shell.warning('Identical columns cannot be sorted topologically.');
        //@ts-ignore
        this.props.sortBy = 'alphabet';
      } else
        this.data = topSort(this.segments);
    }
    if (this.direction === 'counterclockwise') this.data.reverse();
  }

  _computeChords() {
    this.chords.length = 0;
    const source = this.fromColumnAggr!.getRawData();
    const fromCatList = this.fromColumnAggr!.categories;
    const target = this.toColumnAggr!.getRawData();
    const toCatList = this.toColumnAggr!.categories;

    for (let i = 0; i < this.rowCountAggr!; i++) {
      const sourceLabel = fromCatList[source[i]];
      const targetLabel = toCatList[target[i]];
      const sourceBlock = this.segments[sourceLabel]['datum'];
      const targetBlock = this.segments[targetLabel]['datum'];
      const sourceStep = sourceBlock.len * (this.chordWeights![i] / this.segments[sourceLabel]['aggTotal']);
      const targetStep = targetBlock.len * (this.chordWeights![i] / this.segments[targetLabel]['aggTotal']);

      this.chords.push({
        source: {
          id: sourceBlock.id,
          start: sourceBlock.pos,
          end: sourceBlock.pos + sourceStep,
          label: sourceLabel,
        },
        target: {
          id: targetBlock.id,
          start: targetBlock.pos,
          end: targetBlock.pos + targetStep,
          label: targetLabel,
        },
        value: this.chordWeights![i],
      });

      sourceBlock.pos += sourceStep;
      if (sourceLabel === targetLabel) continue;
      targetBlock.pos += targetStep;
    }

    this.chordConf.events = {
      mouseover: (datum: any, index: number, nodes: any, event: MouseEvent) => {
        select(nodes[index]).attr('opacity', this.highlightedChordOpacity!);
        ui.tooltip.showRowGroup(this.dataFrame, (i) => {
          return this.dataFrame.filter.get(i) &&
            this.fromColumn!.get(i) === datum.source.label &&
            this.toColumn!.get(i) === datum.target.label;
        }, event.x, event.y);
      },
      mouseout: (datum: any, index: number, nodes: any, event: MouseEvent) => {
        select(nodes[index]).attr('opacity', this.chordOpacity!);
        ui.tooltip.hide();
      },
      mousedown: (datum: any, index: number, nodes: any, event: MouseEvent) => {
        this.dataFrame.selection.handleClick((i) => {
          return this.dataFrame.filter.get(i) &&
            this.fromColumn!.get(i) === datum.source.label &&
            this.toColumn!.get(i) === datum.target.label;
        }, event);
      },
    };
  }

  _rotateLabels(labels: any) {
    labels.filter((d: any, i: number, nodes: any) => {
      return +(select(nodes[i]).attr('transform').match(/\d+\.?\d*/g)![0]) >= 180;
    }).selectAll('text')
      .attr('transform', (d: any, i: number, nodes: any) => select(nodes[i]).attr('transform') + ' rotate(180) ')
      .attr('text-anchor', 'end');
  }

  _cropLabels(labels: any) {
    labels.selectAll('text').each((d: any, i: number, nodes: any) => {
      const el = select(nodes[i]);
      let textLength = el.node().getComputedTextLength();
      let text = el.text();
      while (text.length && textLength > (this.outerRadiusMargin! - 15)) {
        text = text.slice(0, -1);
        el.text(text + '\u2026');
        textLength = el.node().getComputedTextLength();
      }
    });
  }

  _showErrorMessage(msg: string) {this.root.appendChild(ui.divText(msg, 'd4-viewer-error'));}

  render(computeData = true) {
    $(this.root).empty();

    if (!this._testColumns()) {
      this._showErrorMessage('Not enough data to produce the result.');
      return;
    }

    if (computeData) {
      // TODO: change when columnFilter setter is available
      this.fromColumn = this.dataFrame.getCol(this.fromColumnName!);
      this.toColumn = this.dataFrame.getCol(this.toColumnName!);
    }

    if (this.fromColumn!.type !== 'string' || this.toColumn!.type !== 'string') {
      this._showErrorMessage('Data of a non-string type cannot be plotted.');
      return;
    }
    if (this.fromColumn!.categories.length * this.toColumn!.categories.length > this.maxLinkNumber!) {
      this._showErrorMessage('Too many categories to render.');
      return;
    }

    if (computeData) {
      this._generateData();
      if (this.distinctCols) this._computeChords();
    }

    const width = this.root.parentElement!.clientWidth;
    const height = this.root.parentElement!.clientHeight;
    const size = Math.min(width, height);

    const circos = new Circos({
      container: this.root,
      width: size,
      height: size,
    });

    this.conf.innerRadius = Math.max(0, size/2 - this.innerRadiusMargin!);
    this.conf.outerRadius = Math.max(0, size/2 - this.outerRadiusMargin!);
    // this.chordConf.radius = d => (d.source.id === d.target.id) ? this.conf.outerRadius : null;

    circos.layout(this.data, this.conf);

    // let smallBlocks = this.data.some(d => d.end - d.start < this.minLengthThreshold);
    // if (smallBlocks) {
    //   this.root.appendChild(ui.divText('Too many categories to render.', 'd4-viewer-error'));
    //   return;
    // }

    if (this.distinctCols)
      circos.chords('chords-track', this.chords, this.chordConf);


    const avgNameLen = this.categories!.reduce((a, b) => a + b.length, 0) / this.categories!.length;
    const showLabels = avgNameLen < this.maxAvgNameLen!;

    if (showLabels) {
      circos.text('labels', this.data.map((d: any) => ({
        block_id: d.id,
        position: this.freqMap![d.label] / 2,
        value: d.label,
      })), this.labelConf);
    }

    circos.render();

    if (showLabels) {
      const labels = select(this.root).selectAll('.block');
      this._rotateLabels(labels); // fix label rotation past 180
      this._cropLabels(labels);
    }

    (this.root.firstChild! as HTMLElement).setAttribute('style',
      'position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);');
  }
}
