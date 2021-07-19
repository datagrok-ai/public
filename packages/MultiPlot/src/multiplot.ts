import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
import {initTimeLine} from './timeLinesRender';
// import * as deb from "./../debug.js";

export class MultiPlotViewer extends DG.JsViewer {
  // properties
  private defaultTitleHeight: string = this.string('defaultTitleHeight', '25px');
  private backColor: number = this.int('backColor', 0xffffff);
  private verticalLinesColor: number = this.int('verticalLinesColor', 0xf0f0f0);
  private verticalLines: boolean = this.bool('verticalLines', true);
  private showControls: boolean = this.bool('showControls', true);
  private timeLineWidth: number = this.int('timeLineWidth', 1);
  private timeLineRadius: number = this.int('timeLineRadius', 3);

  private paramA = this.string('paramA', 'string inside');
  private typeComboElements = [];
  private showHideElements = [];
  private closeElements = [];
  private plots = [];
  private box: any;
  private echart: any;
  private echarts: any = echarts;
  private echartOptions: any;
  private visiblePlotsCount: number;
  private initTimeLine = initTimeLine;
  private timeLineSeries: any;
  private timeLineIndex: number = 11;
  private statusChartIndex: number = 10;
  private timeLineOverlapShift = 4;
  // private circleRange: number = 3;
  // private timeLineWidth: number = 1;
  private timeLinesData: any = [];
  private statusChartData: any = [];
  private categoryColors: any = [];
  private currentRowCount: number = 0;
  private selection: any = [];
  private plotTitleHighMargin: number = 10;
  private plotTitleLowMargin: number = 5;
  private controlsTopShift: number = -4;
  private headerHiddenShift: number = -10;
  private tables = {};
  private isTablesLoaded: number = 0;
  private tooltipOffset: number = 10;
  private visibleIndexes: number[];
  private marker = this.string('marker', 'ring', {choices: ['circle', 'rect', 'ring', 'diamond']});
  private markerSize = this.int('markerSize', 6);
  private lineWidth = this.int('lineWidth', 2);
  private markerPosition = this.string('markerPosition', 'main line',
      {choices: ['main line', 'above main line', 'scatter']});
  private selectionColor = DG.Color.toRgb(DG.Color.selectedRows);
  private options = {
    series: [
      {
        table: 'ae__2__lb__2_',
        title: 'title1',
        type: 'scatter',
        x: 'AESTDY',
        //    y: 'LBTEST',
        y: 'LBSTRESN',

        yType: 'value',
        color: 'red',
        markerShape: 'square',
        height: '1flex',
        show: 1,
      },
      {
        table: 'ae__2__lb__2_',
        title: 'title11',
        type: 'timeLine',
        x: 'LBTEST',
        y: ['AESTDY', 'LBDY'],
        yType: 'category',
        color: 'red',
        markerShape: 'square',
        height: '2flex',
        show: 1,
      },
      {
        table: 'ae__2__lb__2_',
        title: 'title2',
        type: 'scatter',
        x: 'AESTDY',
        y: 'LBSTRESN',
        yType: 'value',
        color: 'red',
        markerShape: 'square',
        height: '20%',
        show: 1,
      },
      /*   {
        table: 'ae__2__lb__2_',
        title: 'title1',
        type: 'timeLine',
        x: 'AETERM',
        y: ['AESTDY', 'LBDY'],
        yType: 'category',
        color: 'red',
        markerShape: 'square',
        height: '2flex',
        show: 1,
      },*/
    ],
  }

  constructor() {
    super();
    console.log('------------------------------------- MULTIPLOT ------------------------------');
    console.log('this.root: ', this.root);

    this.categoryColors = DG.Color.categoricalPalette.map(DG.Color.toRgb);
    this.plots = this.options.series;
    const allTypes = ['scatter', 'line', 'bar', 'timeLine'];
    // init plots with some data y=x^2, later can be replaced
    this.plots.map((e) => {
      e.series = {
        type: 'scatter',
        data: [],
      };
    });

    this.init();

    this.box = this.root.getBoundingClientRect();
    ui.onSizeChanged(this.root).subscribe((e) => {
      if (!this.echartOptions) return 0;
      this.updateOptionsPositions();
      if (this.echart) {
        this.updatePlots();
        this.render();
        this.echart.resize();
      }
    });
  }

  onPropertyChanged(property: DG.Property): void {
    const name = property.name;
    const val = property.get(this);
    if (name === 'defaultTitleHeight') {
      this.defaultTitleHeight = val;
      this.updateHeight();
      this.render();
    }
    if (name === 'backColor') {
    }
    if (name === 'verticalLines') {
    }
    if (name === 'showControls') {
      this.setControlsVisibility();
    }

    this.updatePlots();
    this.render();
    // super.onPropertyChanged(property);
  }

  normalize100(column: DG.Column): DG.Row[] {
    const r = [];
    const delta = column.max - column.min;
    const rawData = column.getRawData();
    for (let i = 0; i < column.length; i++) {
      // avoid case when max == min
      const t = delta ? (rawData[i] - column.min) * 100 / delta : rawData[i];
      r.push(t);
    }
    return r;
  }

  // Assign size and position of viewer's parts (charts, controls)
  // input: height of viewer's container
  // also takes data from this.plots:
  //    [{height: '10px', title: 'title0'},{height: '20%'}, {height: 'free'}]
  // output: [{top: '20px', height: '10px', topTitle: '0px', titleHeight}, ...]
  parseHeight(height: number): any[] {
    const r = [];
    let floatHeight = height;
    let totalFlexPoints = 0;

    for (let i = 0; i < this.plots.length; i++) { // first cycle
      const titleHeight = this.plots[i].title ? parseFloat(this.defaultTitleHeight) : 0;
      let plotHeight = 0;
      let flexPoints = 0;
      if (this.plots[i].height.includes('px')) {
        plotHeight = parseFloat(this.plots[i].height);
      }
      if (this.plots[i].height.includes('%')) {
        plotHeight = parseFloat(this.plots[i].height) * height / 100;
      }
      if (this.plots[i].height.includes('flex') && this.plots[i].show) {
        flexPoints = parseFloat(this.plots[i].height);
        totalFlexPoints += flexPoints;
      };
      if (this.plots[i].show == 0) plotHeight = 0;

      r.push({
        flexPoints: flexPoints,
        height: plotHeight,
        titleHeight: titleHeight,
        titleText: this.plots[i].title,
      });
      floatHeight -= plotHeight + titleHeight +
        1 * (this.plotTitleHighMargin + this.plotTitleLowMargin);
    } // first cycle

    // weight (in pixels) of one flex point
    const flexPointHeight = floatHeight / totalFlexPoints;

    let currentTop = 0;
    let echartIndex = 0; // index in echarts library
    // second cycle (set positions and height)
    for (let j = 0; j < r.length; j++) {
      currentTop += this.plotTitleHighMargin;
      r[j].titleTop = currentTop;
      currentTop += r[j].titleHeight;
      if (this.plots[j].show) {
        currentTop += this.plotTitleLowMargin;
        r[j].top = currentTop;
        if (r[j].flexPoints > 0.00001) {
          r[j].height = r[j].flexPoints * flexPointHeight;
        };
        currentTop += r[j].height;
        this.plots[j].i = echartIndex;
        echartIndex++;
      };
    }

    // correct height of the last visible plot to avoid overlap of controls
    let j = this.plots.length - 1;
    while (j > -1 && !this.plots[j].show) {
      j--;
    }
    if (j < this.plots.length - 1 && j >= 0) {
      r[j].height += this.headerHiddenShift;
    }
    return r;
  } // parseHeight

  // output: filled echart options array only for height and position
  updateOptionsPositions(): void {
    this.box = this.root ? this.root.getBoundingClientRect() : {height: 300};
    const heightData = this.parseHeight(this.box.height - 30);
    let visibleIndex = 0;
    for (let i = 0; i < this.plots.length; i++) {
      this.echartOptions.title[i].top = (heightData[i].titleTop) + 'px';
      this.echartOptions.title[i].left = '10px';
      this.echartOptions.title[i].text = heightData[i].titleText;
      this.echartOptions.title[i].textStyle = {
        'height': heightData[i].titleHeight,
        'fontSize': (parseFloat(this.defaultTitleHeight) * .6),
      };
      if (this.typeComboElements[i]) {
        this.typeComboElements[i].style.top = heightData[i].titleTop + this.controlsTopShift + 'px';
      }

      if (this.showHideElements[i]) {
        this.showHideElements[i].style.top = heightData[i].titleTop + 7 + this.controlsTopShift + 'px';
      }

      if (this.closeElements[i]) {
        this.closeElements[i].style.top = heightData[i].titleTop + 7 + this.controlsTopShift + 'px';
      }

      if (!this.plots[i].show) continue;

      // only for visible plots
      this.echartOptions.grid[visibleIndex].y2 = 22;
      this.echartOptions.grid[visibleIndex].top = heightData[i].top + 'px';
      this.echartOptions.grid[visibleIndex].bottom = '';
      this.echartOptions.grid[visibleIndex].show = heightData[i].show ? true : false;
      this.echartOptions.grid[visibleIndex].height = heightData[i].height + 'px';
      visibleIndex++;
    } // for i
  } // updateOptionsPositions

  // takes data from this.plots and fills options for echart library
  updatePlots(): void {
    function toColor(num: number): string {
      num >>>= 0;
      const b = num & 0xFF;
      const g = (num & 0xFF00) >>> 8;
      const r = (num & 0xFF0000) >>> 16;
      // a = ( (num & 0xFF000000) >>> 24 ) / 255 ;
      const a = 1;
      return 'rgba(' + [r, g, b, a].join(',') + ')';
    }
    this.clearPlots();
    if (this.isTablesLoaded === 0) return;
    this.echartOptions.backgroundColor = toColor(this.backColor);

    // update positions
    this.updateOptionsPositions();
    this.echartOptions.series = [];
    let visibleIndex = 0;
    this.visibleIndexes = [];
    for (let i = 0; i < this.plots.length; i++) {
      this.echartOptions.title[i].left = '10px';
      this.echartOptions.title[i].text = this.plots[i].title;
      if (!this.plots[i].show) continue;

      this.echartOptions.grid[visibleIndex].left = '19%';
      this.echartOptions.grid[visibleIndex].right = '4%';
      this.echartOptions.grid[visibleIndex].show = this.plots[i].show;
      this.echartOptions.xAxis[visibleIndex].gridIndex = visibleIndex;
      this.echartOptions.xAxis[visibleIndex].type = 'value';
      //    this.echartOptions.xAxis[i].axisTick = { inside: true, length: 1000 };
      this.echartOptions.xAxis[visibleIndex].show = false;
      this.echartOptions.xAxis[visibleIndex].min = 0;
      this.echartOptions.xAxis[visibleIndex].max = 200;
      this.echartOptions.yAxis[visibleIndex].gridIndex = visibleIndex;
      this.echartOptions.yAxis[visibleIndex].type = this.plots[i].yType || 'value';
      this.echartOptions.yAxis[visibleIndex].show = this.plots[i].show;
      this.echartOptions.yAxis[visibleIndex].triggerEvent = true;

      let currentSeries = {
        type: this.plots[i].series.type,
        //    show: this.plots[i].show,
        large: 'true',
        gridIndex: visibleIndex,
        symbol: this.plots[i].series.symbol || (() => { }),
        itemStyle: this.plots[i].series.itemStyle ?
          {color: this.plots[i].series.itemStyle} : {},
        yAxis: {type: this.plots[i].yType || 'value'},
        xAxisIndex: visibleIndex,
        yAxisIndex: visibleIndex,
        data: this.plots[i].series.data,
        //   coordinateSystem: 'cartesian2d',
        encode: {x: 0, y: 1},
      };

      if (this.plots[i].type === 'timeLine') {
        this.plots[i].timeLinesSeries = this.initTimeLine(visibleIndex, this);
        this.plots[i].subjectCol = this.tables[this.plots[i].table].getCol(this.plots[i].x);
        this.plots[i].subjects = this.plots[i].subjectCol.categories;
        this.plots[i].subjBuf = this.plots[i].subjectCol.getRawData();
        currentSeries = this.plots[i].timeLinesSeries;
        currentSeries.xAxisIndex = visibleIndex;
        currentSeries.yAxisIndex = visibleIndex;
        currentSeries.gridIndex = visibleIndex;
        this.echartOptions.xAxis[visibleIndex].type = 'value';
        this.echartOptions.yAxis[visibleIndex].type = 'category';
      }

      this.echartOptions.series.push(currentSeries);
      this.visibleIndexes.push(i);
      visibleIndex++;
    } // for i<this.plots.length

    this.echartOptions.responsive = false;
    if (visibleIndex === 0) {
      return;
    }

    this.echartOptions.xAxis[visibleIndex - 1].axisTick = {
      inside: true,
      length: this.verticalLines ? 2000 : 5,
      lineStyle: {color: toColor(this.verticalLinesColor)},
    };
    this.echartOptions.xAxis[visibleIndex - 1].show = true;
    this.echartOptions.xAxis[visibleIndex - 1].type = 'value';
  } // updatePlots

  // only updates heights
  updateHeight() : void {
    this.updateOptionsPositions();
    this.render();
  }

  // fill echart options with arrays of empty objects
  clearPlots() : void {
    this.visiblePlotsCount = 0; // number of visible charts
    for (let i = 0; i < this.plots.length; i++) {
      if (this.plots[i].show) {
        this.visiblePlotsCount++;
      }
    }
    function createEmptyObjects(n: number) : any {
      return new Array(n).fill(null).map(() => {
        return {};
      });
    }

    this.echartOptions = {
      title: createEmptyObjects(this.plots.length),
      /*
      'tooltip2': {
        trigger: 'axis',
        axisPointer: {type: 'shadow'},
      },
      */
      tooltip: {
        trigger: 'axis',
        showContent: false,
        axisPointer: {type: 'shadow'},
      },
      textStyle: {
        fontFamily: 'Roboto'
      },
      grid: createEmptyObjects(this.visiblePlotsCount),
      xAxis: createEmptyObjects(this.visiblePlotsCount),
      yAxis: createEmptyObjects(this.visiblePlotsCount),
      dataZoom: [{
        type: 'inside',
        xAxisIndex: Array.from(Array(this.visiblePlotsCount).keys()),
        start: 0,
        end: 100,
      }],
      animation: false,
      series: createEmptyObjects(this.visiblePlotsCount),
    };
    for (let i = 0; i < this.echartOptions.series.length; i++) {
      this.echartOptions.series[i].type = 'scatter';
      this.echartOptions.series[i].data = [];
    }
  } // clearPlots

  onEvent(e: DG.Events): void {
  }

  addMenu(): void {
    grok.events.onContextMenu.subscribe((args) => {
      //    if (!(args.args.context instanceof DG.Viewer)) { return 0; };

      // get opened tables (names, tabs);
      const tabs = {};
      const names = grok.shell.tables.map((t) => t.name);
      for (let i = 0; i < names.length; i++) {
        tabs[names[i]] = grok.shell.tables[i];
      }

      const callback = (item) => {
        const table = tabs[item];
        const nCols = Array.from(table.columns.numerical);
        const colNames = nCols.map((e: { name: string }) => e.name);
        const r0raw = table.getCol(colNames[0]).getRawData();
        const column0 = this.normalize100(table.getCol(colNames[0]));
        const column1 = table.getCol(colNames[1]).getRawData();
        const data = [];
        for (let i = 0; i < column0.length; i++) {
          data.push([column0[i], column1[i]]);
        };
        this.plots.push({
          height: '1flex', title: item + ' ' + colNames[1] + '( ' + colNames[0] + ' )',
          series: {
            type: 'scatter',
            large: true,
            data: data,
          },
          show: 1,
        });
        this.createElements();
        this.updatePlots();
        this.render();
      }; // callback

      // add context menu items (right click);
      const menu = args.args.menu.group('New plot');
      menu.items(names, callback);
    }); // onContextMenu subscribe
  }

  getStatusChartData(table: DG.DataFrame, indexes: Int32Array): any[] {
    const rez = [];
    for (let i = 0; i < table.rowCount; i++) {
      const row = table.row(i);

      const valTestName = row['LBTEST'];
      const valStartTime = row['AESTDY'];

      rez.push([valStartTime, valTestName]);
    }
    this.statusChartData = rez;
    return rez;
  }

  // build 2d array for series.data of echart
  getUniversalData(table: DG.DataFrame, fieldsNames: string[], indexes: Int32Array): any[] {
    const r = [];

    function getRowFields(row: DG.Row): any[] {
      const fields = [];
      for (let i = 0; i < fieldsNames.length; i++) {
        fields.push(row[fieldsNames[i]]);
      }
      return fields;
    };

    if (indexes) {
      for (let ind = 0; ind < indexes.length; ind++) {
        const row = table.row(indexes[ind]);
        const fields = getRowFields(row);
        r.push(fields);
      }
    } else {
      for (let i = 0; i < table.rowCount; i++) {
        const row = table.row(i);
        const fields = getRowFields(row);
        r.push(fields);
      }
    }
    return r;
  }

  onTableAttached(): void {
    this.addMenu();
    const tableArray = grok.shell.tables;
    this.tables = {};
    for (let i=0; i<tableArray.length; i++) {
      this.tables[tableArray[i].name] = tableArray[i];
    }

    if (this.isTablesLoaded) {
      return;
    }
    this.isTablesLoaded = 1;
    this.updateFilter();
    this.render();

    // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    // this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    // @ts-ignore
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => {
      this.selection = this.dataFrame.selection.getSelectedIndexes();
      this.plots.map((plot) => {
        const defaultColor = this.categoryColors[0];
        const selectionColor = this.categoryColors[1];
        plot.series.itemStyle = (e, i) => {
          let color = defaultColor;
          if (this.selection.includes(e.dataIndex)) {
            color = selectionColor;
          }
          return color;
        };
      });
      this.updatePlots();
      this.render();
    }));   
    // @ts-ignore
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => {
      this.updateFilter();
      this.render();
    }));

  } // table attached

  updateFilter(): void {
    for (let i=0; i< this.plots.length; i++) {
      const tableName = this.plots[i].table;
      const table = this.tables[tableName];
      const indexes = table.filter.getSelectedIndexes();
      const x = this.plots[i].x;
      const y = this.plots[i].y;
      const xArray = Array.isArray(x) ? x : [x];
      const yArray = Array.isArray(y) ? y : [y];

      const data = this.getUniversalData(
          table,
          // [this.plots[i].x, this.plots[i].y],
          xArray.concat(yArray),
          indexes,
      );
      this.plots[i].series.data = data;
      if (this.plots[i].type != 'timeLine') continue;
      this.timeLinesData = data;
      this.timeLineIndex = i;
    }
    this.updatePlots();
  }

  isGroup(componentIndex: number) : boolean {
    const type = this.plots[this.visibleIndexes[componentIndex]].type;
    if (type === 'scatter' || type === 'line') return false;
    return true;
  }

  init(): void {
    this.plots.map((plot) => {
      const defaultColor = this.categoryColors[0];
      const selectionColor = this.categoryColors[1];
      plot.series.itemStyle = (e, i) => {
        let color = defaultColor;
        if (this.selection.includes(e.dataIndex)) {
          color = selectionColor;
        }
        return color;
      };
    });
    if (!this.echart) this.echart = echarts.init(this.root, null, {renderer: 'canvas'});

    this.echart.on('dataZoom', (e) => {
      console.log('zoom ', e);
    });

    this.echart.on('click', (params) => {
      console.log('click params ', params);
      const iPlot : number = this.visibleIndexes[params.componentIndex];
      const table : DG.DataFrame = this.tables[this.plots[iPlot].table];
      const subjBuf = this.plots[iPlot];
      const x = params.event.event.x + this.tooltipOffset;
      const y = params.event.event.y + this.tooltipOffset;
      const xColName = this.plots[iPlot].x;
      const yColName = this.plots[iPlot].y;
      const colNames = [xColName, yColName];
      table.selection.handleClick( (i) => {
        console.log('params value');
        if (params.componentType === 'yAxis') {
          return this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]] === params.value;
        }
        if (params.componentType === 'series') {
          if (this.isGroup(params.componentIndex)) {
            return params.value[0] ===
              this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]];
          } else {
            return params.dataIndex === i;
          }
        }
      }, params.event.event);
    });

    this.echart.on('mouseover', (params) => {
      const iPlot : number = this.visibleIndexes[params.componentIndex];
      const table : DG.DataFrame = this.tables[this.plots[iPlot].table];
      const subjBuf = this.plots[iPlot];
      const x = params.event.event.x + this.tooltipOffset;
      const y = params.event.event.y + this.tooltipOffset;
      const xColName = this.plots[iPlot].x;
      const yColName = this.plots[iPlot].y;
      const colNames = [xColName, yColName];
      const val = params.value[1];

      if (params.componentType === 'yAxis') {
        if (this.isGroup(params.componentIndex)) {
          ui.tooltip.showRowGroup(table, (i) => {
            return params.value === this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]];
          }, x, y);
        }
      }

      if (params.componentType === 'series') {
        if (!this.isGroup(params.componentIndex)) {
          ui.tooltip.show(ui.divV(
              params.data.map((e, i) => ui.div([colNames[i] + ': ' + e + ''])),
          ), x, y);
        } else {
          ui.tooltip.showRowGroup(table, (i) => {
            return params.value[0] === this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]]; // &&
            //            params.value[1] === (this.startCol.isNone(i) ? null : this.startBuf[i]) &&
            //          params.value[2] === (this.endCol.isNone(i) ? null : this.endBuf[i]);
          }, x, y);
        }
      } // series
    }); // mouseover

    this.echart.on('mouseout', () => ui.tooltip.hide());

    this.clearPlots();
    this.createElements();
    this.updateOptionsPositions();
    //   this.updatePlots();
    //  this.render();
    this.render();
  } // init

  deleteElements(): void {
    this.typeComboElements.map((e) => e.remove());
    this.showHideElements.map((e) => e.remove());
    this.closeElements.map((e) => e.remove());
  }

  setControlsVisibility() : void {
    this.typeComboElements.map((e) => {
      e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
    this.showHideElements.map((e) => {
      e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
    this.closeElements.map((e) => {
      e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
  }

  createElements(): void {
    this.deleteElements();

    // create comboboxes to choose plot types
    this.typeComboElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      const inputPlotType: any = ui.choiceInput('', 'scatter', ['scatter', 'line', 'bar'], ((i) => (event) => {
        this.plots[i].series.type = event;
        this.updatePlots();
        this.render();
      })(i));
      this.typeComboElements.push(inputPlotType.root);
      this.root.appendChild(inputPlotType.root);
      inputPlotType.root.style.position = 'absolute';
      inputPlotType.root.style.right = '28px';
      inputPlotType.root.style['flex-direction'] = 'row';
      inputPlotType.root.style.top = (40 * i) + 'px';
    }

    // create checkboxes for show/hide plots
    this.showHideElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      const inputPlotType: any = ui.div([ui.iconFA('angle-right'), ui.iconFA('angle-down')]);
      const showHideIcons = inputPlotType.querySelectorAll('i');
      showHideIcons[0].style.display = 'none';
      inputPlotType.showSwitch = 1;
      inputPlotType.addEventListener('click', ((i) => (e) => {
        const div = e.target.parentNode;
        div.showSwitch = 1 - div.showSwitch;
        const displays = ['', 'none'];
        const els = div.querySelectorAll('i');
        const isShown = div.showSwitch;
        els[0].style.display = displays[isShown];
        els[1].style.display = displays[1 - isShown];
        this.plots[i].show = isShown;
        this.updateFilter();
        this.updatePlots();
        this.render();
      })(i));
      this.typeComboElements.push(inputPlotType);
      this.root.appendChild(inputPlotType);
      inputPlotType.style.position = 'absolute';
      inputPlotType.style.left = '3px';
      inputPlotType.style['flex-direction'] = 'row';
      inputPlotType.style.top = (40 * i) + 'px';
      this.showHideElements.push(inputPlotType);
    }

    // create close 'X' icons
    this.closeElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      const inputClose: any = ui.icons.close(((i) => () => {
        //    grok.shell.info('click' + i);
        this.plots.splice(i, 1);
        this.createElements();
        this.updateHeight();
        this.updatePlots();
        this.render();
      })(i), 'Close');
      inputClose.style.position = 'absolute';
      inputClose.style.right = '15px';
      inputClose.style.top = (40 * i) + 'px';
      inputClose.style.flexDirection = 'row';
      this.closeElements.push(inputClose);
      this.root.appendChild(inputClose);
    }
  } // createElements

  render(): void {
    console.log('set echart options: ', this.echartOptions);
    if (this.echart) {
      this.echart.clear();
      // this.echart.setOption(this.echartOptions, true);
      setTimeout(() => this.echart.setOption(this.echartOptions, true), 200);
    }
  }
}
