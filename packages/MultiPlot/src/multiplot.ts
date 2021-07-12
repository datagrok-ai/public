import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
import {DataFrame, Property, Viewer} from 'datagrok-api/dg';
// import * as deb from "./../debug.js";

export class MultiPlotViewer extends DG.JsViewer {
  // properties
  private defaultTitleHeight: string = this.string('defaultTitleHeight', '25px');
  private BackColor: string = this.string('BackColor', 'white');
  private verticalLinesColor: string = this.string('verticalLinesColor', '#cccccc');
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
  private echartOptions: any;
  private visiblePlotsCount: number;
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
  private tables = {}
  private isTablesLoaded: number = 0;
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
        title: 'title1',
        type: 'timeLine',
        x: 'LBTEST',
        y: ['AESTDY', 'LBDY'],
        yType: 'category',
        color: 'red',
        markerShape: 'square',
        height: '1flex',
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
    ],
  }

  constructor() {
    super();
    console.log('------------------------------------- MULTIPLOT ------------------------------');
    console.log('this.root: ', this.root);
    console.log('this a ', this.paramA);

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
        //   console.log(JSON.stringify(this.plots, null, 2));
        this.updatePlots();
        this.setEchartOptions();
        this.echart.resize();
        console.log('this this.paramA -------------------------', this.paramA);
      }
    });
  }

  setEchartOptions(): void {
    console.log('set echart options: ', this.echartOptions);
    if (this.echart) {
      // this.echart.clear();
      // this.echart.setOption(this.echartOptions, true);
      setTimeout(() => this.echart.setOption(this.echartOptions, true), 200);
    }
  }

  onPropertyChanged(property: DG.Property): void {
    const name = property.name;
    const val = property.get(this);
    console.log('property changed: ', name, val, this.verticalLinesColor);
    if (name === 'defaultTitleHeight') {
      this.defaultTitleHeight = val;
      this.updateHeight();
      this.setEchartOptions();
    }
    if (name === 'BackColor') {
      console.log('BackColor', this.BackColor);
    }
    if (name === 'verticalLines') {
      console.log(this.verticalLines);
    }
    if (name === 'showControls') {
      this.setControlsVisibility();
    }

    this.updatePlots();
    this.setEchartOptions();
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

    console.log('currentTop  :', currentTop, r);
    return r;
  } // parseHeight

  // output: filled echart options array only for height and position
  updateOptionsPositions(): void {
    this.box = this.root ? this.root.getBoundingClientRect() : {height: 300};
    const heightData = this.parseHeight(this.box.height - 30);
    console.log('heightData:  ', heightData);

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
    this.clearPlots();

    this.echartOptions.backgroundColor = this.BackColor;

    // update positions
    console.log('update plots positions', this.plots);
    this.updateOptionsPositions();
    this.echartOptions.series = [];
    let visibleIndex = 0;
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
        currentSeries = this.plots[i].timeLinesSeries;
        currentSeries.xAxisIndex = visibleIndex;
        currentSeries.yAxisIndex = visibleIndex;
        currentSeries.gridIndex = visibleIndex;
        this.echartOptions.xAxis[visibleIndex].type = 'value';
        this.echartOptions.yAxis[visibleIndex].type = 'category';
      }

      this.echartOptions.series.push(currentSeries);
      visibleIndex++;
    } // for i<this.plots.length

    this.echartOptions.xAxis[visibleIndex - 1].axisTick = {
      inside: true,
      length: this.verticalLines ? 2000 : 5,
      lineStyle: {color: this.verticalLinesColor},
    };
    this.echartOptions.xAxis[visibleIndex - 1].show = true;
    this.echartOptions.xAxis[visibleIndex - 1].type = 'value';
    this.echartOptions.responsive = false;

    console.log('echart options: ', this.echartOptions);
  } // updatePlots

  // only updates heights
  updateHeight() : void {
    this.updateOptionsPositions();
    this.setEchartOptions();
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
      'title': createEmptyObjects(this.plots.length),
      'tooltip': {
        trigger: 'axis',
        axisPointer: {type: 'shadow'},
      },
      'grid': createEmptyObjects(this.visiblePlotsCount),
      'xAxis': createEmptyObjects(this.visiblePlotsCount),
      'yAxis': createEmptyObjects(this.visiblePlotsCount),
      'dataZoom': [{
        type: 'inside',
        xAxisIndex: Array.from(Array(this.visiblePlotsCount).keys()),
        start: 0,
        end: 100,
      }],
      'animation': false,
      'series': createEmptyObjects(this.visiblePlotsCount),
    };
    for (let i = 0; i < this.echartOptions.series.length; i++) {
      this.echartOptions.series[i].type = 'scatter';
      this.echartOptions.series[i].data = [];
    }

    console.log('clear plots end', this.echartOptions);
  } // clearPlots

  onEvent(e: DG.Events): void {
    console.log('event: ', e);
  }

  addMenu(): void {
    grok.events.onContextMenu.subscribe((args) => {
      console.log('args: ', args.args);
      //    if (!(args.args.context instanceof DG.Viewer)) { return 0; };

      // get opened tables (names, tabs);
      const tabs = {};
      const names = grok.shell.tables.map((t) => t.name);
      for (let i = 0; i < names.length; i++) {
        tabs[names[i]] = grok.shell.tables[i];
      }

      const callback = (item) => {
        console.log('item ', item);
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
        this.setEchartOptions();
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
    console.log('table attached');
    console.log('this ', this);
    console.log('THIS:PLOTS: ', this.plots);
    console.log('tables: ', grok.shell.tables);
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
      this.setEchartOptions();
      console.log(this.selection);
      this.render();
    }));

    // @ts-ignore
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => {
      console.log('filter change');
      this.updateFilter();
      this.setEchartOptions();

      this.render();
    }));
  } // table attached

  updateFilter(): void {
    for (let i=0; i< this.plots.length; i++) {
      const tableName = this.plots[i].table;
      const table = this.tables[tableName];
      const indexes = table.filter.getSelectedIndexes();
      console.log(indexes);
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
      this.plots[i].timeLinesSeries = this.initTimeLine();
    }
    this.updatePlots();
  }

  applyFilter(filter: DG.BitSet): void {
    this.clearPlots();
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
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
    this.clearPlots();
    this.createElements();
    this.updateOptionsPositions();
    //   this.updatePlots();
    //  this.setEchartOptions();
    this.render();
  }

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
        console.log('changed ', event, i);
        this.plots[i].series.type = event;
        this.updatePlots();
        this.setEchartOptions();
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
        //     console.log('click ', i, div.showSwitch);
        div.showSwitch = 1 - div.showSwitch;
        const displays = ['', 'none'];
        const els = div.querySelectorAll('i');
        const isShown = div.showSwitch;
        els[0].style.display = displays[isShown];
        els[1].style.display = displays[1 - isShown];
        this.plots[i].show = isShown;
        this.updatePlots();
        this.setEchartOptions();
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
        this.setEchartOptions();
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

  }

  initTimeLine(): any {
    console.log('init time line :', this.timeLinesData);
    let count = 0;
    const renderFailed = 1234;
    // this.timeLineSeries = {
    const timeLinesSeries = {
      type: 'custom',
      progressive: 0,
      animation: false,
      renderItem: ((t: any, echarts: any) => (params, echartAPI) => {
        const data = this.timeLinesData;
        if (!data || data.length == 0) return renderFailed;
        const customDebug = false;
        if (customDebug) console.log('renderItem data: ', data);
        if (customDebug) console.log('custom render ', params, echartAPI);
        const av0 = echartAPI.value(0);
        const av1 = echartAPI.value(1);
        const av2 = echartAPI.value(2);
        if (customDebug) console.log('values: ', av0, av1, av2);
        const gridTopRaw = this.echartOptions.grid[this.timeLineIndex].top;
        const gridTop = parseFloat(gridTopRaw);
        if (customDebug) console.log('gridtop ', gridTop);
        let overlap = false;
        if (params.dataIndex > 0 && data[params.dataIndex - 1][0] === data[params.dataIndex][0] &&
          echartAPI.value(1) <= data[params.dataIndex - 1][2]) {
          if (customDebug) console.log('Overlap:', echartAPI.value(1), data[params.dataIndex - 1][2]);
          overlap = true;
        }
        const categoryIndex = echartAPI.value(0);
        const start = echartAPI.coord([echartAPI.value(1), categoryIndex + 0]);
        const end = echartAPI.coord([echartAPI.value(2), categoryIndex]);
        if (customDebug) console.log('index start end ', categoryIndex, start, end);
        const height = echartAPI.size([0, 1])[1];
        if (customDebug) console.log('height: ', height);
        const rect0 = {
          x: start[0],
          y: start[1] - this.timeLineWidth / 2,
          width: end[0] - start[0],
          height: this.timeLineWidth,
        };
        const rect1 = {
          x: params.coordSys.x,
          y: params.coordSys.y, // + gridTop,
          width: params.coordSys.width,
          height: params.coordSys.height,
        };
        if (customDebug) console.log('rect0, rect1 ', rect0, rect1);
        let rectShape = echarts.graphic.clipRectByRect(rect0, rect1);
        if (!rectShape) {
          if (customDebug) console.error('no rect shape');
          rectShape = {x: -220, y: 20, width: 100, height: 10};
        }
        let endColor = this.categoryColors[categoryIndex % this.categoryColors.length];
        if (this.selection.includes(params.dataIndex)) endColor = 'red';
        let group = {
          type: 'group',
          children: [{
            type: 'rect',
            transition: ['shape'],
            // shape: rectShape,
            shape: rectShape,
            style: {fill: echartAPI.value(3)},
          },
          {
            type: 'circle',
            shape: {
              cx: start[0], cy: end[1] + 0, r: this.timeLineRadius,
            },
            style: {fill: this.categoryColors[categoryIndex % this.categoryColors.length]},
          },
          {
            type: 'circle',
            shape: {
              cx: end[0], cy: end[1] + 0, r: this.timeLineRadius,
            },
            style: {fill: endColor},
          },
          ],
        };
        if (overlap && (end[0] - start[0] > 20)) {
          // let shift = count ? 20 : -20;
          const shift = (count % 3) ? (count % 3 === 2) ?
            0 : this.timeLineOverlapShift - height / 2 : height / 2 - this.timeLineOverlapShift;
          //    shift = 0;
          count += 1;
          rectShape.y += shift;
          group = {
            type: 'group',
            children: [{
              type: 'rect',
              transition: ['shape'],
              shape: rectShape,
              style: {fill: echartAPI.value(3)},
            },
            {
              type: 'circle',
              shape: {
                cx: start[0], cy: end[1] + shift, r: this.timeLineRadius,
              },
              style: {fill: 'darkgreen'},
            },
            {
              type: 'circle',
              shape: {
                cx: end[0], cy: end[1] + shift, r: this.timeLineRadius,
              },
              style: {fill: 'red'},
            },
            ],
          };
        } // overlap
        if (customDebug) console.log('group: ', group);
        return group;
      })(this, echarts), // render item
      encode: {
        x: [1, 2],
        y: 0,
        tooltip: 3,
      },
      data: this.timeLinesData,
    }; // this.timeLineSeries
    return timeLinesSeries;
  } // initTimeLine
}
