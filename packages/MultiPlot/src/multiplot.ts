

/*
updateHeight: without touch data
updatePlots: complete update
*/

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
//import * as deb from "./../debug.js";

export class MultiPlotViewer extends DG.JsViewer {

  private defaultTitleHeight: any = this.string('defaultTitleHeight', '25px');
  private typeComboElements = [];
  private showHideElements = [];
  private closeElements = [];
  private plots = [];
  private box: any;
  private echart: any;
  private echartOptions: any;
  private visiblePlotsCount: number;
  private timeLineSeries: any;
  private timeLineIndex: number = 1;
  private statusChartIndex: number = 0;
  private timeLineOverlapShift = 4;
  private circleRange: number = 3;
  private timeLineWidth: number = 1;
  private timeLinesData: any = [];
  private statusChartData: any = [];
  private categoryColors: any = [];
  private currentRowCount: number = 0;
  private selection: any = [];
  private plotTitleHighMargin: number = 10;
  private plotTitleLowMargin: number = 5;
  private controlsTopShift: number = -4;

 
  constructor() {
    super();
    console.log('------------------------------------- MULTIPLOT ------------------------------');
    console.log('this.root: ', this.root);
    this.categoryColors = DG.Color.categoricalPalette.map(DG.Color.toRgb);

    this.plots = [
      { height: '1flex', title: 'Events', show: 1, MPtype: 'scatter' },
      { height: '2flex', title: 'Cars length(wheel.base)', show: 1, MPtype: 'timeLine' },
      { height: '80px', title: 'Lab Chemistry', show: 1, MPtype: 'scatter' }
    ];

    var allTypes = ['scatter', 'line', 'bar', 'timeLine'];
    // init plots with some data y=x^2, later can be replaced
    this.plots.map((e, ip) => {
      var r = [];
      for (var i = 0; i < (ip / 1000 + 5) * 2; i++) {
        r.push([i * 20, i * i * 1]);
      }
      e.series = {
        type: 'scatter',
        data: r
      };
    });

    this.initTimeLine();
    this.plots[this.timeLineIndex].series = this.timeLineSeries;
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
      }
    });
  } // ctor 

  setEchartOptions() {

    console.log('set echart options: ', this.echartOptions);

    if (this.echart) {
      // this.echart.clear();
      // this.echart.setOption(this.echartOptions, true);
      setTimeout(() => this.echart.setOption(this.echartOptions, true), 500);
    }
  }

  onPropertyChanged(property) {
    var name = property.name;
    var val = property.get();
    if (name === 'defaultTitleHeight') {
      this.defaultTitleHeight = val;
      this.updateHeight();
    }
    //super.onPropertyChanged(property);
  }

  // input: dataFrame.column
  // output: array of normalized numbers from 0 to 100
  normalize100(column) {
    var r = [];
    var d = column.max - column.min;
    var rawData = column.getRawData();

    for (var i = 0; i < column.length; i++) {
      // avoid case when max == min
      var t = d ? (rawData[i] - column.min) * 100 / d : rawData[i];
      r.push(t);
    }
    return r;
  }

  // input: [{height: '10px', title: 'title0'},{height: '20%'}, {height: 'free'}]
  // output: [{top: '20px', height: '10px', topTitle: '0px', titleHeight}, ...]
  parseHeight(height) {
    var r = [];
    var floatHeight = height;
    var totalFlexPoints = 0;

    for (var i = 0; i < this.plots.length; i++) { // first cycle
      let titleHeight = this.plots[i].title ? parseFloat(this.defaultTitleHeight) : 0;
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
        titleText: this.plots[i].title
      });
      floatHeight -= plotHeight + titleHeight + 
        1 * (this.plotTitleHighMargin + this.plotTitleLowMargin);
    } // first cycle

    console.log('r0 json ', JSON.stringify(r));

    // flexPointHeight - weight (in pixels) of one flex point
    let flexPointHeight = floatHeight / totalFlexPoints;

    var currentTop = 0;
    var echartIndex = 0; // index in echarts library
    for (var j = 0; j < r.length; j++) {
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
    console.log('currentTop  :', currentTop, r);
    return r;
  } // parseHeight

  // output: filled echart options array only for height and position
  updateOptionsPositions() {
    this.box = this.root ? this.root.getBoundingClientRect() : { height: 300 };
    var heightData = this.parseHeight(this.box.height - 30);
    console.log('heightData:  ', heightData);

    var visibleIndex = 0;
    for (var i = 0; i < this.plots.length; i++) {
      this.echartOptions.title[i].top = (heightData[i].titleTop) + 'px';
      this.echartOptions.title[i].left = '10px';
      this.echartOptions.title[i].text = heightData[i].titleText;
      this.echartOptions.title[i].textStyle = {
        height: heightData[i].titleHeight,
        'fontSize': (this.defaultTitleHeight * .6)
      };
      if (this.typeComboElements[i]) {
        this.typeComboElements[i].style.top = heightData[i].titleTop + this.controlsTopShift +'px';
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

  }  // updateOptionsPositions

  // input: this.plots
  // output: filled options for echart library
  updatePlots() {
    this.clearPlots();

    // update positions
    console.log('update plots positions', this.plots);
    this.updateOptionsPositions();
    this.echartOptions.series = [];
    var visibleIndex = 0;
    for (var i = 0; i < this.plots.length; i++) {
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
          { color: this.plots[i].series.itemStyle } : {},
        yAxis: { type: this.plots[i].yType || 'value' },
        xAxisIndex: visibleIndex,
        yAxisIndex: visibleIndex,
        data: this.plots[i].series.data,
        coordinateSystem: 'cartesian2d',
        encode: { x: 0, y: 1 }
      }

      if (this.plots[i].MPtype === 'timeLine') {
        currentSeries = this.plots[i].series;
        currentSeries.xAxisIndex = visibleIndex;
        currentSeries.yAxisIndex = visibleIndex;
        currentSeries.gridIndex = visibleIndex;
        this.echartOptions.xAxis[visibleIndex].type = 'value';
        //this.echartOptions.xAxis[visibleIndex].type = 'value';
        this.echartOptions.yAxis[visibleIndex].type = 'category';
      }

      this.echartOptions.series.push(currentSeries);
      visibleIndex++;
    } // for i<this.plots.length

    this.echartOptions.xAxis[visibleIndex - 1].axisTick = {
      inside: true,
      length: 1000,
      lineStyle: { color: '#cccccc' }
    };
    this.echartOptions.xAxis[visibleIndex - 1].show = true;
    this.echartOptions.xAxis[visibleIndex - 1].type = 'value';
    this.echartOptions.responsive = false;

    console.log('echart options: ', this.echartOptions);
  } // updatePlots

  // only updates heights
  updateHeight() {
    this.updateOptionsPositions();
    this.setEchartOptions();
  }

  clearPlots() {
    this.visiblePlotsCount = 0; // number of visible charts
    for (var i = 0; i < this.plots.length; i++) {
      if (this.plots[i].show) this.visiblePlotsCount++;
    }
    function createEmptyObjects(n) {
      var r = [];
      for (var i = 0; i < n; i++) r.push({});
      return r;
    }

    this.echartOptions = {
      title: createEmptyObjects(this.plots.length),
      'tooltip': {
        trigger: 'axis',
        axisPointer: { type: 'shadow' }
      },
      grid: createEmptyObjects(this.visiblePlotsCount),
      xAxis: createEmptyObjects(this.visiblePlotsCount),
      yAxis: createEmptyObjects(this.visiblePlotsCount),
      dataZoom: [{
        type: 'inside',
        xAxisIndex: Array.from(Array(this.visiblePlotsCount).keys()),
        start: 0,
        end: 100
      }],
      animation: false,
      series: createEmptyObjects(this.visiblePlotsCount)
    };
    for (var i = 0; i < this.echartOptions.series.length; i++) {
      this.echartOptions.series[i].type = 'scatter';
      this.echartOptions.series[i].data = [];
    }

    console.log('clear plots end', this.echartOptions);
  } // clearPlots

  onEvent(e) {
    console.log('event: ', e);
  }

  addMenu() {
    grok.events.onContextMenu.subscribe((args) => {
      console.log('args: ', args.args);
      //    if (!(args.args.context instanceof DG.Viewer)) { return 0; };

      // get opened tables (names, tabs);
      var tabs = {}
      var names: string[] = grok.shell.tables.map(t => t.name);
      for (var i = 0; i < names.length; i++) {
        tabs[names[i]] = grok.shell.tables[i];
      }

      var callback = (item) => {
        console.log('item ', item);
        var table = tabs[item];
        var nCols = Array.from(table.columns.numerical);
        var colNames = nCols.map((e: { name: string }) => e.name);
        var r0raw = table.getCol(colNames[0]).getRawData();
        var column0 = this.normalize100(table.getCol(colNames[0]));
        var column1 = table.getCol(colNames[1]).getRawData();
        var data = [];
        for (var i = 0; i < column0.length; i++) {
          data.push([column0[i], column1[i]]);
        };
        this.plots.push({
          height: '1flex', title: item + ' ' + colNames[1] + '( ' + colNames[0] + ' )',
          series: {
            type: 'scatter',
            large: true,
            data: data
          },
          show: 1
        });
        this.createElements();
        this.updatePlots();
        this.setEchartOptions();
      }; // callback

      // add context menu items (right click);
      let menu = args.args.menu.group('New plot');
      menu.items(names, callback);

    }); // onContextMenu subscribe
  }

  getStatusChartData(table, indexes) {
    var rez = [];
    for (var i = 0; i < table.rowCount; i++) {
      var row = table.row(i);

      var valTestName = row['LBTEST'];
      var valStartTime = row['AESTDY'];

      rez.push([valStartTime, valTestName]);
    }
    this.statusChartData = rez;
  }

  /*  build 2d array for series.data of echart
    @param (table) 
    @param (fieldsNames) 
    @param (indexes) : Array<int>
    @returns (Array<Row>)  */
  getUniversalData(table, fieldsNames, indexes): any[] {
    var r = [];

    function getRowFields(row) {
      var fields = [];
      for (var i = 0; i < fieldsNames.length; i++) {
        fields.push(row[fieldsNames[i]]);
      }
      return fields;
    };

    if (indexes) {
      for (var ind = 0; ind < indexes.length; ind++) {
        let row = table.row(indexes[ind]);
        let fields = getRowFields(row);
        r.push(fields);
      }
    } else {
      for (var i = 0; i < table.rowCount; i++) {
        let row = table.row(i);
        let fields = getRowFields(row);
        r.push(fields);
      }
    }
    return r;
  }


  onTableAttached() {
    if (this.dataFrame.name != 'ae__2__lb__2_') return 0;
    this.addMenu();
    console.log('table attached');
    console.log('this ', this);
    console.log('THIS:PLOTS: ', this.plots);
    var colTimeLineStart: any;
    var colTimeLineEnd: any;

    // load timelines
    for (let tab of grok.shell.tables) {
      if (tab.name === 'ae__2__lb__2_') {
        this.plots[this.statusChartIndex].series.itemStyle = (e, i) => {
          var color = this.categoryColors[0];
          if (this.selection.includes(e.dataIndex)) color = this.categoryColors[1];
          return color;
        }
        this.plots[this.statusChartIndex].series.symbol = (e) => {
          return (e[0] > 22) ? "square" : "triangle";
        }
        this.plots[this.statusChartIndex].yType = 'category';
        console.log('update Filter');
        this.updateFilter(tab);
      }
    } // for tables

    // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    // this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

    // @ts-ignore
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => {
      this.selection = this.dataFrame.selection.getSelectedIndexes();
      this.updatePlots();
      this.setEchartOptions();
      console.log(this.selection);
      this.render();
    }));

    // @ts-ignore
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => {
      console.log('filter change');
      this.updateFilter(0);
      this.setEchartOptions();

      this.render();
    }));

  } // table attached

  updateFilter(tab) {
    if (this.dataFrame.name != 'ae__2__lb__2_') {
      console.error('wrong table');
      return 0;
    }
    var table = tab || this.dataFrame;
    if (table.name != 'ae__2__lb__2_') {
      debugger
    }
    var testNameFieldName = 'LBTEST';
    var startTimeFieldName = 'AESTDY';
    var endTimeFieldName = 'LBDY'
    var testValueFieldName = 'LBSTRESN';
    var testLoFieldName = 'LBSTNRLO';
    var testHiFieldName = 'LBSTNRHI';
    var indexes = table.filter.getSelectedIndexes();
    this.statusChartData = this.getUniversalData(
      table, [startTimeFieldName, endTimeFieldName], indexes
    );
    this.plots[this.statusChartIndex].series.data = this.statusChartData;

    //    return 0;
    this.timeLinesData = this.getUniversalData(
      table, [testNameFieldName, startTimeFieldName, endTimeFieldName], indexes
    );
    this.plots[this.timeLineIndex].series.data = this.timeLinesData;

    this.updatePlots();
    // console.log(new Error());
    //  this.setEchartOptions();
  }

  applyFilter(filter) {
    this.clearPlots();

  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  init() {
    if (!this.echart) this.echart = echarts.init(this.root, null, { renderer: 'canvas' });
    this.clearPlots();
    this.createElements();
    this.updateOptionsPositions();
    this.updatePlots();
    //  this.setEchartOptions();
    this.render();
  }

  deleteElements() {
    this.typeComboElements.map(e => {
      e.remove();
    });
    this.showHideElements.map(e => {
      e.remove();
    });
  }

  createElements() {
    this.deleteElements();

    // create comboboxed to choose plot types
    this.typeComboElements = [];
    for (var i = 0; i < this.plots.length-22; i++) {
      var inputPlotType: any = ui.choiceInput('', 'scatter', ['scatter', 'line', 'bar'], ((i) => (event) => {
        console.log('changed ', event, i);
        this.plots[i].series.type = event;
        this.updatePlots();
        this.setEchartOptions();
      })(i));

 

      this.typeComboElements.push(inputPlotType.root);
      this.root.appendChild(inputPlotType.root);
      inputPlotType.root.style.position = 'absolute';
      inputPlotType.root.style.right = "38px";
      inputPlotType.root.style['flex-direction'] = 'row';
      inputPlotType.root.style.top = (40 * i) + 'px';
    }


    // create checkboxes for show/hide plots
    this.showHideElements = [];
    for (var i = 0; i < this.plots.length; i++) {
        var inputPlotType: any = ui.div([ui.iconFA('angle-right'), ui.iconFA('angle-down')])
        var showHideIcons = inputPlotType.querySelectorAll('i');
      //  debugger
      showHideIcons[0].style.display = 'none';
        inputPlotType.DGswitch = 1;
      inputPlotType.addEventListener('click', ((i) => (e) => {
        var div = e.target.parentNode;
        console.log('click ', i, e.target.parentNode.DGswitch);
        e.target.parentNode.DGswitch = 1 - e.target.parentNode.DGswitch;
        var displays = ['', 'none'];
        var els = div.querySelectorAll('i');
        var sw = e.target.parentNode.DGswitch;
        els[0].style.display = displays[sw];
        els[1].style.display = displays[1 - sw];

        this.plots[i].show = sw;
        this.updatePlots();
        this.setEchartOptions();
      })(i))
        
      this.typeComboElements.push(inputPlotType);
      this.root.appendChild(inputPlotType);
      inputPlotType.style.position = 'absolute';
      inputPlotType.style.left = "3px";
      inputPlotType.style['flex-direction'] = 'row';
      inputPlotType.style.top = (40 * i) + 'px';
      this.showHideElements.push(inputPlotType);
    }




 /*
    // create checkboxes for show/hide plots
    this.showHideElements = [];

    for (var i = 0; i < this.plots.length; i++) {
      let inputShowHide: any = ui.boolInput('', true, ((i) => (event) => {
        console.log('checked: ', event, i);
        this.plots[i].show = event;
        this.updatePlots();
        this.setEchartOptions();
      })(i));

      this.root.appendChild(inputShowHide.root);
      inputShowHide.root.style.position = 'absolute';
      inputShowHide.root.style.right = "26px";
      inputShowHide.root.style.top = (40 * i) + 'px';
      inputShowHide.root.style.flexDirection = 'row';
      inputShowHide.root.style.topMargin = '5px';

      this.showHideElements.push(inputShowHide.root);
      let inp = inputShowHide.root.querySelector('input');
      inp.style.width = '12px';
      inp.style.minWidth = '12px';
      inp.style.opacity = '1';

      //     this.root.appendChild(inputShowHide.root);
    }
*/

    // create close 'X' icons
    this.closeElements = [];
    for (var i = 0; i < this.plots.length; i++) {
      var inputClose: any = ui.icons.close(((i) => () => {
        grok.shell.info('click' + i);
        this.plots.splice(i, 1);
        this.updatePlots();
        this.setEchartOptions();
      })(i), 'Close');
      inputClose.style.position = 'absolute';
      inputClose.style.right = "12px";
      inputClose.style.top = (40 * i) + 'px';
      inputClose.style.flexDirection = 'row';
      this.closeElements.push(inputClose);

      this.root.appendChild(inputClose);
    } // close X

  } // createElements

  render() {

  }

  initTimeLine() {
    console.error('init time line :', this.timeLinesData);
    let thisData = this.plots[this.timeLineIndex].series.data;
    let data = this.timeLinesData;
    var t = this;
    let count = 0;
    this.timeLineSeries = {
      type: 'custom',
      progressive: 0,
      animation: false,
      renderItem: ((t: any, echarts: any) => (params, api) => {
        let data = this.timeLinesData;
        if (!data || data.length == 0) return 1234;

        var customDebug = false;

        if (customDebug) console.log('renderItem data: ', data);
        if (customDebug) console.log('custom render ', params, api);
        var av0 = api.value(0);
        var av1 = api.value(1);
        var av2 = api.value(2);
        if (customDebug) console.log('values: ', av0, av1, av2);

        var gridTopRaw = this.echartOptions.grid[this.timeLineIndex].top;
        var gridTop = parseFloat(gridTopRaw);
        if (customDebug) console.log('gridtop ', gridTop);
        let overlap = false;

        if (params.dataIndex > 0 && data[params.dataIndex - 1][0] === data[params.dataIndex][0]
          && api.value(1) <= data[params.dataIndex - 1][2]) {
          if (customDebug) console.log('Overlap:', api.value(1), data[params.dataIndex - 1][2]);
          overlap = true;
        }

        var categoryIndex = api.value(0);
        var start = api.coord([api.value(1), categoryIndex + 0]);
        var end = api.coord([api.value(2), categoryIndex]);
        if (customDebug) console.log('index start end ', categoryIndex, start, end);

        var height = api.size([0, 1])[1];
        if (customDebug) console.log('height: ', height);

        var rect0 = {
          x: start[0],
          y: start[1] - this.timeLineWidth / 2,
          width: end[0] - start[0],
          height: this.timeLineWidth,
        };
        var rect1 = {
          x: params.coordSys.x,
          y: params.coordSys.y,// + gridTop,
          width: params.coordSys.width,
          height: params.coordSys.height
        };
        if (customDebug) console.log('rect0, rect1 ', rect0, rect1);
        var rectShape = echarts.graphic.clipRectByRect(rect0, rect1);
        if (!rectShape) {
          if (customDebug) console.error('no rect shape');
          rectShape = { x: -220, y: 20, width: 100, height: 10 };
        }
        var endColor = this.categoryColors[categoryIndex % this.categoryColors.length];
        if (this.selection.includes(params.dataIndex)) endColor = 'red';
        let group = {
          type: 'group',
          children: [{
            type: 'rect',
            transition: ['shape'],
            //shape: rectShape,
            shape: rectShape,
            style: { fill: api.value(3) }
          },
          {
            type: 'circle',
            shape: {
              cx: start[0], cy: end[1] + 0, r: this.circleRange
            },
            style: { fill: this.categoryColors[categoryIndex % this.categoryColors.length] }
          },
          {
            type: 'circle',
            shape: {
              cx: end[0], cy: end[1] + 0, r: this.circleRange
            },
            style: { fill: endColor }
          }
          ]
        };
        if (overlap && (end[0] - start[0] > 20)) {
          //let shift = count ? 20 : -20;
          let shift = (count % 3) ? (count % 3 === 2) ?
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
              style: { fill: api.value(3) }
            },
            {
              type: 'circle',
              shape: {
                cx: start[0], cy: end[1] + shift, r: this.circleRange
              },
              style: { fill: 'darkgreen' }
            },
            {
              type: 'circle',
              shape: {
                cx: end[0], cy: end[1] + shift, r: this.circleRange
              },
              style: { fill: 'red' }
            }
            ]
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
      data: this.timeLinesData
    } // this.timeLineSeries
  } // initTimeLine
}
