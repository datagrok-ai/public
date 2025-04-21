import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
import {MPUtils} from './utils';
import {MPLayout} from './layout';
import { BAR, LINE, SCATTER, TIMELINES } from './constants';
import '../../../css/multiplot.css';
// import * as deb from "./../debug.js";

const startCol = 0;
const endCol = 1;
const categoryCol = 2;

@grok.decorators.viewer({
  name: 'Multiplot',
  description: 'Creates a multiplot viewer',
})
export class MultiPlotViewer extends DG.JsViewer {
  // properties
  private defaultTitleHeight: string = this.string('defaultTitleHeight', '25px');
  private backColor: number = this.int('backColor', 0xffffff);
  private verticalLinesColor: number = this.int('verticalLinesColor', 0xf0f0f0);
  private verticalLines: boolean = this.bool('verticalLines', true);
  private showControls: boolean = this.bool('showControls', true);
  private timeLineWidth: number = this.int('timeLineWidth', 1);
  private timeLineRadius: number = this.int('timeLineRadius', 3);
  private marker = this.string('marker', 'ring', {choices: ['circle', 'rect', 'ring', 'diamond']});
  private markerSize = this.int('markerSize', 6);
  private lineWidth = this.int('lineWidth', 2);
  private markerPosition = this.string('markerPosition', 'main line',
    {choices: ['main line', 'above main line', 'scatter']});

  private echart: echarts.EChartsType | null = null;
  // use 'any' type, because of missing echarts.ECharts type
  private echarts: any = echarts;
  // commented and use 'any' because of missing some members in echarts.EChartsOption
  // private echartOptions: echarts.EChartsOption;
  private echartOptions: any;

  private paramA : string = this.string('paramA', 'string inside');

  private timeLineSeries: any;
  private timeLineIndex: number = 11;
  private statusChartIndex: number = 10;
  private timeLineOverlapShift = 4;
  private timeLinesData: any = [];

  private paletteColors: string[] = [];
  private plotTitleHighMargin: number = 10;
  private plotTitleLowMargin: number = 5;
  private controlsTopShift: number = -4;
  private headerHiddenShift: number = -10;
  private tooltipOffset: number = 10;
  private globalMarginTop: number = 15;
  private globalMarginBottom: number = 0;
  private selectionColor : DG.Color = DG.Color.toRgb(DG.Color.selectedRows);
  private categoryLength : number = 9;
  private paramOptions : string = this.string('paramOptions', 'none22');
  private mode : string = 'none'; // 'brushSelected'
  private options: {series: any[], xAxisMinMax: {[key: string]: any}} = {series: [], xAxisMinMax: {}};

  typeComboElements : (HTMLElement | null)[] = [];
  showHideElements : HTMLElement[] = [];
  closeElements : HTMLElement[] = [];
  multiEditElements: (HTMLElement | null)[] = [];
  comboEditElements: (HTMLElement | null)[] = [];
  categoryCombos : HTMLElement[] = [];
  plots: any[] = [];
  box: DOMRect;
  visiblePlotsCount: number = 0;
  visibleIndexes: number[] = [];
  tables: {[key: string]: any} = {};
  isTablesLoaded: number = 0;
  isEchartHandlers : boolean = false;
  selectionIndexes: number[] = [];
  utils : any = new MPUtils();
  layout: any = new MPLayout();
  allTypes = ['scatter', 'line', 'bar', 'timeLine'];

  constructor() {
    super();
    console.log('------------------------------------- MULTIPLOT ------------------------------');
    console.log('this.root: ', this.root);

    this.paletteColors = DG.Color.categoricalPalette.map(DG.Color.toRgb);
    // this.plots = this.options.series;
    // this.init();

    this.box = this.root.getBoundingClientRect();
    ui.onSizeChanged(this.root).subscribe(((this2) => (e: any) => {
      if (!this2.echartOptions) return 0;
      this2.updateOptionsPositions();
      if (this2.echart) {
        this2.updatePlots();
        this2.render();
        this2.echart.resize();
      }
    })(this));
  }

  init(): void {
    console.warn('init');
    if (!this.isTablesLoaded) return;
    this.clearPlots();
    this.plots.map((e) => {
      if (e.series) {
        if (!e.series.type) e.series.type = 'scatter';
      } else
        e.series = {type: 'scatter'};

      e.series.data = [];
      e.selectedIndexes = [];
    });
    if (!this.echart) this.echart = echarts.init(this.root, null, {renderer: 'canvas'});
    this.addEchartHandlers();
    this.updateOptionsPositions();
    this.updatePlots();
    if (this.checkTablesLoaded())
      this.updateFilter();

    this.createElements();

    this.render();
  } // init

  onPropertyChanged(property: DG.Property): void {
    if (!Object.keys(this.tables).length)
      return;
    ;
    const name = property.name;
    const val = property.get(this);
    console.log('property changed: ', name, val);
    if (name === 'defaultTitleHeight') {
      this.defaultTitleHeight = val;
      this.updateOptionsPositions();
      this.render();
    }
    if (name === 'backColor') {
    }
    if (name === 'verticalLines') {
    }
    if (name === 'showControls')
      this.setControlsVisibility();

    // import parameters from external script
    if (name === 'paramOptions') {
      if (val === 'none') return;
      const param = JSON.parse(this.paramOptions);
      if (param.series.length === 0) return;
      this.plots = param.series;


      this.plots = this.utils.getPlotsFromParams(this.tables, param.series);
      this.init();
      return;
    }

    this.updatePlots();
  //  this.render();
    // super.onPropertyChanged(property);
  }

  // output: filled echart options array only for height and position
  updateOptionsPositions(): void {
    this.box = this.root.getBoundingClientRect();
    const heightData = this.layout.parseHeight(this.box.height - 30, this.plots, this);
    let visibleIndex = 0;
    for (let i = 0; i < this.plots.length; i++) {
      this.echartOptions.title[i].top = (heightData[i].titleTop) + 'px';
      this.echartOptions.title[i].left = '10px';
      this.echartOptions.title[i].text = heightData[i].titleText;
      this.echartOptions.title[i].textStyle = {
        'height': heightData[i].titleHeight,
        'fontSize': (parseFloat(this.defaultTitleHeight) * .6),
      };
      this.echartOptions.legend[i].top = (heightData[i].titleTop) + 'px';

      if (this.typeComboElements[i])
        this.typeComboElements[i]!.style.top = heightData[i].titleTop + this.controlsTopShift + 'px';


      if (this.showHideElements[i])
        this.showHideElements[i].style.top = heightData[i].titleTop + 7 + this.controlsTopShift + 'px';


      if (this.closeElements[i])
        this.closeElements[i].style.top = heightData[i].titleTop + 10 + this.controlsTopShift + 'px';


      if (this.multiEditElements[i])
        this.multiEditElements[i]!.style.top = heightData[i].titleTop + 7 + this.controlsTopShift + 'px';


      if (this.comboEditElements[i])
        this.comboEditElements[i]!.style.top = heightData[i].titleTop + this.controlsTopShift + 'px';


      if (this.plots[i].categCombo)
        this.plots[i].categCombo.root.style.top = heightData[i].titleTop + 3 + this.controlsTopShift + 'px';


      if (!this.plots[i].show) continue;

      // only for visible plots
      // this.echartOptions.grid[visibleIndex].y2 = 22;
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
    if (this.isTablesLoaded === 0) return;

    this.echartOptions.backgroundColor = this.utils.toStringColor(this.backColor);

    // update positions
    this.updateOptionsPositions();
    this.echartOptions.series = [];
    let visibleIndex = 0;
    this.visibleIndexes = [];
    this.echartOptions.visualMap = [];
    const leftTitles = [];

    for (let i = 0; i < this.plots.length; i++) {
      const plot = this.plots[i];
      this.echartOptions.title[i].left = '10px';
      this.echartOptions.title[i].text = this.plots[i].title;
      this.echartOptions.legend[i].data = plot.currentCat;
      this.echartOptions.legend[i].type = 'scroll';
      this.echartOptions.legend[i].show = !!plot.showLegend;
      this.echartOptions.legend[i].width = '400px';

      if (!this.plots[i].show) continue;
      if (plot.leftTitle) {
        leftTitles.push({
          left: '3px',
          top: parseFloat(this.echartOptions.title[i].top)+ plot.floatHeight / 2 - 6,
          text: plot.condition.value,
          textStyle: {fontSize: 12, height: 8, width: 50, overflow: 'truncate'},
        });
      }
      this.echartOptions.grid[visibleIndex].left = '82px';
      this.echartOptions.grid[visibleIndex].right = '4%';
      this.echartOptions.grid[visibleIndex].show = this.plots[i].show;
      this.echartOptions.xAxis[visibleIndex].gridIndex = visibleIndex;
      this.echartOptions.xAxis[visibleIndex].type = 'value';
      //    this.echartOptions.xAxis[i].axisTick = { inside: true, length: 1000 };
      this.echartOptions.xAxis[visibleIndex].show = false;
      this.echartOptions.xAxis[visibleIndex].min = this.options.xAxisMinMax['minX'];
      this.echartOptions.xAxis[visibleIndex].max = this.options.xAxisMinMax['maxX'];
      this.echartOptions.yAxis[visibleIndex].gridIndex = visibleIndex;
      this.echartOptions.yAxis[visibleIndex].type = this.plots[i].yType || 'value';
      this.echartOptions.yAxis[visibleIndex].show = this.plots[i].show;
      this.echartOptions.yAxis[visibleIndex].triggerEvent = true;
      this.echartOptions.yAxis[visibleIndex].axisLabel = {width: plot.yLabelWidth};
      this.echartOptions.yAxis[visibleIndex].axisLabel.overflow = plot.yLabelOverflow;
      this.echartOptions.yAxis[visibleIndex].axisLine = {onZero: false};
      this.echartOptions.yAxis[visibleIndex].axisTick = {alignWithLabel: true};

      if (plot.multiLineFieldIndex && plot.series.data.length) {
        let minY: number;
        let maxY: number;
        plot.series.data.forEach((item: any, index: number) => {
          this.echartOptions.series
            .push(this.getCurrentSeries(plot, visibleIndex, i, plot.series.type, plot.currentCat[index], item));
          item.forEach((it: any) => {
            if (!minY || minY > it[1]) minY = it[1]; ;
            if (!maxY || maxY < it[1]) maxY = it[1]; ;
          });
        });
        this.echartOptions.yAxis[visibleIndex].min = plot.additionalParams ?
          plot.additionalParams.min && plot.additionalParams.min < minY! ?
            plot.additionalParams.min : minY! : undefined;
        this.echartOptions.yAxis[visibleIndex].max = plot.additionalParams ?
          plot.additionalParams.max && plot.additionalParams.max > maxY! ?
            plot.additionalParams.max : maxY!: undefined;
      } else {
        this.echartOptions.series.push(
          this.getCurrentSeries(plot, visibleIndex, i, plot.series.type, plot.currentCat, plot.series.data));
      }
      this.visibleIndexes.push(i);
      visibleIndex++;
    } // for i<this.plots.length

    this.echartOptions.title = this.echartOptions.title.concat(leftTitles);

    // this.echartOptions.responsive = false;
    if (visibleIndex === 0)
      return;


    this.echartOptions.xAxis[visibleIndex - 1].axisTick = {
      inside: true,
      length: this.verticalLines ? 2000 : 5,
      lineStyle: {color: this.utils.toStringColor(this.verticalLinesColor)},
    };
    this.echartOptions.xAxis[visibleIndex - 1].show = true;
    this.echartOptions.xAxis[visibleIndex - 1].type = 'value';
  } // updatePlots


  getCurrentSeries(plot: any, visibleIndex: number, i: number, type: string, name: string, data: any) {
    let currentSeries = {
      name: name,
      type: type,
      show: plot.show,
      large: true,
      gridIndex: visibleIndex,
      xAxisIndex: visibleIndex,
      yAxisIndex: visibleIndex,
      data: data,
      coordinateSystem: 'cartesian2d',
      encode: { x: 0, y: 1 },
      selectedMode: 'multiple',
      xAxis: {},
      yAxis: {
        type: plot.yType,
      },
      itemStyle: {},
      markLine: plot.additionalParams ? plot.additionalParams.markLine : undefined,
    };

    if (plot.statusChart) {
      currentSeries.itemStyle = {
        color: this.getItemStyleColorFunc(plot.visualMap, plot),
      };
    }

    if (plot.visualMap) {
      if (plot.visualMap.pieces) {
        currentSeries.itemStyle = {
          color: this.getItemStyleColorFunc(plot.visualMap, plot),
        };
      } else {
        // if (plot.visualMap) {
        // keep it to find is any execution ever happens here
        debugger;
        const map = plot.visualMap;
        if (map.type === 'piecewise') {
          const min = map.pieces[0].min;
          const max = map.pieces[0].max;
          map.pieces.push({ max: min, color: this.paletteColors[0] });
          map.pieces.push({ min: max, color: this.paletteColors[0] });
        }
        map.seriesIndex = visibleIndex;
        this.echartOptions.visualMap.push(map);
        map.show = false;
        // }
      } // if visualMap.pieces
    }

    if (plot.type === 'timeLine') {
      plot.timeLinesSeries = this.initTimeLine(i, this);
      currentSeries = plot.timeLinesSeries;
      currentSeries.xAxisIndex = visibleIndex;
      currentSeries.yAxisIndex = visibleIndex;
      currentSeries.gridIndex = visibleIndex;
      // currentSeries.encode = {x: [1, 2], y: 0};
      this.echartOptions.xAxis[visibleIndex].type = 'value';
      this.echartOptions.yAxis[visibleIndex].type = 'category';
    }

    return currentSeries;
  }


  // create function to use as EChart callback with Datagrok mixins
  // get callback function to define color of marker
  getItemStyleColorFunc(visualMap: any, plot: any) : any {
    let customColorFunc : any = () => {};
    const defaultColor = this.paletteColors[0];
    const selectionColor = this.paletteColors[1];
    const table = this.tables[plot.tableName];
    if (visualMap) {
      const min = visualMap.pieces[0].min;
      const max = visualMap.pieces[0].max;
      const vMapColor = visualMap.pieces[0].color;
      customColorFunc = (e: any) => {
        return e.data[2] > min && e.data[2] < max ? vMapColor : defaultColor;
      };
    }
    if (plot.statusChart) {
      console.error('status chart get color func');
      customColorFunc = (e: any) => {
        let val = e.data[plot.statusChart.valueField];
        let min = e.data[plot.statusChart.minField];
        let max = e.data[plot.statusChart.maxField];
        if (typeof val == 'string') val = parseFloat(val);
        if (typeof min == 'string') min = parseFloat(min);
        if (typeof max == 'string') max = parseFloat(max);

        return val > min && val < max ? 'green' : 'red';
        return val > min && val < max ? defaultColor : plot.statusChart.alertColor;
      };
    }
    const f = (e: any) => {
      const customColor = customColorFunc(e);
      // if (plot.selectedIndexes.includes(e.dataIndex)) {
      if (this.utils.getBitByIndex32(table.selection, e.dataIndex))
        return selectionColor;

      return customColor;
    };
    return f.bind(this);
  }

  // shapes provided by ECharts:
  // 'circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none'
  getMarkerSymbolFunc(customSymbolFunc: any, plot: any) : any {
    const defaultSymbol = 'circle';
    console.warn('getShapeFunc: ', defaultSymbol);
    function f(e: any) {
      // console.log('symbol ', e);
      const customSymbol = customSymbolFunc(e);
      return customSymbol ? customSymbol : defaultSymbol;
    }
    return f;
  }

  // only updates heights
  updateHeight() : void {
    this.updateOptionsPositions();
    this.render();
  }

  // fill echart options with arrays of empty objects
  clearPlots() : void {
    this.visiblePlotsCount = 0; // number of visible charts
    for (let i = 0; i < this.plots.length; i++) {
      if (this.plots[i].show)
        this.visiblePlotsCount++;
    }
    function createEmptyObjects(n: number) : any {
      return new Array(n).fill(null).map(() => {
        return {};
      });
    }

    this.echartOptions = {
      title: createEmptyObjects(this.plots.length),
      legend: createEmptyObjects(this.plots.length),
      tooltip: {
        trigger: 'axis',
        showContent: false,
        axisPointer: {type: 'shadow'},
      },
      textStyle: {
        fontFamily: 'Roboto',
        overflow: 'truncate',
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

      brush: {
        toolbox: ['rect', 'polygon', 'lineX', 'lineY', 'keep', 'clear'],
        xAxisIndex: 'all',
      },

    };
    for (let i = 0; i < this.echartOptions.series.length; i++) {
      this.echartOptions.series[i].type = 'scatter';
      this.echartOptions.series[i].data = [];
    }
  } // clearPlots

  // onEvent(e: DG.Events): void {  }

  checkTablesLoaded() : boolean {
    const plotTablesList : string[] = this.plots.map((e) => e.tableName);
    const openTablesArray: string[] = Object.keys(this.tables);
    const notLoaded : string[] = plotTablesList.filter((e) => openTablesArray.indexOf(e) == -1);
    console.warn('not loaded list', notLoaded);
    return notLoaded.length == 0;
  }

  onTableAttached(): void {
    if (!this.checkTablesLoaded() || Object.keys(this.tables).length === 0)
      return;

    console.warn('all tables loaded');
    this.plots = this.utils.getPlotsFromParams(this.tables, this.options.series);
    this.isTablesLoaded = 1;
    this.init();
    this.updatePlots();
    this.updateFilter();
    this.render();

    // this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    // this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => {
    }));
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => {
      this.updateFilter();
      this.render();
    }));
  } // table attached

  updateFilter(): void {
    this.timeLinesData = [];
    for (let i=0; i< this.plots.length; i++) {
      const plot = this.plots[i];
      const tableName = plot.tableName;
      const table = this.tables[tableName];
      const indexes = table.filter.getSelectedIndexes();
      const x = plot.x;
      const y = plot.y;
      const xArray = Array.isArray(x) ? x : [x];
      const yArray = Array.isArray(y) ? y : [y];
      const visualMapColumnName = [];
      if (plot.visualMap) {
        if (plot.visualMap.column)
          visualMapColumnName.push(plot.visualMap.column);
        else
          visualMapColumnName.push(y);
      }

      const data = this.utils.getUniversalData(
        table,
        // xArray.concat(yArray).concat(visualMapColumnName),
        xArray.concat(yArray).concat(plot.extraFields ?? []),
        indexes,
        plot.condition,
      );
      /*       if (plot.categoryColumnIndex) {
        for (let i=0; i<data.length; i++) {
          data[i][plot.categoryColumnIndex] = this.utils.trimCategoryString(
              data[i][plot.categoryColumnIndex],
              this.categoryLength,
          );
        }
      } */
      //break data into several arrays according to category(required to draw multiple lines on linechart)
      if (plot.multiLineFieldIndex && plot.currentCat) {
        const multipleData: any[] = [];
        plot.currentCat.forEach((cat: any) => {
          const catData: any[] = [];
          data.forEach((item: any) =>{
            if (item[plot.multiLineFieldIndex] === cat)
              catData.push(item);
          });
          multipleData.push(catData);
        });
        plot.series.data = multipleData;
      } else
        plot.series.data = data;


      if (plot.type != 'timeLine') {
        this.timeLinesData.push([]);
        continue;
      };
      this.timeLinesData.push(data);
      // this.plots[i].series.data = data;
      this.timeLineIndex = i;
    }
    this.updatePlots();
  }

  isGroup(componentIndex: number, componentType: string) : boolean {
    const type = this.plots[componentIndex].type;
    const yType = this.plots[componentIndex].yType;

    if ((type === 'scatter' || type === 'line') &&
      yType == 'category' && componentType == 'yAxis') return true;
    if (type === 'scatter' || type === 'line') return false;
    return true;
  }

  addEchartHandlers() : void {
    if (this.isEchartHandlers) return;

    this.echart!.on('dataZoom', (e: any) => {
      console.log('zoom ', e);
    });

    this.echart!.on('brushSelected', (e: any) => {
      this.mode = 'brushSelected';
      console.log('brush selected ', e);
      if (e.batch.length === 0) return;
      const batchSelected = e.batch[0].selected;
      if (batchSelected.length === 0) return;
      const index = batchSelected[0].dataIndex;
      for (let i=0; i<batchSelected.length; i++) {
        const table = this.tables[this.plots[i].tableName];
        table.selection.setAll(0);
        for (let j=0; j<index.length; j++)
          table.selection.set(index[j], true);
      };
      // this.updatePlots();
      // this.render();
    }); // onBrushSelected

    this.echart!.on('brushEnd', (e: any) => {
      console.error('brushEnd, ', e);
    });

    this.echart!.on('mousedown', {datatype: ''}, (e: any) => {
      console.log('mouse down: ', e);
    });

    this.echart!.on('mouseup', {datatype: ''}, (e: any) => {
      console.log('mouse up: ', e);
    });

    this.echart!.on('click', {datatype: 'all'}, (params : any) => {
      console.log('echart click', params);
      const iPlot : number = this.visibleIndexes[params.componentIndex];
      const table : DG.DataFrame = this.tables[this.plots[iPlot].tableName];
      const indexes = table.filter.getSelectedIndexes();
      table.currentRowIdx = indexes[params.dataIndex];
      table.selection.handleClick((i) => {
        if (params.componentType === 'yAxis')
          return this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]] === params.value;

        if (params.componentType === 'series') {
          if (this.isGroup(params.componentIndex, '')) {
            return params.value[0] ===
              this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]];
          } else
            return params.dataIndex === i;
        }
        return false;
      }, params.event.event);
    });

    this.echart!.on('mouseover', (params: any) => {
      const iPlot : number = this.visibleIndexes[this.echartOptions.series[params.componentIndex].gridIndex];
      const x = (params.event.event as MouseEvent).x + this.tooltipOffset;
      const y = (params.event.event as MouseEvent).y + this.tooltipOffset;
      const xColName = this.utils.getArrayOfColumnNames(this.plots[iPlot].x);
      const yColName = this.utils.getArrayOfColumnNames(this.plots[iPlot].y);
      let colNames = xColName.concat(yColName);
      if (this.plots[iPlot].extraFields) colNames = colNames.concat(this.plots[iPlot].extraFields);

      if (params.componentType === 'yAxis') {
        if (this.isGroup(iPlot, params.componentType)) {
          /* ui.tooltip.showRowGroup(table, (i) => {
            return params.value === this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]];
          }, x, y); */
          ui.tooltip.show(ui.div([params.value as string]), x, y);
        }
      }

      if (params.componentType === 'series') {
        if (!this.isGroup(iPlot, '')) {
          ui.tooltip.show(ui.divV(
            (params.data as any[]).map((e, i) => ui.div([colNames[i] + ': ' + e + ''])),
          ), x, y);
        } else {
          /*  ui.tooltip.showRowGroup(table, (i) => {
            return params.value[0] === this.plots[iPlot].subjects[this.plots[iPlot].subjBuf[i]]; // &&
            //            params.value[1] === (this.startCol.isNone(i) ? null : this.startBuf[i]) &&
            //          params.value[2] === (this.endCol.isNone(i) ? null : this.endBuf[i]);
          }, x, y); */
          ui.tooltip.show(ui.divV(
            (params.data as any[]).map((e, i) => ui.div([colNames[i] + ': ' + e + ''])),
          ), x, y);
        }
      } // series
    }); // mouseover

    this.echart!.on('mouseout', () => ui.tooltip.hide());
    this.isEchartHandlers = true;
  }

  deleteElements(): void {
    this.typeComboElements.map((e) => {if (e) e.remove();});
    this.showHideElements.map((e) => e.remove());
    this.closeElements.map((e) => e.remove());
    this.multiEditElements.map((e) => {if (e) e.remove();});
    this.comboEditElements.map((e) => {if (e) e.remove();});
  }

  setControlsVisibility() : void {
    this.typeComboElements.map((e) => {
      if (e)
        e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
    this.showHideElements.map((e) => {
      e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
    this.closeElements.map((e) => {
      e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
    this.multiEditElements.map((e) => {
      if (e)
        e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
    this.comboEditElements.map((e) => {
      if (e)
        e.style.visibility = this.showControls ? 'visible' : 'hidden';
    });
  }

  createElements(): void {
    this.deleteElements();

    // create comboboxes to choose plot types
    this.typeComboElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      if (this.plots[i].type !== TIMELINES) {
        const inputPlotType: DG.InputBase = ui.input.choice('', {
          value: this.plots[i].series.type, items: [SCATTER, LINE, BAR],
          onValueChanged: (event) => {
            this.plots[i].series.type = event;
            this.updatePlots();
            this.render();
          }});
        this.typeComboElements.push(inputPlotType.root);
        this.root.appendChild(inputPlotType.root);
        inputPlotType.root.classList.add('multiplot-plot-type-input');
        inputPlotType.root.style.top = (40 * i) + 'px';
        inputPlotType.root.querySelector('select')!.style.borderBottom = '0px';
      } else
        this.typeComboElements.push(null);
    }

    // create checkboxes for show/hide plots
    this.showHideElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      const inputPlotType: any = ui.div([ui.iconFA('angle-right'), ui.iconFA('angle-down')]);
      const showHideIcons = inputPlotType.querySelectorAll('i');
      showHideIcons[this.plots[i].show ? 0 : 1].style.display = 'none';
      inputPlotType.showSwitch = this.plots[i].show ? 1 : 0;
      inputPlotType.addEventListener('click', (e: any) => {
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
      });
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
      const inputClose: any = ui.icons.close(() => {
        this.plots.splice(i, 1);
        this.createElements();
        // this.updateHeight();
        this.updatePlots();
        this.render();
      }, 'Close');
      inputClose.style.position = 'absolute';
      inputClose.style.right = '15px';
      inputClose.style.top = (40 * i) + 'px';
      inputClose.style.flexDirection = 'row';
      this.closeElements.push(inputClose);
      this.root.appendChild(inputClose);
    }

    // create edit icons
    this.multiEditElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      if (this.plots[i].multiEdit) {
        const inputEdit = ui.icons.settings(() => {
          const valuesMultiChoices = ui.input.multiChoice('', {
            value: this.plots[i].multiEdit.selectedValues,
            items: this.plots[i].multiEdit.options,
            onValueChanged: (v) => {
              this.plots[i].multiEdit.selectedValues = valuesMultiChoices.value;
            },
          });
          //@ts-ignore
          valuesMultiChoices.input.style.maxWidth = '100%';
          //@ts-ignore
          valuesMultiChoices.input.style.maxHeight = '100%';
          ui.dialog({ title: 'Select values' })
            .add(ui.div([valuesMultiChoices], { style: { width: '400px', height: '300px' } }))
            .onOK(() => {
              if (this.plots[i].multiEdit.editValue === 'category') {
                this.updatePlotByCategory(i,
                  this.plots[i].multiEdit.selectedValues, this.plots[i].multiEdit.updateTitle);
              }
            })
            .show();
        }, 'Edit values');
        inputEdit.style.right = '100px';
        inputEdit.style.position = 'absolute';
        inputEdit.style.top = (40 * i) + 'px';
        inputEdit.style.flexDirection = 'row';
        this.root.appendChild(inputEdit);
        this.multiEditElements.push(inputEdit);
      } else
        this.multiEditElements.push(null);
    }


    // create edit comboboxes
    this.comboEditElements = [];
    for (let i = 0; i < this.plots.length; i++) {
      if (this.plots[i].comboEdit) {
        const inputEdit = ui.input.choice(this.plots[i].comboEdit.editName, {
          value: this.plots[i].comboEdit.selectedValue,
          items: this.plots[i].comboEdit.options,
          onValueChanged: (v) => {
            this.plots[i].comboEdit.selectedValue = inputEdit.value;
            this.setAdditionalParams(i,
              this.plots[i].comboEdit.additionalParams[this.plots[i].comboEdit.selectedValue]);
            if (this.plots[i].comboEdit.editValue === 'category')
              this.updatePlotByCategory(i, this.plots[i].comboEdit.selectedValue, this.plots[i].multiEdit.updateTitle);

            if (this.plots[i].comboEdit.editValue === 'y')
              this.updatePlotByYAxis(i, this.plots[i].comboEdit.values[this.plots[i].comboEdit.selectedValue]);
          },
        });
        inputEdit.root.style.top = (40 * i) + 'px';
        inputEdit.root.classList.add('edit-plot-input');
        this.root.appendChild(inputEdit.root);
        this.comboEditElements.push(inputEdit.root);
      } else
        this.comboEditElements.push(null);
    }

    // create combobox with categories
    /*     this.categoryCombos = [];
    for (let i=0; i<this.plots.length; i++) {
      const plot = this.plots[i];
      if (plot.allCats) {
        const categCombo: any = ui.choiceInput('', plot.currentCat, plot.allCats, (e) => {
          plot.currentCat = e;
          plot.condition.value = e;
          this.updateFilter();
          this.render();
        });
        this.categoryCombos.push(categCombo);
        plot.categCombo = categCombo;
        this.root.appendChild(categCombo.root);
        categCombo.root.style.position = 'absolute';
        categCombo.root.style.right = '78px';
        categCombo.root.style['flex-direction'] = 'row';
        categCombo.root.style.top = (40 * i) + 'px';
        categCombo.root.querySelector('select').style.borderBottom = '0px';
      }
    } */
  } // createElements


  updatePlotByCategory(plotIndex: number, category: any, updateTitle: boolean) {
    if (this.plots[plotIndex].type === 'scatter' && category.length)
      this.plots[plotIndex].height = `${category.length*20}px`;

    this.plots[plotIndex].currentCat = category;
    this.plots[plotIndex].condition.value = category;
    if (updateTitle) this.plots[plotIndex].title = category;
    this.updateFilter();
    this.render();
  }

  updatePlotByYAxis(plotIndex: number, yAxis: any) {
    this.plots[plotIndex].y = yAxis;
    this.updateFilter();
    this.render();
  }

  setAdditionalParams(plotIndex: number, params: any) {
    this.plots[plotIndex].additionalParams = params;
  }

  render(): void {
    console.log('set echart options: ', this.echartOptions);
    if (this.echart) {
      this.echart.clear();
      // this.echart.setOption(this.echartOptions, true);
      setTimeout(() => this.echart!.setOption(this.echartOptions, true), 200);
    }
  }


  // initTimeLine

  timeLinesD(params: any, api: any) {
    const overlap = false;

    const categoryIndex = api.value(categoryCol);
    const start = api.coord([api.value(startCol), categoryIndex]);
    const end = api.coord([api.value(endCol), categoryIndex]);
    const width = end[0] - start[0];

    const group: {type: string, children: any[]} = {
      type: 'group',
      children: [],
    };

    // if no reason to draw line
    if (isNaN(api.value(startCol)) || isNaN(api.value(endCol)) || this.markerSize > width) {
      const xPos = (shift: any) => isNaN(start[0]) ? end[0] : start[0] - shift;
      const yPos = (shift: any) => end[1] - (this.markerPosition === 'main line' ? shift :
        this.markerPosition === 'above main line' ? Math.max(this.markerSize, this.lineWidth) + shift :
          ((params.dataIndex % 2) * 2 - 1)*(this.markerSize * 3));

      const marker = {
        type: this.marker,
        shape: this.marker === 'circle' ? {
          cx: xPos(0),
          cy: yPos(0),
          r: this.markerSize / 2,
        } : this.marker === 'ring' ? {
          cx: xPos(0),
          cy: yPos(0),
          r: this.markerSize / 2,
          r0: this.markerSize / 4,
        } : {
          x: xPos(this.markerSize / 2),
          y: yPos(this.markerSize / 2),
          width: this.markerSize,
          height: this.markerSize,
        },
        style: {
          fill: 'red',
        /*
        fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
          this.data[params.dataIndex][3][0] : api.value(3)]
          */
        },
        x: 0,
        y: 0,
        rotation: 0,
      };

      if (this.marker === 'diamond') {
        marker.type = 'rect';
        marker.x = xPos(0);
        marker.y = yPos(0);
        marker.shape.x = -this.markerSize / 2;
        marker.shape.y = -this.markerSize / 2;
        marker.shape.r = this.markerSize / 4;
        marker.rotation = 0.785398;
      } else if (this.marker === 'rect') {
        marker.x = 0;
        marker.y = 0;
        marker.shape.x = xPos(this.markerSize / 2);
        marker.shape.y = yPos(this.markerSize / 2);
        marker.shape.r = 0;
        marker.rotation = 0;
      }

      group.children.push(marker);
    } else {
    // draw line
      if (!this.echarts)
        debugger;

      const rectShape = this.echarts.graphic.clipRectByRect({
        x: start[0],
        y: start[1] - this.lineWidth / 2,
        width: width,
        height: this.lineWidth,
        //  type: 'rect'
      }, {
        x: params.coordSys.x,
        y: params.coordSys.y,
        width: params.coordSys.width,
        height: params.coordSys.height,
      });

      if (overlap) {
        const height = api.size([0, 1])[1];
        const offset = Math.max(this.markerSize * 2, this.lineWidth);
        // Shift along the Y axis
        //@ts-ignore
        rectShape.y += (this.count % 3) ? (this.count % 3 === 2) ?
          0 : offset-height/2 : height/2-offset; //@ts-ignore
        this.count += 1;
      }

      group.children.push({
        type: 'rect',
        transition: ['shape'],
        shape: rectShape,
        // style: { fill: api.value(4) ? this.selectionColor : this.colorMap[isNaN(api.value(3)) ?
        // this.data[params.dataIndex][3][0] : api.value(3)] }
        style: {fill: 'green'},
      });
    }

    return group;
  }

  initTimeLine(visibleIndex: number, viewerThis: any) {
    const r = {
      type: 'custom',
      progressive: 0,
      animation: false,
      renderItem: this.timeLinesD.bind(viewerThis),
      encode: {
        x: [startCol, endCol],
        y: categoryCol,
        tooltip: 3,
      },
      data: this.timeLinesData[visibleIndex],
    };
    return r;
  }
}
