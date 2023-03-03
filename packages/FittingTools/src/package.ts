/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {fit, sigmoid, FitErrorModel} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';
import {PeService} from "./parameter-estimation-service/pe-service";
import {getConfidence} from "@datagrok-libraries/statistics/src/confidence-intervals";
export const _package = new DG.Package();


//tags: app
//name: Neo Fit
export async function NeoFit() {
  DG.Utils.openFile({
    accept: '.csv',
    open: async(file) => {
      const str = await file.text();
      const dframe = grok.data.parseCsv(str);
      const view = grok.shell.addTableView(dframe);
      const grid = view.grid;
      processNeoFit(grid);
    }
  });
}

export function processNeoFit(grid: DG.Grid) : void {
  const dframe = grid.dataFrame;
  const colPlot = dframe.columns.byName('PLOT');
  if (colPlot === null)
    return;
  let validCount = 0;
  const columnFit = DG.Column.fromType(DG.TYPE.OBJECT, 'PLOT_DF', dframe.rowCount).init((n) => {
    const str = colPlot.get(n);
    let json = null;
    try {
      json = JSON.parse(str);
    } catch (e) {
      return;
    }

    ++validCount;
    const xAxisName = json.xAxisName + ', ' +  json.xAxisUnit;
    const yAxisName = json.yAxisName;

    const arCurves: Array<object> = json.curves;
    if( arCurves.length === 0)
      return;

    let pointCount = 0;
    for (let j = 0; j < arCurves.length; ++j) {
      const curve: any = arCurves[j];
      const points: Array<object> = curve.points;
      pointCount += points.length;
    }

    const arFLs = new Array<any>(arCurves.length);
    const arCfs = new Array<any>();
    const arX = new Array<number>(pointCount);
    const arY = new Array<number>(pointCount);
    const arCurveId = new Array<string>(pointCount);
    let yMaxOverall = -2000.0;
    pointCount = 0;
    for (let j = 0; j < arCurves.length; ++j) {
      const curve: any = arCurves[j];
      const points: Array<object> = curve.points;

      const arXFit = new Array<number>(points.length);
      const arYFit = new Array<number>(points.length);
      for (let k = 0; k < points.length; ++k) {
        arCurveId[k + pointCount] = 'a'+j.toString();
        arX[k + pointCount] = parseFloat((points[k] as any).X.toString());
        arX[k + pointCount] = Math.log10(arX[k+ pointCount]);
        arXFit[k] = arX[k + pointCount];

        arY[k + pointCount] = parseFloat((points[k] as any).Y.toString());
        arYFit[k] = arY[k + pointCount];
      }

      pointCount += points.length;
      const segYCopy = arYFit.slice(0);
      segYCopy.sort((a, b) => a - b);
      //console.log(segYCopy);
      const minY = segYCopy[0];
      const maxY = segYCopy[segYCopy.length -1];
      if (maxY > yMaxOverall)
        yMaxOverall = maxY;

      const idxMed = Math.ceil(segYCopy.length / 2);
      const medY = segYCopy[idxMed];
      let xAtMedY = -1;
      for (let k = 0; k < arY.length; ++k) {
        if (arYFit[k] === medY) {
          xAtMedY = arXFit[k];
          break;
        }
      }

      function binByX(x: number[], y: number[]) : Map<number, number[]> {
        const map = new Map();
        let ar = null;
        for (let n = 0; n < x.length; ++n) {
          ar = map.get(x[n]);
          if (ar === undefined) {
            ar = new Array<number>();
            map.set(x[n], ar);
          }
          ar.push(y[n]);
        }
        return map;
      }

      const cr = DG.Color.toHtml(DG.Color.getCategoricalColor(j));
      const fitRes = fit({x: arXFit, y: arYFit}, [maxY, 1.2, xAtMedY, minY], sigmoid, FitErrorModel.Constant);
      const params = fitRes.parameters;
      const iceptY = fitRes.fittedCurve(params[2]);
      const map = binByX(arXFit, arYFit);
      const xKeys = Array.from(map.keys());
      let conf = null;
      let ar = null;

      arCfs.push({
        title: ' ',
        description: ' ',
        formula: '${' + xAxisName + '} = ' + params[2].toString(),
        min: 0,
        max: iceptY,
        zIndex: -29,
        color: cr,
        opacity: 100,
        width: 1,
        style: 'dotted',
        visible: true,
      });

      arCfs.push({
        title: ' ',
        description: ' ',
        formula: '${' + yAxisName + '} = ' + iceptY.toString(),
        max: params[2],
        zIndex: -29,
        color: cr,
        opacity: 100,
        width: 1,
        style: 'dotted',
        visible: true,
      });

      for (let m = 0; m < xKeys.length; ++m) {
        ar = map.get(xKeys[m]);
        if (ar === undefined)
          continue;

        conf = getConfidence(ar);
        arCfs.push({
          title: ' ',
          description: ' ',
          formula: '${' + xAxisName + '} = ' + xKeys[m] + (0.001 * j).toString(),
          min: conf.bottom,
          max: conf.top,
          zIndex: -29,
          color: cr,
          opacity: 40,
          width: 10,
          visible: true,
        });
      }

      arFLs[j] = {
        title: ' ',
        description: ' ',
        formula: '${' + yAxisName + '} = ' + params[3] + ' + (' + params[0] + ' - ' + params[3] + ')/(1 + Pow(10, ' + params[1] + '*(${' + xAxisName + '} - ' + params[2] + ')))',
        zIndex: -30,
        color: cr,
        width: 2,
        visible: true,
      };
    }

    const dfOut = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array(xAxisName, new Float32Array(arX)),
      DG.Column.fromFloat32Array(yAxisName, new Float32Array(arY)),
      DG.Column.fromStrings('CurveId', arCurveId)]);

    const helper = dfOut.plot.scatter({x: xAxisName, y: yAxisName, markerDefaultSize: 4, markersColumnName: 'CurveId', colorColumnName: 'CurveId'});
    helper.setOptions({yMin: 0, yMax: Math.max(100, yMaxOverall)});

    for (let i = 0; i < arCurves.length; ++i) {
      helper!.meta.formulaLines.addLine(arFLs[i]);
    }

    for (let i = 0; i < arCfs.length; ++i) {
      helper!.meta.formulaLines.addLine(arCfs[i]);
    }
    return helper;
  });
  dframe.columns.add(columnFit);
  dframe.columns.remove('PLOT');

  grid.setOptions({rowHeight: 150});

  const columnGridFit = grid.columns.byName('PLOT_DF');
  columnGridFit!.width = 200;
  columnGridFit!.cellType = 'scatterplot';
}


//tags: app
//name: Demo PE Workers
export async function runPEWorkers() {
  const x : number[] = [];
  const y : number[] = [];
  const columnCount = 50000;
  const counts : number[] = new Array<number>(columnCount);

  for (let nCol=0; nCol<columnCount; ++nCol) {
    const dff = grok.data.demo.doseResponse();
    const xRand = dff.columns.byName('concentration [ug/ml]').toList();
    const yRand = dff.columns.byName('viability [%]').toList();
    x.push(...xRand);
    y.push(...yRand);
    counts[nCol] = xRand.length;
  }

 const service = await PeService.create();
 const timeStart = new Date().getTime();
 const params = await service.fit(x, y, counts);
 const timeEnd = new Date().getTime();
 console.log('Spent: ' + Math.floor((timeEnd - timeStart)));
 //console.log(params);
}


//tags: app
//name: Demo Curve Fitting Service
export async function curveFitRendererAsync() {
  const X: number[] = [];
  const Y: number[] = [];

  const columnCount = 10;
  const rowCount = 10;
  const totalFitCount = columnCount * rowCount;
  const counts: number[] = new Array<number>(totalFitCount);

  for (let cellIdx = 0; cellIdx < totalFitCount; ++cellIdx) {
    const dff = grok.data.demo.doseResponse();
    const xRand = dff.columns.byName('concentration [ug/ml]').toList();
    const yRand = dff.columns.byName('viability [%]').toList();
    X.push(...xRand);
    Y.push(...yRand);
    counts[cellIdx] = xRand.length;
  }

  const service = await PeService.create();
  const timeStart = new Date().getTime();
  const allParams = await service.fit(X, Y, counts);
  const timeEnd = new Date().getTime();
  const paramCount = 4;
  console.log('Check output: ' + (allParams.length / paramCount === totalFitCount));

  let shift = 0;
  let fitIdx = 0;
  const strIC50ColName = 'FittedCurve';
  const columnsFit = new Array(columnCount);
  for (let nCol=0; nCol<columnCount; ++nCol) {
    const columnFit = DG.Column.fromType(DG.TYPE.DATA_FRAME, strIC50ColName + nCol.toString(), rowCount).init((i) => {
      const params = allParams.slice(paramCount*fitIdx, paramCount*(fitIdx +1));
      const x = X.slice(shift, shift + counts[fitIdx]);
      const y = Y.slice(shift, shift + counts[fitIdx]);
      shift += counts[fitIdx];
      ++fitIdx;
      const dfOut = DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array('X', new Float32Array(x)),
        DG.Column.fromFloat32Array('Y', new Float32Array(y))]);
      const scatterOptions = {
        x: 'X',
        y: 'Y',
        markerDefaultSize: 4
      };

      const sp = dfOut.plot.scatter(scatterOptions);
      sp!.meta.formulaLines.addLine({
        title: 'Fitted',
        formula: '${Y} = ' + params[3] + ' + (' + params[0] + ' - ' + params[3] + ')/(1 + Pow(10, ' + params[1] + '*(${X} - ' + params[2] + ')))',
        zIndex: -30,
        color: "#FF7F00",
        width: 2,
        visible: true,
      });

      return sp;
    });

    columnsFit[nCol] = columnFit;
  }

  const t = DG.DataFrame.fromColumns(columnsFit);
  const view: DG.TableView = grok.shell.addTableView(t);
  view.grid.setOptions({rowHeight: 150});

  for (let n=1; n<view.grid.columns.length; ++n) {
    view.grid.columns.byIndex(n)!.width = 200;
    view.grid.columns.byIndex(n)!.cellType = 'scatterplot';
  }
}


//tags: app
//name: Demo Curve Fittings
export async function curveFitRenderer() {
  const strIC50ColName = 'FittedCurve';
  const columnCount = 10;
  const rowCount = 100;
  const columnsFit = new Array(columnCount);
  for (let nCol=0; nCol<columnCount; ++nCol) {
    const columnFit = DG.Column.fromType(DG.TYPE.DATA_FRAME, strIC50ColName + nCol.toString(), rowCount).init((i) => {
      const dff = grok.data.demo.doseResponse();
      const x = dff.columns.byName('concentration [ug/ml]').toList();
      const y = dff.columns.byName('viability [%]').toList();

      //fit
      const minY = dff.columns.byName('viability [%]').min;
      const maxY = dff.columns.byName('viability [%]').max;
      const medY = dff.columns.byName('viability [%]').stats.med;
      let xAtMedY = -1;
      for (let n = 0; n < y.length; ++n) {
        if (y[n] === medY) {
          xAtMedY = x[n];
          break;
        }
      }

      //Example of fit curve usage
      const fitRes = fit({x: x, y: y}, [maxY, 1.2, xAtMedY, minY], sigmoid, FitErrorModel.Constant);
      const params = fitRes.parameters;

      const dfOut = DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array('X', new Float32Array(x)),
        DG.Column.fromFloat32Array('Y', new Float32Array(y))]);
      const scatterOptions = {
        x: 'X',
        y: 'Y',
        markerDefaultSize: 4
      };

      const helper = dfOut.plot.scatter(scatterOptions);
      helper!.meta.formulaLines.addLine({
        title: 'Fitted',
        formula: '${Y} = ' + params[3] + ' + (' + params[0] + ' - ' + params[3] + ')/(1 + Pow(10, ' + params[1] + '*(${X} - ' + params[2] + ')))',
        zIndex: -30,
        color: "#FF7F00",
        width: 2,
        visible: true,
      });

      //dfOut.temp = helper;
      return helper;
    });

    columnsFit[nCol] = columnFit;
  }

  const t = DG.DataFrame.fromColumns(columnsFit);
  const view: DG.TableView = grok.shell.addTableView(t);
  view.grid.setOptions({rowHeight: 150});

  for (let n=1; n<view.grid.columns.length; ++n) {
    view.grid.columns.byIndex(n)!.width = 200;
    view.grid.columns.byIndex(n)!.cellType = 'scatterplot';
  }
}

//name: Curves
//tags: app
export async function sim() {
  const dataX = new Float32Array([4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8]);
  const dataY = new Float32Array([1.005, 0.981, 0.921, 0.987,
                                  0.948, 0.810, 0.730, 0.558,
                                  0.501, 0.381, 0.271, 0.148,
                                  0.139, 0.078, 0.095, -0.037, 0.01]);

  let t = DG.DataFrame.fromColumns([
    DG.Column.fromFloat32Array('X', dataX),
    DG.Column.fromFloat32Array('Y', dataY)
  ]);

  let lc = DG.Viewer.scatterPlot(t);

  const x = [4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8];
  const y = [1.005, 0.981, 0.921, 0.987, 0.948, 0.810, 0.730, 0.558, 0.501, 0.381, 0.271, 0.148, 0.139, 0.078, 0.095, -0.037, 0.01];
  const data = {x: x, y: y};

  //Example of fit curve usage
  let fitRes = fit(data, [0.9, 1.2, 4.95, 0.1], sigmoid, FitErrorModel.Constant);
  //const curveFitter = new CurveFitter(sigmoid, data, [0.9, 1.2, 4.95, 0.1], "constant");

  let btn = ui.bigButton('FIT', () => {

    let params = fitRes.parameters;
    let fittedFunction = fitRes.fittedCurve;
    let confidenceTop = fitRes.confidenceTop;
    let confidenceBottom = fitRes.confidenceBottom;

    for (let i = 0; i < x.length; ++i) {
      console.log("y = " + y[i]);
      console.log(fittedFunction(x[i]));
      console.log(confidenceTop(x[i]));
      console.log(confidenceBottom(x[i]));
    }

    //Example of fit curve usage ends



    lc.meta.formulaLines.removeAt(-30);
    lc.meta.formulaLines.addLine({
      title: 'Fitted',
      formula: '${Y} = ' + params[3] +' + (' + params[0] +' - ' + params[3] + ')/(1 + Pow(10, ' + params[1] + '*(${X} - ' + params[2] +')))',
      zIndex: -30,
      color: "#FF7F00",
      width: 2,
      visible: true,
    });

    // const top =  1.4*out + params[3];
    // const bottom =  -1.4*out + params[3];

    // lc.meta.formulaLines.addLine({
    //   title: 'Fitted',
    //   formula: '${Y} = ' + top +' + (' + params[0] +' - ' + params[3] + ')/(1 + Pow(10, ' + params[1] + '*(${X} - ' + params[2] +')))',
    //   zIndex: -30,
    //   color: "#1e5945",
    //   width: 2,
    //   visible: true,
    // });

    // lc.meta.formulaLines.addLine({
    //   title: 'Fitted',
    //   formula: '${Y} = ' + bottom +' + (' + params[0] +' - ' + params[3] + ')/(1 + Pow(10, ' + params[1] + '*(${X} - ' + params[2] +')))',
    //   zIndex: -30,
    //   color: "#1e5945",
    //   width: 2,
    //   visible: true,
    // });
  });


  lc.meta.formulaLines.addLine({
    title: 'Original',
    formula: '${Y} = 0.1 + (0.9 - 0.3)/(0.7 + Pow(10, 1.2*(${X} - 4.95)))',
    zIndex: -30,
    color: "#FFA500",
    width: 2,
    visible: true,
  });

  const v = grok.shell.newView('Curve fit', [
    ui.panel([
      ui.div([
        ui.div([
          ui.divH([ui.h1('Inputs')]),
          ui.divV([
            btn
          ], 'ui-form'),
        ], 'ui-form'),
      ], 'ui-form'),
    ])
  ]);
  v.box = true;

  v.append(lc);
}
