/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {fit, sigmoid, FitErrorModel} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';
export const _package = new DG.Package();


//tags: app
//name: Demo Curve Fit
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
      const data = {x: x, y: y};
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
      const fitRes = fit(data, [maxY, 1.2, xAtMedY, minY], sigmoid, FitErrorModel.Constant);
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