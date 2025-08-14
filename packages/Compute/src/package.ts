/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OutliersSelectionViewer} from './outliers-selection/outliers-selection-viewer';
import {
  CustomFunctionView as CustomFunctionViewInst,
} from '@datagrok-libraries/compute-utils';
import {ModelCatalogView, ModelHandler, startModelCatalog, makeModelTreeBrowser, renderRestPanel, setModelCatalogEventHandlers, setModelCatalogHandler} from '@datagrok-libraries/compute-utils/model-catalog';

import {FittingView} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';

export const _package = new DG.Package();

//name: openModelFromFuncall
//input: funccall funccall
export function openModelFromFuncall(funccall: DG.FuncCall) {
  ModelHandler.openModelFromFunccall(funccall);
}

//name: OutliersSelectionViewer
//description: Creates an outliers selection viewer
//tags: viewer
//output: viewer result
export function OutliersSelection() {
  return new OutliersSelectionViewer();
}

//name: renderRestPanel
//input: func func
//output: widget panel
export async function renderPanel(func: DG.Func): Promise<DG.Widget> {
  return renderRestPanel(func);
}

let startUriLoaded = false;
let initCompleted = false;

const options = {
  _package,
  ViewClass: ModelCatalogView,
  segment: 'Modelhub',
  viewName: 'Model Hub',
  funcName: 'modelCatalog',
  setStartUriLoaded: () => startUriLoaded = true,
  getStartUriLoaded: () => startUriLoaded,
};

//tags: init
export function init() {
  if (initCompleted)
    return;

  setModelCatalogHandler();
  setModelCatalogEventHandlers(options);

  initCompleted = true;
}

//name: Model Hub
//tags: app
//output: view v
//meta.browsePath: Compute
export function modelCatalog() {
  return startModelCatalog(options);
}

//input: dynamic treeNode
//input: view browseView
export function modelCatalogTreeBrowser(treeNode: DG.TreeViewGroup) {
  makeModelTreeBrowser(treeNode);
}

////
// Compute-utils API section
///

export const CFV = CustomFunctionViewInst;

//name: fitTestFunc
//description: Test for optimization: multiple scalars output
//input: double x1 = 1 {caption: param1; min: -3; max: 3}
//input: double x2 = -1 {caption: param2; min: -3; max: 3}
//input: dataframe y {caption: table}
//input: bool bool
//output: int integer
//output: double float1
//output: double float2
//output: dataframe table1 {viewer: Line chart(block: 60) | Grid(block: 40) }
//output: dataframe table2 {viewer: Line chart(block: 60) | Grid(block: 40) }
//editor: Compute2:RichFunctionViewEditor
//meta.features: {"fitting": true, "sens-analysis": true}
//meta.runOnInput: true
export function fitTestFunc(x1: number, x2: number, y: DG.DataFrame, bool: boolean) {
  return {
    integer: x1**3 * (x1 - 1) * x2**3 * (x2 - 1),
    float1: (x2 - 1)**2 + (x1 - 1)**2,
    float2: (x2 - 1)**4 + (x1 - 1)**4,
    table1: DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [x1 + x2 + 1, 4, 9, 16, 25]),
    ]),
    table2: DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [x1 + x2 + 1, 8, 27, 64, 125]),
    ]),
  };
}

//name: testFittingOutputs
//description: Test for optimization: multiple scalars output
export async function testFittingOutputs() {
  const func = await grok.functions.find('Compute:fitTestFunc');

  if (func === null) {
    grok.shell.error('The function "Compute:fitTestFunc" not found!');
    return;
  }

  const targetDf1 = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [1, 4.3, 9.1, 16, 25]),
  ]);
  targetDf1.name = 'test-df1';

  const targetDf2 = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'arg', [1, 2, 3, 4, 5]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'func', [1, 8.5, 27.6, 64.9, 125]),
  ]);
  targetDf2.name = 'test-df2';

  await FittingView.fromEmpty(func, {
    targets: {
      integer: {default: 123, enabled: true},
      float1: {default: 456.789, enabled: true},
      table1: {
        default: targetDf1,
        enabled: true,
        argumentCol: 'arg',
      },
      table2: {
        default: targetDf2,
        enabled: true,
        argumentCol: 'arg',
      },
    },
  });
}
