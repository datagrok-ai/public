/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {testEvent} from '@datagrok-libraries/test/src/test';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {App, HweWindow} from '@datagrok-libraries/bio/src/helm/types';
import {HelmInputBase, IHelmHelper, IHelmInputInitOptions} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';

import {getPropertiesWidget} from './widgets/properties-widget';
import {HelmGridCellRenderer, HelmGridCellRendererBack} from './utils/helm-grid-cell-renderer';
import {_getHelmService, HelmPackage} from './package-utils';
import {HelmHelper} from './helm-helper';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import type {JSDraw2Module, OrgHelmModule, ScilModule} from './types';

export * from './package.g';
export const _package = new HelmPackage(/*{debug: true}/**/);

/** Temporary polyfill */

function getDecoratorFunc() {
  return function(args: any) {
    return function(
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor
    ) { };
  };
}

// Ensure decorators object exists and polyfill missing decorators
if (!grok.decorators)
  (grok as any).decorators = {};

const decorators = [
  'func', 'init', 'param', 'panel', 'editor', 'demo', 'app',
  'appTreeBrowser', 'fileHandler', 'fileExporter', 'model', 'viewer', 'filter', 'cellRenderer', 'autostart',
  'dashboard', 'folderViewer', 'semTypeDetector', 'packageSettingsEditor', 'functionAnalysis', 'converter',
  'fileViewer', 'model', 'treeBrowser', 'polyfill'
];

decorators.forEach((decorator) => {
  if (!(grok.decorators as any)[decorator])
    (grok.decorators as any)[decorator] = getDecoratorFunc();
});

/** End temporary polyfill */


/*
  Loading modules:
  Through Helm/package.json/sources section
    dojo from ajax.googleapis.com
    HELMWebEditor
      JSDraw.Lite is embedded into HELMWebEditor bundle (dist/package.js)
 */

declare const window: Window & HweWindow;
declare const scil: ScilModule;
declare const JSDraw2: JSDraw2Module;
declare const org: OrgHelmModule;

let initHelmPromise: Promise<void> | null = null;

async function initHelmInt(): Promise<void> {
  const logPrefix: string = 'Helm: _package.initHelm()';
  _package.logger.debug(`${logPrefix}, start`);

  try {
    await getRdKitModule();
    const seqHelper = await getSeqHelper();
    const libHelper = await getMonomerLibHelper();
    await _package.initHELMWebEditor();
    const helmHelper: IHelmHelper = new HelmHelper(seqHelper, _package.logger);
    _package.completeInit(helmHelper, libHelper);
  } catch (err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.error(`Package \'Helm\' init error:\n${errMsg}`);
    const errRes = new Error(`${logPrefix} error:\n  ${errMsg}\n${errStack}`);
    errRes.stack = errStack;
    throw errRes;
  } finally {
    _package.logger.debug(`${logPrefix}, finally`);
  }
  _package.logger.debug(`${logPrefix}, end`);
}

function openWebEditor(cell: DG.GridCell, value?: string, units?: string) {
  const view = ui.div();
  const col = cell.cell.column as DG.Column<string>;
  const sh = _package.seqHelper.getSeqHandler(col);
  const app: App =
    _package.helmHelper.createWebEditorApp(view, !!cell && units === undefined ? cell.cell.value : value!);
  const dlg = ui.dialog({showHeader: false, showFooter: true});
  dlg.add(view)
    .onOK(() => {
      const helmValue: string = app.canvas!.getHelm(true).replace(/<\/span>/g, '')
        .replace(/<span style='background:#bbf;'>/g, '');
      if (!!cell) {
        if (units === undefined)
          cell.setValue(helmValue);
        else {
          const convertedRes = sh.convertHelmToFastaSeparator(helmValue, units!, sh.separator);
          cell.setValue(convertedRes);
        }
      }
    }).show({modal: true, fullScreen: true});

  const dlgCntDiv = $(dlg.root).find('div').get()[0] as HTMLDivElement;
  dlgCntDiv.className = dlgCntDiv.className.replace('dlg- ui-form', 'dlg-ui-form');
}

function checkMonomersAndOpenWebEditor(cell: DG.GridCell, value?: string, units?: string) {
  openWebEditor(cell, value, units);
}

export class PackageFunctions {
  @grok.decorators.init()
  static async initHelm(): Promise<void> {
    if (initHelmPromise === null)
      initHelmPromise = initHelmInt();
    return initHelmPromise;
  }


  @grok.decorators.func({
    'description': 'Helm renderer service',
    'outputs': [{'type': 'object', 'name': 'result'}]
  })
  static async getHelmService(): Promise<HelmServiceBase> {
    return _getHelmService();
  }

  @grok.decorators.func({
    'meta': {
      'columnTags': 'quality=Macromolecule, units=helm',
      'cellType': 'helm',
      'role': 'cellRenderer'
    },
    'outputs': [{name: 'result', type: 'grid_cell_renderer'}]
  })
  static helmCellRenderer(): DG.GridCellRenderer {
    const logPrefix = `Helm: _package.getHelmCellRenderer()`;
    _package.logger.debug(`${logPrefix}, start`);
    // return new HelmCellRenderer(); // old
    return new HelmGridCellRenderer(); // new
  }


  @grok.decorators.func({
    'meta': {
      'columnTags': 'quality=Macromolecule, units=helm',
      'role': 'cellEditor'
    },
    'description': 'Macromolecule'
  })
  static editMoleculeCell(
    @grok.decorators.param({'type': 'grid_cell'}) cell: DG.GridCell): void {
    checkMonomersAndOpenWebEditor(cell, undefined, undefined);
  }


  @grok.decorators.func({
    'meta': {'action': 'Edit Helm...'},
    'name': 'Edit Helm...',
    'description': 'Adds editor'
  })
  static openEditor(
    @grok.decorators.param({'options': {'semType': 'Macromolecule'}}) mol: DG.SemanticValue): void {
    const df = grok.shell.tv.grid.dataFrame;
    const col = df.columns.bySemType('Macromolecule')! as DG.Column<string>;
    const colSh = _package.seqHelper.getSeqHandler(col);
    const colUnits = col.meta.units;
    if (df.currentRowIdx === -1)
      return;
    const gCell = DG.GridCell.fromColumnRow(grok.shell.tv.grid, col.name, df.currentRowIdx);
    if (colUnits === NOTATION.HELM)
      checkMonomersAndOpenWebEditor(gCell, undefined, undefined);
    else {
      const convert = colSh.getConverter(NOTATION.HELM);
      const helmMol = convert(mol.value);
      checkMonomersAndOpenWebEditor(gCell, helmMol, col.meta.units!);
    }
  }


  @grok.decorators.panel({
    'name': 'Properties',
    'meta': {role: 'widgets', domain: 'bio'},
  })
  static propertiesWidget(
    @grok.decorators.param({'options': {'semType': 'Macromolecule'}}) sequence: DG.SemanticValue): DG.Widget {
    return getPropertiesWidget(sequence);
  }

  @grok.decorators.func()
  static getMolfiles(
    @grok.decorators.param({'type': 'column', 'options': {'semType': 'Macromolecule'}}) col: DG.Column<string>): DG.Column<string> {
    const helmStrList = col.toList();
    const molfileList = _package.helmHelper.getMolfiles(helmStrList);
    const molfileCol = DG.Column.fromStrings('mols', molfileList);
    return molfileCol;
  }


  @grok.decorators.func({
    'meta': {
      'propertyType': 'string',
      'semType': 'Macromolecule',
      'role': 'valueEditor'
    },
    'outputs': [{'type': 'object', 'name': 'result'}]
  })
  static helmInput(
    @grok.decorators.param({'options': {'optional': true}}) name: string,
    @grok.decorators.param({'type': 'object', 'options': {'optional': true}}) options: IHelmInputInitOptions): HelmInputBase {
    // TODO: Annotate for semType = 'Macromolecule' AND units = 'helm'
    return _package.helmHelper.createHelmInput(name, options);
  }


  @grok.decorators.func({'outputs': [{'type': 'object', 'name': 'result'}]})
  static getHelmHelper(): IHelmHelper {
    return _package.helmHelper;
  }


  @grok.decorators.func()
  static async measureCellRenderer(): Promise<void> {
    const grid = grok.shell.tv.grid;
    const gridCol = grid.columns.byName('sequence')!;
    const back = gridCol.temp.rendererBack as HelmGridCellRendererBack;

    let etSum: number = 0;
    let etCount: number = 0;
    for (let i = 0; i < 20; ++i) {
      const t1 = window.performance.now();
      let t2: number;
      if (!back.cacheEnabled) {
        await testEvent(back.onRendered, () => {
          t2 = window.performance.now();
          _package.logger.info(`measureCellRenderer() cache disabled , ET: ${t2 - t1} ms`);
        }, () => {
          back.invalidate(); // grid.invalidate();
        }, 5000);
      } else {
        await testEvent(grid.onAfterDrawContent, () => {
          t2 = window.performance.now();
          _package.logger.info(`measureCellRenderer() cache enabled, ET: ${t2! - t1} ms`);
        }, () => {
          grid.invalidate();
        }, 5000);
      }

      etSum += (t2! - t1);
      etCount++;
    }
    _package.logger.info(`measureCellRenderer(), avg ET: ${etSum / etCount} ms`);
  }


  @grok.decorators.func({
    'name': 'Highlight Monomers'
  })
  static async highlightMonomers(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Test app highlight monomers');
    try {
      const df = await _package.files.readCsv('tests/peptide.csv');
      df.name = 'peptide.csv';
      await initHelmPromise;

      const seqCol = df.columns.byName('helm')!;
      await grok.data.detectSemanticTypes(df);
      await grok.functions.call('Bio:toAtomicLevel', {table: df, seqCol: seqCol, nonlinear: true, highlight: true});

      // const resView = DG.TableView.create(df, false);
      // return resView;
      grok.shell.addTableView(df);
    } finally {
      pi.close();
    }
  }
}
