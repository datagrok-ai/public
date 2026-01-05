/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PhylocanvasGlViewer} from './viewers/phylocanvas-gl-viewer';
import {PhylocanvasGlViewerApp} from './apps/phylocanvas-gl-viewer-app';
import {TreeToGridApp} from './apps/tree-to-grid-app';
import {injectTreeToGridUI} from './viewers/inject-tree-to-grid';
import {TreeInGridCellApp} from './apps/tree-in-grid-cell-app';
import {PhylocanvasGlService} from './utils/phylocanvas-gl-service';
import {TreeCutAsTreeApp} from './apps/tree-cut-as-tree-app';
import {findNewick} from './scripts-api';
import {PhylocanvasGlServiceBase} from '@datagrok-libraries/bio/src/viewers/phylocanvas-gl-viewer';
import {getTreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';

export * from './package.g';
export const _package = new DG.Package();

type PtvWindowType = Window & { $phylocanvasGlService?: PhylocanvasGlService };
declare const window: PtvWindowType;

export class PackageFunctions {
  @grok.decorators.func({
    meta: {icon: 'files/icons/phylocanvasgl-viewer.svg', role: 'viewer'},
    name: 'PhylocanvasGL',
    description: 'Phylogenetic tree visualization',
    outputs: [{type: 'viewer', name: 'result'}]
  })
  static phylocanvasGlViewer(): PhylocanvasGlViewer { 
    return new PhylocanvasGlViewer();
  }



  @grok.decorators.func({
    description: 'Test/demo app for PhylocanvasGlViewer'
  })
  static async phylocanvasGlViewerApp(): Promise<void> {
  
    const pi = DG.TaskBarProgressIndicator.create('open PhylocanvasGlViewer app');
    try {
      const app = new PhylocanvasGlViewerApp('phylocanvasGlViewerApp');
      await app.init();
    } catch (err: unknown) {
      const msg: string = 'PhyloTreeViewer phylocanvasGlViewerApp() error: ' +
        `${err instanceof Error ? err.message : (err as Object).toString()}`;
      grok.shell.error(msg);
      console.error(msg);
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    name: 'TreeToGrid',
    description: 'Test/demo app for TreeToGrid (PhylocanvasGL based)',
    outputs: [{type: 'object', name: 'result'}]
  })
  static async treeToGridApp(): Promise<TreeToGridApp> {
  
    const pi = DG.TaskBarProgressIndicator.create('open treeInGrid app');
    try {
      const app = new TreeToGridApp();
      await app.init();
      return app;
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({
    name: 'TreeCutAsTree',
    description: 'Test/demo app for TreeCutAsTree'
  })
  static async treeCutAsTreeApp(): Promise<void> {
  
    const pi = DG.TaskBarProgressIndicator.create('open treeCutAsTree app');
    try {
      const app = new TreeCutAsTreeApp();
      await app.init();
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    name: 'TreeInGridCell',
    description: 'Test/demo app for TreeInGridCell'
  })
  static async treeInGridCellApp(): Promise<void> {
  
    const pi = DG.TaskBarProgressIndicator.create('open TreeInGridCell app');
    try {
      const app = new TreeInGridCellApp();
      await app.init();
    } catch (err: unknown) {
      const msg: string = 'PhyloTreeViewer treeInGridCellApp() error: ' +
        `${err instanceof Error ? err.message : (err as Object).toString()}`;
      grok.shell.error(msg);
      //@ts-ignore
      console.error(err);
      // if ('stack' in err)
      //   console.error(err['stack']);
    } finally {
      pi.close();
    }
  }


  @grok.decorators.func({
    name: 'injectTree',
    description: 'Opens Newick file'
  })
  static async injectTreeToGrid(
    @grok.decorators.param({'type':'viewer'})  grid: DG.Grid,
    newickText: string,
    @grok.decorators.param({'options':{'optional':true}})   leafColName?: string) {
  
    const colNameList: string[] = grid.dataFrame.columns.names();
    leafColName = leafColName ??
      grid.dataFrame.getTag('.newickLeafColumn') ??
      colNameList.find((colName) => colName.toLowerCase() == 'node') ??
      colNameList.find((colName) => colName.toLowerCase() == 'leaf') ??
      colNameList.find((colName) => colName.toLowerCase() == 'id');
    if (!leafColName)
      throw new Error('The leaf column name can not be inferred. Specify it as an argument.');
    await injectTreeToGridUI(grid, newickText, leafColName!);
  }


  @grok.decorators.func({outputs: [{type: 'object', name: 'result'}]})
  static getPhylocanvasGlService(): PhylocanvasGlServiceBase {
  
    if (!(window.$phylocanvasGlService)) {
      const svc: PhylocanvasGlService = new PhylocanvasGlService();
      window.$phylocanvasGlService = svc;
    }

    return window.$phylocanvasGlService;
  }
}