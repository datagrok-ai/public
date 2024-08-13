import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, before, after, awaitCheck, delay, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {createTableView, readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ScaffoldTreeViewer} from '../widgets/scaffold-tree';

category('scaffold tree', () => {
  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  // check that scaffold viewer openes without errors
  test('scaffoldTreeViewerOpens', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('Demo:Files/chem/smiles_50K.zip') :
      await readDataframe('tests/sar-small_test.csv');
    const tv = grok.shell.addTableView(df);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    tv.addViewer(ScaffoldTreeViewer.TYPE);
    await awaitCheck(() => Array.from(tv.viewers).filter((it) => it.type === ScaffoldTreeViewer.TYPE).length > 0,
      'cannot create viewer', 3000);
    const generateLink = document.querySelector('.chem-scaffold-tree-generate-hint') as HTMLElement;
    if (generateLink)
      generateLink.click();
    const stviewer = Array.from(tv.viewers).filter((it) => it.type === ScaffoldTreeViewer.TYPE)[0] as ScaffoldTreeViewer;
    await awaitCheck(() => stviewer.root.getElementsByClassName('d4-tree-view-group-host')[0].children.length > 0,
      'scaffold tree has not been generated', DG.Test.isInBenchmark ? 3600000 : 60000);
    await delay(2000); //need to scaffold to finish generation
    tv.close();
  }, {timeout: 70000, benchmark: true, stressTest: true});

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });

  test('parent node contains H atom', async () => {
    const tv = await createTableView('mol1K.csv');
    const scaffoldTree = new ScaffoldTreeViewer();
    const table = tv.dataFrame;
    
    scaffoldTree.dataFrame = table;
    scaffoldTree.molCol = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
    tv.dockManager.dock(scaffoldTree.root);
  
    const molStr = 'Nc1ccccc1';
    const rootGroup = await scaffoldTree.createGroup(molStr, scaffoldTree.tree);
    const childGroup = await scaffoldTree.createGroup(molStr, rootGroup!);
    
    const pencilIcon = childGroup?.root.getElementsByClassName('fa-pencil')[0] as HTMLElement;
    pencilIcon.click();
  
    const dialog = Array.from(document.querySelectorAll('.d4-dialog')).find(dialog => 
      dialog.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Edit Scaffold...'
    );
    
    const label = dialog?.querySelector('.ui-label') as HTMLElement;
    const warningExists = label?.textContent?.trim() === 'The edited molecule is not a superstructure of its parent';
    const isRedLabel = label?.style.color === 'red';
  
    expect(warningExists, false);
    expect(isRedLabel, false);
  });
  
});
