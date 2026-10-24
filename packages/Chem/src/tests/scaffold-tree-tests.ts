import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, before, awaitCheck, delay, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {SCAFFOLD_TREE_MAX_MOLECULES, ScaffoldTreeViewer, scaffoldTreeSampleRows, updateVisibleNodes, value}
  from '../widgets/scaffold-tree';
import {ChemicalEngine, setChemicalEngine} from '../engine/chem-engine';
import {PackageFunctions} from '../package';

category('scaffold tree', () => {
  let mol1K: DG.DataFrame;
  let sarSmall: DG.DataFrame;

  function openScaffoldTree(df: DG.DataFrame):
    {tv: DG.TableView, scaffoldTree: ScaffoldTreeViewer, table: DG.DataFrame} {
    const table = df.clone();
    const tv = grok.shell.addTableView(table);
    const scaffoldTree = new ScaffoldTreeViewer();
    scaffoldTree.dataFrame = table;
    scaffoldTree.molCol = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
    tv.dockManager.dock(scaffoldTree.root);
    return {tv, scaffoldTree, table};
  }

  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    mol1K = await readDataframe('mol1K.csv');
    sarSmall = await readDataframe('tests/sar-small_test.csv');
    await grok.data.detectSemanticTypes(mol1K);
    await grok.data.detectSemanticTypes(sarSmall);
  });

  // check that scaffold viewer openes without errors
  test('scaffoldTreeViewerOpens', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('smiles.csv') :
      await readDataframe('tests/sar-small_test.csv');
    if (!DG.Test.isInBenchmark) //leave only first 5 rows for scaffold tree analysis
      df.rows.removeAt(5, df.rowCount - 5);
    await grok.data.detectSemanticTypes(df);
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
      'scaffold tree has not been generated', DG.Test.isInBenchmark ? 3600000 : 180000);
    await delay(2000); //need to scaffold to finish generation
    tv.close();
  }, {timeout: 190000, benchmark: true, stressTest: true, benchmarkTimeout: 300000, skipReason: 'GROK-17648'});

  test('parent node contains H atom', async () => {
    const {tv, scaffoldTree} = openScaffoldTree(mol1K);
    await delay(1000);

    const molStr = 'Nc1ccccc1';
    const rootGroup = await scaffoldTree.createGroup(molStr, scaffoldTree.tree);
    const childGroup = await scaffoldTree.createGroup(molStr, rootGroup!);

    const pencilIcon = childGroup?.root.getElementsByClassName('fa-pencil')[0] as HTMLElement;
    pencilIcon.click();

    const dialog = Array.from(document.querySelectorAll('.d4-dialog')).find((dialog) =>
      dialog.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Edit Scaffold...',
    );

    const label = dialog?.querySelector('.ui-label') as HTMLElement;
    const warningExists = label?.textContent?.trim() === 'The edited molecule is not a superstructure of its parent';
    const isRedLabel = label?.style.color === 'red';

    expect(warningExists, false);
    expect(isRedLabel, false);
    tv.close();
  });

  test('edit invalid structure', async () => {
    const {tv, scaffoldTree, table} = openScaffoldTree(mol1K);

    // Add valid structure and validate
    const molStr = 'c1ccccc1';
    const rootGroup = await scaffoldTree.createGroup(molStr, scaffoldTree.tree);
    await updateVisibleNodes(scaffoldTree);
    rootGroup?.checkBox?.click();
    expect(table.filter.trueCount, 928);

    // Edit to invalid structure
    const invalidSmiles = 'CC(C)(C)(C)(C)c1ccccc1';
    await editStructure(scaffoldTree, rootGroup!, invalidSmiles);
    await scaffoldTree.updateFilters();
    expect(table.filter.trueCount, 0);

    // Revert to valid structure
    await editStructure(scaffoldTree, rootGroup!, molStr);
    await scaffoldTree.updateFilters();
    expect(table.filter.trueCount, 928);
    tv.close();
  }, {skipReason: 'GROK-16714'});

  test('bitset staleness ignores the valueless root group', async () => {
    const {tv, scaffoldTree} = openScaffoldTree(sarSmall);

    const group = scaffoldTree.createGroup('c1ccccc1', scaffoldTree.tree);
    scaffoldTree.tree.currentItem = scaffoldTree.tree;
    expect(scaffoldTree.tree.currentItem.value, null);
    expect(scaffoldTree.isBitsetStale(scaffoldTree.tree.currentItem), false);
    expect(scaffoldTree.isBitsetStale(group), false);
    tv.close();
  });

  test('opens when the current view is not a table view', async () => {
    grok.shell.addView(DG.View.create());
    const table = sarSmall.clone();
    const viewer = DG.Viewer.fromType(ScaffoldTreeViewer.TYPE, table);
    await awaitCheck(() => viewer.props.moleculeColumnName === 'smiles', 'molecule column is not set', 3000);
    expect(viewer.props.Table, table.name);
    grok.shell.closeAll();
  });

  /** smiles.csv's 1,000 molecules `times` times over, every 100th row empty: a large molecule table. */
  async function largeTable(times: number): Promise<DG.DataFrame> {
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
    const values = Array.from({length: smiles.length * times},
      (_, i) => i % 100 === 99 ? '' : smiles[i % smiles.length]);
    const col = DG.Column.fromStrings('smiles', values);
    col.semType = DG.SEMTYPE.MOLECULE;
    return DG.DataFrame.fromColumns([col]);
  }

  test('a large column is sampled to 750 molecules spread over it', async () => {
    const values = (await largeTable(5)).col('smiles')!.toList() as string[];
    const readable = (m: string) => !!m;
    const rows = scaffoldTreeSampleRows(values, SCAFFOLD_TREE_MAX_MOLECULES, readable);
    expect(rows.length, 750);
    expect(rows.every((r, i) => i === 0 || r > rows[i - 1]), true, 'rows in order, none twice');
    expect(rows.every((r) => values[r] !== ''), true, 'no empty molecule');
    // the first and the last stretch (each values.length / 750 rows) give the first and the last row
    const stretch = Math.ceil(values.length / SCAFFOLD_TREE_MAX_MOLECULES);
    expect(rows[0] < stretch && rows[rows.length - 1] >= values.length - stretch, true,
      `spread from ${rows[0]} to ${rows.at(-1)}`);
    expect(scaffoldTreeSampleRows(values, SCAFFOLD_TREE_MAX_MOLECULES, readable).join(), rows.join(), 'deterministic');
    const small = values.slice(0, 700);
    expect(scaffoldTreeSampleRows(small, SCAFFOLD_TREE_MAX_MOLECULES, readable).length, 700,
      'a small column: all rows');
  });

  test('Generate is open on a large table, from a 750-molecule sample its tooltip names', async () => {
    const {tv, scaffoldTree} = openScaffoldTree(await largeTable(5));
    await delay(500);
    expect(scaffoldTree.generateBlockedReason(), null);
    expect(scaffoldTree.generateHint(), 'Generate from a sample of 750 of 4,950 molecules');
    expect(scaffoldTree._iconGenerate?.classList.contains('inactive'), false, 'Generate is enabled');
    // the tree is generated from the sample (the script is not run here: the generated table is caught)
    let generatedRows = -1;
    scaffoldTree.finishGenerateTree = async (df: DG.DataFrame) => {
      generatedRows = df.rowCount;
    };
    await scaffoldTree.generateTree();
    expect(generatedRows, 750);
    expect(scaffoldTree.generatedFrom, 'a sample of 750 of 4,950 molecules');
    tv.close();
    // a table of at most 750 molecules: as before, all of them
    const small = openScaffoldTree(sarSmall);
    await delay(500);
    expect(small.scaffoldTree.generateBlockedReason(), null);
    expect(small.scaffoldTree.generateHint(), 'Generate from molecular column');
    small.tv.close();
  });

  test('the Python script on a 750-molecule sample', async () => {
    // timing (the stand's Python): the tree the RDKit engine builds from the sample of a 4,950-molecule table. A
    // measurement: with --benchmark only (chemical engine's "scaffold tree: Crux against the Python script" checks
    // the script's tree node by node)
    if (!DG.Test.isInBenchmark)
      return;
    const values = (await largeTable(5)).col('smiles')!.toList() as string[];
    const rows = scaffoldTreeSampleRows(values, SCAFFOLD_TREE_MAX_MOLECULES, (m) => !!m);
    const sample = DG.Column.fromStrings('smiles', rows.map((r) => values[r]));
    sample.semType = DG.SEMTYPE.MOLECULE;
    const df = DG.DataFrame.fromColumns([sample]);
    setChemicalEngine(ChemicalEngine.RDKit);
    try {
      const started = performance.now();
      const tree = JSON.parse(await PackageFunctions.getScaffoldTree(df, 10, false));
      const seconds = (performance.now() - started) / 1000;
      console.log(`scaffold tree | the Python script on ${df.rowCount} molecules: ${seconds.toFixed(1)} s, ` +
        `${tree.length} root scaffolds`);
      expect(tree.length > 0, true);
    } finally {
      setChemicalEngine(null);
    }
  }, {timeout: 900000, benchmark: true});

  async function editStructure(scaffoldTree: ScaffoldTreeViewer, group: DG.TreeViewGroup, smiles: string) {
    await scaffoldTree.openEditSketcher(group);
    scaffoldTree.wrapper?.sketcher.setSmiles(smiles);
    const saveButton = scaffoldTree.wrapper?.dialog.root.querySelector('button[name="button-Save"]') as HTMLElement;
    const valueGroup = value(group);
    valueGroup.smiles = smiles;
    saveButton.click();
    await updateVisibleNodes(scaffoldTree);
  }
});
