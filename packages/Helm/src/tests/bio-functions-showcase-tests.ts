import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

import {initHelmMainPackage} from './utils';

import {_package} from '../package-test';

const ROWS = 53;

category('Bio functions: HELM showcase', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    await initHelmMainPackage();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  async function openShowcase(): Promise<{df: DG.DataFrame, tv: DG.TableView, helm: DG.Column<string>}> {
    const df = await _package.files.readCsv('samples/helm-showcase.csv');
    const tv = grok.shell.addTableView(df);
    await grok.data.detectSemanticTypes(df);
    expectHelmColumn(df);
    return {df, tv, helm: df.getCol('HELM')};
  }

  function expectHelmColumn(df: DG.DataFrame): void {
    const helm = df.getCol('HELM');
    expect(df.rowCount, ROWS);
    expect(helm.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(helm.meta.units, NOTATION.HELM);
    expect(helm.getTag(DG.TAGS.CELL_RENDERER), 'helm');
  }

  function findViewer(tv: DG.TableView, type: string): DG.Viewer {
    const viewer = Array.from(tv.viewers).find((v) => v.type === type);
    expect(viewer != null, true, `${type} viewer not added`);
    return viewer!;
  }

  async function awaitReading(widget: DG.Widget, name: string, value: any): Promise<void> {
    await awaitCheck(() => widget.getWidgetStatus().values?.[name] === value,
      `"${name}" is ${widget.getWidgetStatus().values?.[name]}, expected ${value}`, 20000);
  }

  test('Composition', async () => {
    const {df, tv} = await openShowcase();
    await grok.functions.call('Bio:compositionAnalysis');
    const wl = findViewer(tv, 'WebLogo');
    expect(wl.props['sequenceColumnName'], 'HELM');
    await awaitReading(wl, 'rows shown', ROWS);
    await awaitCheck(() => wl.getWidgetStatus().hitAreas?.['position 1'] != null,
      'WebLogo has no "position 1" area', 10000);
    expectHelmColumn(df);
  });

  test('Sequence Space', async () => {
    const {df, tv, helm} = await openShowcase();
    await grok.functions.call('Bio:sequenceSpaceTopMenu', {table: df, molecules: helm,
      methodName: DimReductionMethods.UMAP, similarityMetric: MmDistanceFunctionsNames.HAMMING,
      plotEmbeddings: true, clusterEmbeddings: true});
    const x = df.getCol('Embed_X_1');
    const y = df.getCol('Embed_Y_1');
    expect(x.stats.missingValueCount, 0);
    expect(y.stats.missingValueCount, 0);
    expect(x.stats.uniqueCount >= 10, true, `Embed_X_1 has ${x.stats.uniqueCount} distinct values`);
    expect(findViewer(tv, DG.VIEWER.SCATTER_PLOT).props['xColumnName'], 'Embed_X_1');
    expectHelmColumn(df);
  });

  test('Convert Sequence Notation', async () => {
    const {df, helm} = await openShowcase();
    const res: DG.Column<string> = await grok.functions.call('Bio:convertNotation',
      {table: df, sequence: helm, targetNotation: NOTATION.SEPARATOR, separator: '-'});
    expect(df.col(res.name) != null, true);
    expect(res.meta.units, NOTATION.SEPARATOR);
    expect(res.get(1), 'A-C-D-E-F-G-H-I-K-L');
    expectHelmColumn(df);
  });

  test('Extract Region', async () => {
    const {df, helm} = await openShowcase();
    await grok.functions.call('Bio:getRegionTopMenu',
      {table: df, sequence: helm, start: '1', end: '2', name: 'region 1-2'});
    const region = df.getCol('region 1-2');
    expect(region.meta.units, NOTATION.HELM);
    expect(region.get(1), 'PEPTIDE1{A.C}$$$$');
    expect(region.get(4), 'PEPTIDE1{A.A}$$$$');
    expectHelmColumn(df);
  });

  test('Similarity Search', async () => {
    const {df, tv} = await openShowcase();
    df.currentRowIdx = 0;
    await grok.functions.call('Bio:similaritySearchTopMenu');
    const viewer = findViewer(tv, 'Sequence Similarity Search');
    await awaitReading(viewer, 'neighbours', 11);
    expect(viewer.getWidgetStatus().values?.['source column'], 'HELM');
    expect(viewer.getWidgetStatus().values?.['target row'], 0);
    expectHelmColumn(df);
  });

  test('Diversity Search', async () => {
    const {df, tv} = await openShowcase();
    await grok.functions.call('Bio:diversitySearchTopMenu');
    const viewer = findViewer(tv, 'Sequence Diversity Search');
    await awaitReading(viewer, 'subset size', 10);
    expect(viewer.getWidgetStatus().values?.['source column'], 'HELM');
    expect(viewer.getWidgetStatus().values?.['distinct sequences'], 10);
    expectHelmColumn(df);
  });

  test('Subsequence Search', async () => {
    const {df, tv, helm} = await openShowcase();
    await grok.functions.call('Bio:SubsequenceSearchTopMenu', {macromolecules: helm});
    const fg = tv.getFiltersGroup({createDefaultFilters: false});
    await awaitCheck(() => fg.filters.length === 1, `${fg.filters.length} filters, expected 1`, 10000);
    expect(fg.getStates('HELM', 'Bio:bioSubstructureFilter').length, 1);
    expectHelmColumn(df);
  });

  test('Split to Monomers', async () => {
    const {df, helm} = await openShowcase();
    await grok.functions.call('Bio:splitToMonomersTopMenu', {table: df, sequence: helm});
    expect(df.getCol('1').semType, 'Monomer');
    expect(df.getCol('1').get(1), 'A');
    expect(df.getCol('10').get(1), 'L');
    expectHelmColumn(df);
  });
});
