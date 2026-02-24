import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {scripts} from '../package-api';

const synthonLibraryCache: Map<string, DG.FileInfo> = new Map();
const synthonLookupCache: Map<string, Map<string, string>> = new Map();
let availableLibraries: string[] | null = null;

export async function listSynthonLibraries(): Promise<string[]> {
  if (!availableLibraries) {
    try {
      const files = await _package.files.list('synthon-data/', false, 'csv');
      availableLibraries = files.map((f) => f.name);
    } catch (_e) {
      availableLibraries = [];
    }
  }
  return availableLibraries;
}

async function getSynthonLibrary(fileName: string): Promise<DG.FileInfo> {
  if (!synthonLibraryCache.has(fileName)) {
    const csv = await _package.files.readAsText(`synthon-data/${fileName}`);
    synthonLibraryCache.set(fileName, DG.FileInfo.fromString(fileName, csv));
  }
  return synthonLibraryCache.get(fileName)!;
}

async function getSynthonLookup(fileName: string): Promise<Map<string, string>> {
  if (!synthonLookupCache.has(fileName)) {
    const df = await _package.files.readCsv(`synthon-data/${fileName}`);
    await grok.data.detectSemanticTypes(df);
    const lookup = new Map<string, string>();
    const smilesCol = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
    if (!smilesCol)
      throw new Error(`${fileName} doesn't contain molecule column`);
    const idCol = df.col('syntnon #') ?? df.getCol('synton_id');
    for (let i = 0; i < df.rowCount; i++)
      lookup.set(idCol.get(i), smilesCol.get(i));
    synthonLookupCache.set(fileName, lookup);
  }
  return synthonLookupCache.get(fileName)!;
}

function buildResultDf(
  searchDf: DG.DataFrame, isSimilarity: boolean, lookup: Map<string, string>,
): DG.DataFrame {
  const rowCount = searchDf.rowCount;
  const nameCol = searchDf.getCol('name');

  let maxPrecursors = 0;
  for (let i = 0; i < rowCount; i++) {
    const name: string = nameCol.get(i) ?? '';
    const parts = name.split(';');
    const count = parts.length > 1 ? parts.length - 1 : 0;
    if (count > maxPrecursors)
      maxPrecursors = count;
  }

  const cols: DG.Column[] = [];
  const structureCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'product', rowCount);
  structureCol.semType = DG.SEMTYPE.MOLECULE;
  cols.push(structureCol);

  const reactionIdCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'reaction_id', rowCount);
  cols.push(reactionIdCol);

  let simCol: DG.Column | null = null;
  if (isSimilarity) {
    simCol = DG.Column.fromType(DG.COLUMN_TYPE.FLOAT, 'similarity', rowCount);
    cols.push(simCol);
  }

  const precursorCols: DG.Column[] = [];
  const precursorIdCols: DG.Column[] = [];
  for (let p = 0; p < maxPrecursors; p++) {
    const molCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, `precursor_${p + 1}`, rowCount);
    molCol.semType = DG.SEMTYPE.MOLECULE;
    const idCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, `precursor_${p + 1}_id`, rowCount);
    precursorCols.push(molCol);
    precursorIdCols.push(idCol);
    cols.push(molCol);
    cols.push(idCol);
  }

  const smilesCol = searchDf.getCol('smiles');
  const srcSimCol = isSimilarity ? searchDf.getCol('similarity') : null;

  for (let i = 0; i < rowCount; i++) {
    structureCol.set(i, smilesCol.get(i));
    const name: string = nameCol.get(i) ?? '';
    const parts = name.split(';');

    if (parts.length > 1) {
      reactionIdCol.set(i, parts[parts.length - 1]);
      for (let p = 0; p < parts.length - 1; p++) {
        const synthonId = parts[p].trim();
        precursorIdCols[p].set(i, synthonId);
        precursorCols[p].set(i, lookup.get(synthonId) ?? '');
      }
    }

    if (srcSimCol)
      simCol!.set(i, srcSimCol.get(i));
  }

  return DG.DataFrame.fromColumns(cols);
}

export async function synthonSearch(
  spaceName: string, molecule: string, maxHits: number,
  searchType: string, similarityCutoff?: number,
): Promise<DG.DataFrame> {
  const fileName = spaceName.endsWith('.csv') ? spaceName : `${spaceName}.csv`;
  const lib = await getSynthonLibrary(fileName);
  const lookup = await getSynthonLookup(fileName);

  let df: DG.DataFrame;
  if (searchType === 'similarity')
    df = await scripts.synthonSimilaritySearch(molecule, lib, fileName, maxHits, similarityCutoff ?? 0.5);
  else
    df = await scripts.synthonSubstructureSearch(molecule, lib, fileName, maxHits);

  if (!df || df.rowCount === 0) {
    return buildResultDf(DG.DataFrame.fromColumns([
      DG.Column.fromStrings('smiles', []),
      DG.Column.fromStrings('name', []),
      ...(searchType === 'similarity' ? [DG.Column.fromFloat32Array('similarity', new Float32Array(0))] : []),
    ]), searchType === 'similarity', lookup);
  }

  return buildResultDf(df, searchType === 'similarity', lookup);
}

// Widget

async function synthonSearchWidget(
  molecule: string, searchType: 'substructure' | 'similarity', tableName: string,
): Promise<DG.Widget> {
  const isSimilarity = searchType === 'similarity';

  const resultsHost = ui.div([]);
  resultsHost.style.position = 'relative';
  resultsHost.style.minHeight = '50px';
  const headerHost = ui.div([]);

  const libraries = await listSynthonLibraries();
  if (libraries.length === 0)
    return new DG.Widget(ui.divText('No synthon spaces found in synthon-data/'));

  const libraryInput = ui.input.choice('Space', {
    value: libraries[0],
    items: libraries,
    onValueChanged: () => runSearch(),
  });

  const maxHitsInput = ui.input.int('Max hits', {value: 100, onValueChanged: () => runSearch()});

  const controlsDiv = ui.divV([libraryInput.root, maxHitsInput.root]);

  let cutoffInput: DG.InputBase<number | null> | null = null;
  if (isSimilarity) {
    const cutoffValueLabel = ui.divText('0.50', 'ui-input-description');
    cutoffInput = ui.input.slider('Cutoff', {value: 0.5, min: 0, max: 1, step: 0.01, onValueChanged: () => {
      cutoffValueLabel.textContent = (cutoffInput!.value ?? 0.5).toFixed(2);
    }});
    cutoffInput.value = 0.5;
    cutoffInput.root.querySelector('.ui-input-editor')?.addEventListener('mouseup', () => runSearch());
    cutoffInput.addOptions(cutoffValueLabel);
    controlsDiv.appendChild(cutoffInput.root);
  }

  const panel = ui.divV([controlsDiv, headerHost, resultsHost]);
  const searchLabel = isSimilarity ? 'Similarity search is running...' : 'Substructure search is running...';

  async function runSearch(): Promise<void> {
    ui.empty(headerHost);
    ui.setUpdateIndicator(resultsHost, true, searchLabel);

    try {
      const resultDf = await synthonSearch(
        libraryInput.value!, molecule, maxHitsInput.value ?? 100,
        searchType, cutoffInput?.value ?? 0.5,
      );
      resultDf.name = tableName;

      ui.setUpdateIndicator(resultsHost, false);
      ui.empty(resultsHost);

      if (resultDf.rowCount === 0) {
        resultsHost.appendChild(ui.divText('No matches'));
        return;
      }

      const grid = resultDf.plot.grid();
      grid.root.style.width = '100%';
      grid.root.style.height = '400px';
      resultsHost.appendChild(grid.root);

      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        grok.shell.addTableView(resultDf);
      }, 'Open compounds as table'));
    } catch (e: any) {
      ui.setUpdateIndicator(resultsHost, false);
      ui.empty(resultsHost);
      resultsHost.appendChild(ui.divText(`Error: ${e.message ?? 'Unknown error'}`));
    }
  }

  runSearch();
  return new DG.Widget(panel);
}

export function _synthonSubstructureSearchWidget(molecule: string): Promise<DG.Widget> {
  return synthonSearchWidget(molecule, 'substructure', 'Synthon Substructure Search Results');
}

export function _synthonSimilaritySearchWidget(molecule: string): Promise<DG.Widget> {
  return synthonSearchWidget(molecule, 'similarity', 'Synthon Similarity Search Results');
}
