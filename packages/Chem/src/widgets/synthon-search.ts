import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {scripts} from '../package-api';

const synthonLookupCache: Map<string, Map<string, string>> = new Map();

export async function getSynthonSpaces(): Promise<string[]> {
  try {
    const files = await _package.files.list('synthon-data/', false, 'csv');
    return files.map((f) => f.name);
  } catch (_e) {
    return [];
  }
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
      lookup.set(idCol.get(i).toString(), smilesCol.get(i));
    synthonLookupCache.set(fileName, lookup);
  }
  return synthonLookupCache.get(fileName)!;
}

function buildResultDf(
  searchDf: DG.DataFrame, isSimilarity: boolean,
  returnSynthons: boolean, lookup?: Map<string, string>,
): DG.DataFrame {
  const requiredCols = isSimilarity ? ['smiles', 'similarity'] : ['smiles'];
  const resultDf = searchDf.clone(undefined, requiredCols);
  resultDf.getCol('smiles').semType = DG.SEMTYPE.MOLECULE;

  if (!returnSynthons)
    return resultDf;

  const rowCount = searchDf.rowCount;
  const nameCol = searchDf.getCol('name');

  // Parse names and find max synthon count in a single pass
  const parsed: {reactionId: string; synthonIds: string[]}[] = new Array(rowCount);
  let maxSynthons = 0;
  for (let i = 0; i < rowCount; i++) {
    const name: string = nameCol.get(i) ?? '';
    const parts = name.split(';');
    if (parts.length > 1) {
      const synthonIds = parts.slice(0, -1).map((s) => s.trim());
      parsed[i] = {reactionId: parts[parts.length - 1], synthonIds};
      if (synthonIds.length > maxSynthons)
        maxSynthons = synthonIds.length;
    } else
      parsed[i] = {reactionId: '', synthonIds: []};
  }

  // Create and populate synthon columns
  const reactionIdCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'reaction_id', rowCount);
  const synthonCols: DG.Column[] = [];
  const synthonIdCols: DG.Column[] = [];
  for (let p = 0; p < maxSynthons; p++) {
    const molCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, `synthon_${p + 1}`, rowCount);
    molCol.semType = DG.SEMTYPE.MOLECULE;
    synthonCols.push(molCol);
    synthonIdCols.push(DG.Column.fromType(DG.COLUMN_TYPE.STRING, `synthon_${p + 1}_id`, rowCount));
  }

  for (let i = 0; i < rowCount; i++) {
    const {reactionId, synthonIds} = parsed[i];
    reactionIdCol.set(i, reactionId);
    for (let p = 0; p < synthonIds.length; p++) {
      synthonIdCols[p].set(i, synthonIds[p]);
      synthonCols[p].set(i, lookup!.get(synthonIds[p]) ?? '');
    }
  }

  // Add columns to resultDf
  resultDf.columns.add(reactionIdCol);
  for (let p = 0; p < maxSynthons; p++) {
    resultDf.columns.add(synthonCols[p]);
    resultDf.columns.add(synthonIdCols[p]);
  }

  return resultDf;
}

export async function synthonSearch(
  spaceName: string, molecule: string, maxHits: number,
  searchType: string, similarityCutoff?: number, returnSynthons?: boolean,
): Promise<DG.DataFrame> {
  const isExact = searchType === 'exact';
  const isSimilarity = searchType === 'similarity' || isExact;
  const includeSimilarity = isSimilarity && !isExact;
  const cutoff = isExact ? 1.0 : (similarityCutoff ?? 0.5);

  const fileName = spaceName.endsWith('.csv') ? spaceName : `${spaceName}.csv`;
  const lib = DG.FileInfo.fromString(fileName, await _package.files.readAsText(`synthon-data/${fileName}`));

  let df: DG.DataFrame;
  if (isSimilarity)
    df = await scripts.synthonSimilaritySearch(molecule, lib, fileName, maxHits, cutoff);
  else
    df = await scripts.synthonSubstructureSearch(molecule, lib, fileName, maxHits);

  if (!df || df.rowCount === 0) {
    const emptyDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('smiles', []),
      DG.Column.fromStrings('name', []),
      ...(isSimilarity ? [DG.Column.fromFloat32Array('similarity', new Float32Array(0))] : []),
    ]);
    return buildResultDf(emptyDf, includeSimilarity, !!returnSynthons);
  }

  const lookup = returnSynthons ? await getSynthonLookup(fileName) : undefined;
  return buildResultDf(df, includeSimilarity, !!returnSynthons, lookup);
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

  const libraries = await getSynthonSpaces();
  if (libraries.length === 0)
    return new DG.Widget(ui.divText('No synthon spaces found in synthon-data/'));

  const libraryInput = ui.input.choice('Space', {
    value: libraries[0],
    items: libraries,
    onValueChanged: () => runSearch(),
  });

  const maxHitsInput = ui.input.int('Max hits', {value: 100, onValueChanged: () => runSearch()});

  const returnSynthonsInput = ui.input.bool('Return synthons', {value: false, onValueChanged: () => runSearch()});

  let cutoffInput: DG.InputBase<number | null> | null = null;
  if (isSimilarity) {
    const cutoffValueLabel = ui.divText('0.50', 'ui-input-description');
    cutoffInput = ui.input.slider('Cutoff', {value: 0.5, min: 0, max: 1, step: 0.01, onValueChanged: () => {
      cutoffValueLabel.textContent = (cutoffInput!.value ?? 0.5).toFixed(2);
    }});
    cutoffInput.value = 0.5;
    cutoffInput.root.querySelector('.ui-input-editor')?.addEventListener('mouseup', () => runSearch());
    cutoffInput.addOptions(cutoffValueLabel);
  }

  const inputs: DG.InputBase[] = [libraryInput, maxHitsInput];
  if (cutoffInput)
    inputs.push(cutoffInput);
  inputs.push(returnSynthonsInput);

  const form = ui.form(inputs);
  const panel = ui.divV([form, headerHost, resultsHost]);
  const searchLabel = isSimilarity ? 'Similarity search is running...' : 'Substructure search is running...';

  async function runSearch(): Promise<void> {
    ui.empty(headerHost);
    ui.setUpdateIndicator(resultsHost, true, searchLabel);

    try {
      const resultDf = await synthonSearch(
        libraryInput.value!, molecule, maxHitsInput.value ?? 100,
        searchType, cutoffInput?.value ?? 0.5, returnSynthonsInput.value,
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
