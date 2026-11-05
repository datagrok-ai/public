import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package, synthonSearchFunc} from '../package';
import {byEngineAsync} from '../engine/chem-engine';
import {getCruxService} from '../crux/crux-searches';
import {isTrap} from '../crux/crux-service';

/** A library's key in the Crux workers (they keep the libraries they read): its name, length and text's FNV-1a. */
function synthonLibraryKey(name: string, text: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < text.length; i++)
    h = Math.imul(h ^ text.charCodeAt(i), 0x01000193) >>> 0;
  return `${name}:${text.length}:${h.toString(16)}`;
}

/** How many workers build a library's fingerprints at once: each reads the library (~100 MB of its heap for Chem's
 * 5,567 synthons), so a few. */
const SYNTHON_FP_PARTS = 3;
/** The browser's cache of built fingerprints (the script keeps its database on the server likewise). */
const SYNTHON_FP_CACHE = 'chem-synthon-fingerprints';
const synthonFpUrl = (key: string): string => `/chem-synthon-fingerprints/v1/${encodeURIComponent(key)}`;

async function cachedFingerprints(key: string): Promise<Uint8Array | null> {
  try {
    const response = await (await caches.open(SYNTHON_FP_CACHE)).match(synthonFpUrl(key));
    return response ? new Uint8Array(await response.arrayBuffer()) : null;
  } catch {
    return null;
  }
}

async function cacheFingerprints(key: string, bytes: Uint8Array): Promise<void> {
  try {
    await (await caches.open(SYNTHON_FP_CACHE)).put(synthonFpUrl(key), new Response(bytes.slice().buffer));
  } catch {
    // no cache (a private window, no Cache Storage): built again in the next session
  }
}

/** The library's similarity fingerprints on the searching worker: cached, or built across workers (each part after
 * its length, u32 little-endian) and cached. */
async function ensureSynthonFingerprints(key: string, text: string): Promise<void> {
  const service = await getCruxService();
  if (await service.synthonJob<boolean>('synthonHasFingerprints', [key, text]))
    return;
  const cached = await cachedFingerprints(key);
  if (cached) {
    try {
      await service.synthonJob<void>('synthonInstallFingerprints', [key, text, cached]);
      return;
    } catch (e) {
      if (isTrap(e))
        throw e; // a trap: the caller's fallback
      // a cache of another build's format: built again
    }
  }
  const parts = await service.synthonFingerprintParts(key, text, SYNTHON_FP_PARTS);
  const all = new Uint8Array(parts.reduce((sum, p) => sum + 4 + p.length, 0));
  let at = 0;
  for (const p of parts) {
    new DataView(all.buffer).setUint32(at, p.length, true);
    all.set(p, at + 4);
    at += 4 + p.length;
  }
  await service.synthonJob<void>('synthonInstallFingerprints', [key, text, all]);
  await cacheFingerprints(key, all);
}

/**
 * Chem:SynthonSearch: scripts/synthon_search.py's inputs and dataframe. "Chemical engine": in Crux mode on the Crux
 * workers (crux-core's exact port of the script's RDKit 2025.09.3 search: the same hits, order, SMILES and
 * similarities), the script (Chem:SynthonSearchPython) behind it and in RDKit mode.
 */
export async function synthonSearch(molecule: string, synthonLibrary: DG.FileInfo, libraryName: string,
  maxHits: number, searchType: string, similarityCutoff: number, includeSynthons: boolean): Promise<DG.DataFrame> {
  const python = (): Promise<DG.DataFrame> => grok.functions.call('Chem:SynthonSearchPython',
    {molecule, synthonLibrary, libraryName, maxHits, searchType, similarityCutoff, includeSynthons});
  return byEngineAsync('synthonSearch', async () => {
    const text = await synthonLibrary.readAsString();
    const key = synthonLibraryKey(libraryName, text);
    if (searchType === 'similarity' || searchType === 'exact')
      await ensureSynthonFingerprints(key, text);
    const json = await (await getCruxService()).synthonJob<string>('synthonSearch',
      [key, text, molecule, searchType, maxHits, similarityCutoff, includeSynthons]);
    const columns = (JSON.parse(json).columns as {name: string, type: string, values: any[]}[]).map((c) =>
      c.type === 'double' ? DG.Column.fromFloat64Array(c.name, Float64Array.from(c.values)) :
        DG.Column.fromStrings(c.name, c.values));
    const df = DG.DataFrame.fromColumns(columns);
    await grok.data.detectSemanticTypes(df); // as the script's dataframe gets them
    return df;
  }, python);
}

export async function getSynthonSpaces(): Promise<string[]> {
  try {
    const files = await _package.files.list('synthon-data/', false, 'csv');
    return files.map((f) => f.name);
  } catch (_e) {
    return [];
  }
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

  const returnSynthonsInput = ui.input.bool('Include synthons', {value: false, onValueChanged: () => runSearch()});

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
      const resultDf = await synthonSearchFunc(
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
