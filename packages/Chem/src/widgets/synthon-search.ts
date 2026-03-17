import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package, synthonSearchFunc} from '../package';

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
