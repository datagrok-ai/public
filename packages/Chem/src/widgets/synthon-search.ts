import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {_package, drawMolecule} from '../package';
import {scripts} from '../package-api';

const WIDTH = 200;
const HEIGHT = 100;
const MAX_MOLECULES = 20;

let cachedSynthonLibrary: DG.FileInfo | null = null;

async function getDefaultSynthonLibrary(): Promise<DG.FileInfo> {
  if (!cachedSynthonLibrary) {
    const csv = await _package.files.readAsText('synthon-data/Syntons_5567.csv');
    cachedSynthonLibrary = DG.FileInfo.fromString('Syntons_5567.csv', csv);
  }
  return cachedSynthonLibrary;
}

async function synthonSearchWidget(
  searchFn: (lib: DG.FileInfo) => Promise<DG.DataFrame>,
  tableName: string,
  labelFn: (df: DG.DataFrame, i: number) => string,
): Promise<DG.Widget> {
  const headerHost = ui.div([]);
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
  const panel = ui.divV([headerHost, compsHost]);

  try {
    const synthonLibrary = await getDefaultSynthonLibrary();
    const df: DG.DataFrame = await searchFn(synthonLibrary);

    ui.empty(compsHost);

    if (!df || df.rowCount === 0) {
      compsHost.appendChild(ui.divText('No matches'));
      return new DG.Widget(panel);
    }

    const moleculeCol = df.getCol('smiles');
    const molCount = Math.min(df.rowCount, MAX_MOLECULES);

    for (let i = 0; i < molCount; i++) {
      const molHost = ui.divV([], {style: {marginBottom: '12px'}});
      molHost.append(drawMolecule(moleculeCol.get(i), WIDTH, HEIGHT, true));
      const label = labelFn(df, i);
      if (label)
        molHost.appendChild(ui.divText(label, {style: {textAlign: 'center', fontSize: '11px'}}));

      compsHost.appendChild(molHost);
    }

    headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
      moleculeCol.semType = DG.SEMTYPE.MOLECULE;
      df.name = tableName;
      grok.shell.addTableView(df);
    }, 'Open compounds as table'));

    compsHost.style.overflowY = 'auto';
  } catch (e: any) {
    ui.empty(compsHost);
    compsHost.appendChild(ui.divText(`Error: ${e.message ?? 'Unknown error'}`));
  }

  return new DG.Widget(panel);
}

export function _synthonSubstructureSearchWidget(molecule: string): Promise<DG.Widget> {
  return synthonSearchWidget(
    (lib) => scripts.synthonSubstructureSearch(molecule, lib, 100),
    'Synthon Substructure Search Results',
    (df, i) => df.getCol('name').get(i) ?? '',
  );
}

export function _synthonSimilaritySearchWidget(molecule: string): Promise<DG.Widget> {
  return synthonSearchWidget(
    (lib) => scripts.synthonSimilaritySearch(molecule, lib, 100, 0.5),
    'Synthon Similarity Search Results',
    (df, i) => {
      const name = df.getCol('name').get(i);
      const sim = df.getCol('similarity').get(i);
      return [name, sim != null ? sim.toFixed(2) : null].filter(Boolean).join(' | ');
    },
  );
}
