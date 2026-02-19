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

export async function synthonSearchWidget(molecule: string): Promise<DG.Widget> {
  const headerHost = ui.div([]);
  const compsHost = ui.div([ui.loader()], 'd4-flex-wrap chem-viewer-grid chem-search-panel-wrapper');
  const panel = ui.divV([headerHost, compsHost]);

  try {
    const synthonLibrary = await getDefaultSynthonLibrary();
    const df: DG.DataFrame = await scripts.synthonSubstructureSearch(molecule, synthonLibrary, 100);

    ui.empty(compsHost);

    if (!df || df.rowCount === 0) {
      compsHost.appendChild(ui.divText('No matches'));
      return new DG.Widget(panel);
    }

    console.log(df);
    const moleculeCol = df.getCol('smiles');
    const nameCol = df.getCol('name');
    const molCount = Math.min(df.rowCount, MAX_MOLECULES);

    for (let i = 0; i < molCount; i++) {
      const molHost = ui.divV([], {style: {marginBottom: '12px'}});
      molHost.append(drawMolecule(moleculeCol.get(i), WIDTH, HEIGHT, true));
      const name = nameCol.get(i);
      if (name)
        molHost.appendChild(ui.divText(name, {style: {textAlign: 'center', fontSize: '11px'}}));

      compsHost.appendChild(molHost);
    }

    headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
      moleculeCol.semType = DG.SEMTYPE.MOLECULE;
      df.name = 'Synthon Search Results';
      grok.shell.addTableView(df);
    }, 'Open compounds as table'));

    compsHost.style.overflowY = 'auto';
  } catch (e: any) {
    ui.empty(compsHost);
    compsHost.appendChild(ui.divText(`Error: ${e.message ?? 'Unknown error'}`));
  }

  return new DG.Widget(panel);
}
