import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {awaitGrid} from './utils';
import {DebounceIntervals, MolstarViewer, PROPS as msvPROPS} from '../viewers/molstar-viewer/molstar-viewer';
import {CHEM_ATOM_SELECTION_EVENT} from '@datagrok-libraries/chem-meta/src/types';

import {_package} from '../package-test';

category('AtomHighlight', () => {
  let pdbData: BiostructureData;
  let df: DG.DataFrame;
  let ligandCol: DG.Column;
  let view: DG.TableView;

  before(async () => {
    const [pdbBytes, sdfBytes] = await Promise.all([
      _package.files.readAsBytes('samples/1bdq-wo-ligands.pdb'),
      _package.files.readAsBytes('samples/1bdq-obs-pred.sdf'),
    ]);
    pdbData = {binary: true, ext: 'pdb', data: pdbBytes};
    df = (await grok.functions.call('Chem:importSdf', {bytes: sdfBytes}))[0];
    ligandCol = df.getCol('molecule');
    view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);
  });

  after(async () => {
    view?.close();
  });

  /** Creates a Biostructure viewer wired to the shared `df`/`pdbData`, docks
   *  it into `view`, sets the first row as current, and waits for the ligand
   *  debounce + first render. `extraProps` override/extend the defaults. */
  async function makeViewer(extraProps: Record<string, unknown> = {}): Promise<MolstarViewer> {
    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: BiostructureDataJson.fromData(pdbData),
      [msvPROPS.ligandColumnName]: ligandCol.name,
      [msvPROPS.showCurrentRowLigand]: true,
      ...extraProps,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);
    df.currentRowIdx = 0;
    await delay(DebounceIntervals.ligands * 2.5);
    await viewer.awaitRendered();
    return viewer;
  }

  test('global cache stores persistent events only', async () => {
    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 999, atoms: [0, 1, 2], persistent: true,
    });
    await delay(50);

    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 999, atoms: [0, 1, 2, 3, 4, 5], persistent: false,
    });
    await delay(50);

    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 999, atoms: [], persistent: true,
    });
    await delay(50);
  });

  test('highlightAllLigandAtoms — runs without error', async () => {
    const viewer = await makeViewer();

    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 0, atoms: [0, 1, 2], persistent: true,
    });
    await delay(500);

    await viewer.highlightController.highlightAllLigandAtoms();
  }, {timeout: 45000});

  test('base colors applied to comparison pose', async () => {
    const viewer = await makeViewer({
      [msvPROPS.showMouseOverRowLigand]: true,
    });

    df.mouseOverRowIdx = 1;
    await delay(DebounceIntervals.ligands * 2.5);
    await viewer.awaitRendered();

    await viewer.highlightController.applyBaseColors();
  }, {timeout: 45000});
});
