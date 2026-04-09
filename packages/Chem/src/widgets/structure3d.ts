import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Chem} from '../scripts-api';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {CHEM_INTERACTIVE_SELECTION_EVENT} from '../constants';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import type NGL from 'ngl';

const WIDTH = 300;
const HEIGHT = 300;

export async function structure3dWidget(molecule: string): Promise<DG.Widget> {
  const rdKitModule = getRdKitModule();
  try {
    if (!DG.chem.isMolBlock(molecule))
      molecule = _convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  let sdf: string;
  try {
    sdf = MolfileHandler.getInstance(molecule).z.every((coord) => coord === 0) ?
      (await Chem.smilesTo3DCoordinates(molecule)).replaceAll('\\n', '\n') : molecule;
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule has no atoms or malformed'));
  }
  const stringBlob = new Blob([sdf], {type: 'text/plain'});

  const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});
  nglHost.style.width = `${WIDTH}px`;
  nglHost.style.height = `${HEIGHT}px`;
  nglHost.style.backgroundColor = 'white';

  await DG.Utils.loadJsCss(['/js/common/ngl_viewer/ngl.js']);
  //@ts-ignore
  const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'}) as NGL;

  // The NGL component is captured inside the .then() callback so that
  // loadFile completes AFTER the widget host is returned and mounted in
  // the DOM (NGL needs a valid-size container to render). Using await
  // here would block the widget return and NGL would render into a
  // zero-size detached element — breaking the 3D view entirely.
  let comp: NGL.StructureComponent | null = null;
  let highlightRepr: any = null;
  // Mapping from 2D heavy-atom index (0-based, matching RDKit) to 3D NGL
  // atom serial number. The 3D SDF typically includes explicit hydrogens
  // that the 2D depiction treats as implicit, so 2D atom index i does NOT
  // correspond to NGL serial i+1. We build this mapping once at load time
  // by iterating the NGL structure and recording the serial numbers of
  // non-hydrogen atoms in order.
  let heavyAtomSerials: number[] = [];

  stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(c: NGL.StructureComponent) {
    comp = c;
    stage.setSize(WIDTH, HEIGHT);
    comp.addRepresentation('ball+stick');
    comp.autoView();

    // Build the heavy-atom mapping: iterate all atoms in NGL's structure
    // and collect the serial numbers of non-H atoms in order. The i-th
    // entry corresponds to 2D atom index i.
    heavyAtomSerials = [];
    (comp as any).structure.eachAtom(function(ap: any) {
      if (ap.element !== 'H')
        heavyAtomSerials.push(ap.serial);
    });

    // ---- Restore existing picks when the widget is recreated -----------
    // The NGL context-panel widget is destroyed and recreated every time
    // the user clicks a different row. If the user previously picked atoms
    // in the 2D grid for THIS row, the picks are still stored on the
    // molecule column's temp store (via the ISubstructProvider). Read them
    // and apply the highlight immediately so the 3D view stays in sync
    // with the 2D grid, which does persist its highlights across row
    // switches.
    try {
      const df = grok.shell.tv?.dataFrame;
      const rowIdx = df?.currentRowIdx;
      if (df && rowIdx !== undefined && rowIdx >= 0) {
        const molCol = df.columns.toList().find(
          (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE);
        if (molCol) {
          const providers = ((molCol.temp as any)?.[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as any[];
          const picker = providers.find(
            (p: any) => p.__atomPicker && p.__rowIdx === rowIdx);
          if (picker?.__atoms?.size > 0) {
            const atoms = [...picker.__atoms] as number[];
            const serials = atoms
              .filter((i) => i >= 0 && i < heavyAtomSerials.length)
              .map((i) => heavyAtomSerials[i]);
            if (serials.length > 0) {
              const sele = '@' + serials.join(',');
              highlightRepr = comp.addRepresentation('ball+stick', {
                sele: sele,
                colorScheme: 'uniform',
                colorValue: 0xFFFF00,
                radiusScale: 1.8,
                opacity: 0.85,
              });
            }
          }
        }
      }
    } catch (err) {
      console.error('structure3d: initial highlight restore failed', err);
    }
  });

  // ---- 2D ↔ 3D bridge: interactive atom highlighting --------------------
  // When the Chem package's in-grid atom picker fires a selection event,
  // overlay a yellow highlight on the corresponding atoms in this NGL
  // viewer. NGL's atom serial numbers are 1-based (matching SDF atom
  // order), while the event carries 0-based RDKit indices — offset by +1.
  //
  // The highlight is rendered as a second 'ball+stick' representation with
  // a uniform yellow color and a slightly larger radius, layered on top of
  // the base representation. Clearing the selection removes the overlay.
  const eventSub = grok.events.onCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT)
    .subscribe((args: any) => {
      if (!comp) return; // NGL hasn't finished loading yet
      try {
        // Remove the previous highlight overlay (if any)
        if (highlightRepr) {
          comp.removeRepresentation(highlightRepr);
          highlightRepr = null;
        }
        // Apply new highlight if atoms were selected
        if (args?.atoms?.length > 0 && heavyAtomSerials.length > 0) {
          // Map 2D heavy-atom indices to 3D NGL serial numbers.
          // If the 3D SDF has explicit H's (which smilesTo3DCoordinates
          // typically adds), the serial numbers are NOT just index+1 —
          // H's are interleaved and shift everything. The heavyAtomSerials
          // array built at load time gives us the correct mapping.
          const serials = args.atoms
            .filter((i: number) => i >= 0 && i < heavyAtomSerials.length)
            .map((i: number) => heavyAtomSerials[i]);
          if (serials.length === 0) return;
          // NGL selection language: @serial1,serial2,...
          const sele = '@' + serials.join(',');
          highlightRepr = comp.addRepresentation('ball+stick', {
            sele: sele,
            colorScheme: 'uniform',
            colorValue: 0xFFFF00, // yellow, matching the 2D picker
            radiusScale: 1.8, // slightly larger to stand out
            opacity: 0.85,
          });
        }
      } catch (err) {
        console.error('structure3d: highlight failed', err);
      }
    });

  // Clean up the event subscription when the host element is removed from
  // the DOM (e.g. when the user clicks a different row and the panel
  // recreates the widget).
  const observer = new MutationObserver(() => {
    if (!document.body.contains(nglHost)) {
      eventSub.unsubscribe();
      observer.disconnect();
    }
  });
  // Observe the host's parent (or body) for child-list changes
  ui.tools.waitForElementInDom(nglHost).then(() => {
    const parent = nglHost.parentElement ?? document.body;
    observer.observe(parent, {childList: true});

    if (nglHost.closest('.dialog-floating')) {
      const accPanel = nglHost.closest('.panel-content') as HTMLElement;
      if (accPanel) {
        ui.onSizeChanged(accPanel).subscribe((_) => {
          const w = Math.max(accPanel.clientWidth, 300);
          const h = Math.max(accPanel.clientHeight, 300);
          nglHost.style.width = `${w}px`;
          nglHost.style.height = `${h}px`;
          const waitParentEl = nglHost.parentElement;
          if (waitParentEl?.classList.contains('grok-wait')) {
            waitParentEl.style.width = `${w}px`;
            waitParentEl.style.height = `${h}px`;
          }
          stage.handleResize();
        });
      }
    }
  });

  return new DG.Widget(nglHost);
}
