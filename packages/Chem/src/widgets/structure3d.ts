import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Chem} from '../scripts-api';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {cruxReady, isCruxRequested} from '../engine/chem-engine';
import {ChemWasm, cruxModuleSync, noteCruxError} from '../engine/crux-module';
import {warnCruxFallback} from '../engine/fallback';
import {getCruxService} from '../crux/crux-searches';
import {getMacromoleculeAssembler} from '../engine/biopolymer';
import {nglStageOnPage} from '@datagrok-libraries/utils/src/ngl-stage';
import type NGL from 'ngl';
import type {Subscription} from 'rxjs';

const WIDTH = 300;
const HEIGHT = 300;
export const NO_3D_VIEW = 'No 3D view: the viewer could not read the structure';

/** crux's `chemEmbed3d` options (its JSON). */
export interface Embed3dOptions {
  etkdg?: 'v1' | 'v3' | 'srv3';
  legacy?: boolean | 'auto';
  randomCoords?: boolean | 'auto';
  maxIterations?: number;
  timeoutMs: number;
  forceField?: 'none' | 'uff';
  seed?: number;
}

/**
 * Crux mode's embedding options, crux's own recommendation (`chemEmbed3dRecommendedOptions`, crux-core decision
 * chem-props, "3D coordinates"): ETKDGv3 without Hs, random starting coordinates for large molecules, at most 200
 * attempts, a UFF clean-up and a 10 s deadline, so a molecule ETKDG cannot embed in time gives up instead of holding
 * its worker.
 */
export function embed3dOptions(): Embed3dOptions {
  return JSON.parse((cruxModuleSync() as unknown as ChemWasm).chemEmbed3dRecommendedOptions());
}

/** Whether `e` is an embedding's timeout (crux's `chemEmbed3d` throws it; the worker stays usable). */
export const isEmbedTimeout = (e: unknown): boolean =>
  String((e as Error)?.message ?? e).startsWith('chemEmbed3d: timed out');

/**
 * In Crux mode, `molecule`'s 3D MOL block by crux's `chemEmbed3d` on a Crux worker ({@link embed3dOptions}); null
 * otherwise, where Crux makes no block and where the call fails, so the script runs ("Chemical Engine",
 * engine/README.md). A timeout throws: the script would not finish either.
 */
export async function cruxEmbed3d(molecule: string, options?: Embed3dOptions): Promise<string | null> {
  // Crux requested: wait for its module (a fresh page has not loaded it yet) rather than run the script meanwhile
  if (!isCruxRequested() || !await cruxReady())
    return null;
  try {
    return await (await getCruxService()).embed3d(molecule, options ?? embed3dOptions());
  } catch (e) {
    if (isEmbedTimeout(e))
      throw e;
    noteCruxError(e);
    warnCruxFallback('3D Structure', e);
    return null;
  }
}

/**
 * `molecule` (SMILES or a MOL block) with 3D coordinates, a MOL block, by the Chemical Engine in effect: Crux's
 * embedding on a Crux worker in Crux mode ({@link cruxEmbed3d}), the SmilesTo3DCoordinates script in RDKit mode and
 * for a molecule Crux cannot read. Throws when the embedding times out. Chem's one entry point for 3D coordinates
 * (`Chem:embed3d`: the 3D Structure pane, Bio's macromolecule pane).
 */
export async function embed3d(molecule: string): Promise<string> {
  return await cruxEmbed3d(molecule) ?? (await Chem.smilesTo3DCoordinates(molecule)).replaceAll('\\n', '\n');
}

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
    sdf = MolfileHandler.getInstance(molecule).z.every((coord) => coord === 0) ? await embed3d(molecule) : molecule;
  } catch (e) {
    return new DG.Widget(ui.divText(isEmbedTimeout(e) ?
      `No 3D coordinates: the embedding took longer than ${embed3dOptions().timeoutMs / 1000} s` :
      'Molecule has no atoms or malformed'));
  }
  // NGL reads V2000 only: a V3000 block (enhanced stereo, as a CXSMILES' `|a:...|`, is written so) goes to V2000
  // where it fits, its 3D coordinates kept (the stereo groups go: the coordinates show the stereo)
  if (sdf.split('\n', 4)[3]?.includes('V3000')) {
    try {
      sdf = await getMacromoleculeAssembler().toV2000(sdf);
    } catch {
      // more than 999 atoms or bonds: the message below
    }
  }
  const stringBlob = new Blob([sdf], {type: 'text/plain'});

  const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});
  nglHost.style.width = `${WIDTH}px`;
  nglHost.style.height = `${HEIGHT}px`;
  // a square where the pane's stylesheet sets the size to `initial` (an accordion pane): the stage sizes its canvas
  // by the host it is made in, which would otherwise be a line high
  nglHost.style.aspectRatio = '1';
  nglHost.style.backgroundColor = 'white';

  await DG.Utils.loadJsCss(['/js/common/ngl_viewer/ngl.js']);
  // what NGL cannot read: a message, not an empty box (a V3000 block loads as a structure without atoms)
  const unreadable = () => nglHost.replaceChildren(ui.divText(NO_3D_VIEW));
  const subs: Subscription[] = [];
  // the stage while the pane is on the page (WebGL contexts are few)
  nglStageOnPage(nglHost, (host) => {
    //@ts-ignore
    const stage = new NGL.Stage(host, {backgroundColor: 'white'}) as NGL.Stage;
    stage.loadFile(stringBlob, {ext: 'sdf'}).then((loaded) => {
      const comp = loaded as NGL.StructureComponent;
      if (!comp?.structure?.atomCount)
        return unreadable();
      stage.handleResize();
      comp.addRepresentation('ball+stick');
      comp.autoView();
    }, unreadable);
    const accPanel = host.closest('.dialog-floating') ? host.closest('.panel-content') as HTMLElement | null : null;
    if (accPanel) {
      subs.push(ui.onSizeChanged(accPanel).subscribe((_) => {
        const w = Math.max(accPanel.clientWidth, 300);
        const h = Math.max(accPanel.clientHeight, 300);
        host.style.width = `${w}px`;
        host.style.height = `${h}px`;
        const waitParentEl = host.parentElement;
        if (waitParentEl?.classList.contains('grok-wait')) {
          waitParentEl.style.width = `${w}px`;
          waitParentEl.style.height = `${h}px`;
        }
        stage.handleResize();
      }));
    }
    return stage;
  }, () => subs.forEach((sub) => sub.unsubscribe())).catch(unreadable);

  return new DG.Widget(nglHost);
}
