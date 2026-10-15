import * as DG from 'datagrok-api/dg';

import {Structure} from 'molstar/lib/mol-model/structure';
import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {StructureRef} from 'molstar/lib/mol-plugin-state/manager/structure/hierarchy-state';
import {StateElements} from 'molstar/lib/examples/proteopedia-wrapper/helpers';

import {ligandMapItems} from '../ligand-map';
import type {LigandMap} from './molstar-viewer';

/** What the viewer knows that the plugin does not: which state cells are its data, ligands and
 * binding site. */
export type MolstarViewState = {
  structureRefs: string[] | null,
  ligands: LigandMap,
  bindingSiteRefs: string[],
};

/** The representation over the data structure: the one the viewer applied (`updateView`) if it
 * built, else the first one the load preset made. */
function shownRepresentation(plugin: PluginContext, structure: StructureRef): string {
  const applied = plugin.state.data.cells.get(`${StateElements.SequenceVisual}-${structure.cell.transform.ref}`);
  const cell = applied?.status === 'ok' && applied.obj ? applied : structure.components[0]?.representations[0]?.cell;
  return cell?.transform.params?.type?.name ?? '';
}

/** Adds what the Mol* viewer shows to its base status: the WebGL canvas and the readings a test
 * compares. */
export function molstarStatus(
  base: DG.IWidgetStatus, plugin: PluginContext | null, state: MolstarViewState,
): DG.IWidgetStatus {
  const structureRefs = state.structureRefs ?? [];
  const structure = plugin?.managers.structure.hierarchy.current.structures
    .find((s) => structureRefs.includes(s.cell.transform.ref));
  const ligandRows = Array.from(new Set(ligandMapItems(state.ligands)
    .filter((l) => (l.structureRefs?.length ?? 0) > 0).map((l) => l.rowIdx + 1))).sort((a, b) => a - b);
  const cells = plugin?.state.data.cells;
  const bindingSiteShown = state.bindingSiteRefs.length > 0 &&
    state.bindingSiteRefs.every((ref) => !!cells?.get(ref)?.obj);

  const values: {[name: string]: number | string | boolean} = {
    ...base.values,
    'structure loaded': !!structure?.cell.obj,
    'ligands shown': ligandRows.length,
    'ligand rows': ligandRows.join(', '),
    'layout expanded': plugin?.layout.state.isExpanded ?? false,
    'controls shown': plugin?.layout.state.showControls ?? false,
    'binding site shown': bindingSiteShown,
    'representation': plugin && structure ? shownRepresentation(plugin, structure) : '',
  };
  if (bindingSiteShown) {
    values['binding site atoms'] = state.bindingSiteRefs.reduce((n, ref) => {
      const data = cells?.get(ref)?.obj?.data;
      return n + (data instanceof Structure ? data.elementCount : 0);
    }, 0);
  }

  const canvas = plugin?.canvas3dContext?.canvas;
  return {...base, parts: canvas ? {...base.parts, canvas} : base.parts, values};
}
