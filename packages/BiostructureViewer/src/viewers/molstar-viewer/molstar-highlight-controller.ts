/* eslint-disable new-cap */
/**
 * Controller for the 2D↔3D atom-highlighting bridge, extracted from molstar-viewer.ts.
 * Reads viewer state through the narrow HighlightHost interface to avoid a
 * circular compile-time dependency on the full MolstarViewer class.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {MolScriptBuilder} from 'molstar/lib/mol-script/language/builder';
import {Script} from 'molstar/lib/mol-script/script';
import {Structure, StructureSelection, StructureElement} from 'molstar/lib/mol-model/structure';
import {Color} from 'molstar/lib/mol-util/color';
import {
  setStructureOverpaint, clearStructureOverpaint,
} from 'molstar/lib/mol-plugin-state/helpers/structure-overpaint';
import {StructureComponentRef} from 'molstar/lib/mol-plugin-state/manager/structure/hierarchy-state';

import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {
  computeSerials, selectionCacheKey,
} from './molstar-highlight-utils';
import {AtomMapping3D, AtomPickerProvider, CHEM_MOL3D_SELECTION_EVENT, Mol3DHoverEventArgs}
  from '@datagrok-libraries/chem-meta/src/types';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {MolstarViewer, type LigandMap} from './molstar-viewer';

/** Narrow view of `MolstarViewer` that the highlight controller needs. */
export interface HighlightHost {
  getPlugin(): PluginContext | undefined;
  getLigands(): LigandMap;
  getDataFrame(): DG.DataFrame | null;
  getLigandColumnName(): string | undefined;
  readonly viewSyncer: PromiseSyncer;
  readonly logger: ILogger;
}

/** Loaded ligand row bundled with its Molstar state refs. */
type LigandInfo = {rowIdx: number, structureRefs: string[] | null};

/** Internal record built during highlight dispatch. */
type LigandHighlight = {serials: number[], isCurrent: boolean, dataRef: string | null};

export class MolstarHighlightController {
  /** Dedup guard for Molstar's per-pixel hover events. */
  private _lastHoverFired: {
    rowIdx: number, atomSerial: number | null, mode: string,
  } | null = null;

  constructor(private readonly host: HighlightHost) {}

  // -- 3D→2D hover bridge ----------------------------------------------------

  /** Call before re-subscribing Molstar hover so the first post-rebuild event is not deduped. */
  resetHoverDedup(): void {
    this._lastHoverFired = null;
  }

  /** Emits `CHEM_MOL3D_HOVER_EVENT` with dedup. `atomSerial: null` means the
   *  cursor left all ligand atoms — the 2D side uses this to clear its
   *  preview highlight. */
  fireMol3DHover(
    rowIdx: number, atomSerial: number | null, mode: 'preview' | 'paint' | 'erase',
  ): void {
    const last = this._lastHoverFired;
    if (last && last.rowIdx === rowIdx && last.atomSerial === atomSerial && last.mode === mode)
      return;
    this._lastHoverFired = {rowIdx, atomSerial, mode};
    const args: Mol3DHoverEventArgs = {
      mol3DColumnName: this.host.getLigandColumnName()!,
      rowIdx, atom3DSerial: atomSerial, mode,
    };
    grok.events.fireCustomEvent(CHEM_MOL3D_SELECTION_EVENT, args);
  }

  // -- Replay after viewer rebuilds ------------------------------------------

  /** Queues a highlight replay through viewSyncer so it lands after any pending
   *  structure rebuilds (overpaint applied before a rebuild is discarded). */
  replayHighlightIfCached(): void {
    this.host.viewSyncer.sync('replayHighlight', async () => {
      try {
        const plugin = this.host.getPlugin();
        const comps = this._getAllComponents();
        if (comps && plugin)
          await clearStructureOverpaint(plugin, comps);
      } catch {
        /* best-effort */
      }
      await this.applyBaseColors();
      await this.highlightAllLigandAtoms();
    });
  }

  // -- Base colors -----------------------------------------------------------

  /** Tints comparison-pose carbons grey so multiple poses are visually distinct. */
  async applyBaseColors(): Promise<void> {
    const plugin = this.host.getPlugin();
    if (!plugin) return;

    const loadedLigands = this._collectLoadedLigands();
    if (loadedLigands.length < 2) return;

    const allComponents = this._getAllComponents();
    if (!allComponents) return;

    const currentRowIdx = this.host.getDataFrame()?.currentRowIdx ?? -1;

    for (const ligand of loadedLigands) {
      if (ligand.rowIdx === currentRowIdx) continue;

      const structure = this._getStructureForLigand(ligand);
      if (!structure) continue;

      const targetStructure = structure;
      await setStructureOverpaint(
        plugin, allComponents, Color(0x616161),
        async (structureData: Structure) => {
          if (structureData !== targetStructure)
            return StructureElement.Loci.none(structureData);
          const query = MolScriptBuilder.struct.generator.atomGroups({
            'atom-test': MolScriptBuilder.core.rel.eq([
              MolScriptBuilder.acp('elementSymbol'),
              MolScriptBuilder.es('C'),
            ]),
          });
          const sel = Script.getStructureSelection(query, structureData);
          return StructureSelection.toLociWithSourceUnits(sel);
        },
      );
    }
  }

  // -- Multi-ligand highlight (main entry point) -----------------------------

  /** Applies overpaint for ALL loaded ligands. Resolves atoms from (in order):
   *  liveEvent, the persistent selection cache, picker providers on the mol column.
   *  No _highlightInProgress guard — must always run to completion. */
  async highlightAllLigandAtoms(
    liveEvent?: {rowIdx: number, atoms: number[], mapping3D: AtomMapping3D | null} | null,
  ): Promise<void> {
    const plugin = this.host.getPlugin();
    if (!plugin) return;

    const loadedLigands = this._collectLoadedLigands();
    const df = this.host.getDataFrame();
    const currentRowIdx = df?.currentRowIdx ?? -1;
    const molCol = df?.columns.toList().find(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE);
    const dfId = df?.id ?? '';
    const dfName = df?.name ?? '';
    const colName = molCol?.name ?? '';
    const cache = MolstarViewer.selectionCache;

    const highlights: LigandHighlight[] = [];
    for (const ligand of loadedLigands) {
      let atoms: number[] = [];
      let mapping3D: AtomMapping3D | null = null;

      if (liveEvent && liveEvent.rowIdx === ligand.rowIdx && liveEvent.atoms.length > 0) {
        atoms = liveEvent.atoms;
        mapping3D = liveEvent.mapping3D;
      }
      if (atoms.length === 0) {
        const key = selectionCacheKey(dfId, dfName, colName, ligand.rowIdx);
        const cached = cache.get(key);
        if (cached && cached.atoms.length > 0) {
          atoms = cached.atoms;
          mapping3D = cached.mapping3D ?? null;
        }
      }
      if (atoms.length === 0) {
        atoms = this._getHighlightedAtomsFromColumn(ligand.rowIdx, colName);
        const key = selectionCacheKey(dfId, dfName, colName, ligand.rowIdx);
        const cached = cache.get(key);
        mapping3D = cached?.mapping3D ?? null;
      }
      if (atoms.length === 0) continue;

      const structure = this._getStructureForLigand(ligand);
      const serials = computeSerials(atoms, mapping3D, structure);
      if (serials.length === 0) continue;

      // Scope to the ligand's data subtree (structureRefs[0]) to avoid painting the receptor.
      highlights.push({
        serials,
        isCurrent: ligand.rowIdx === currentRowIdx,
        dataRef: ligand.structureRefs?.[0] ?? null,
      });
    }

    const allComponents = this._getAllComponents();
    if (!allComponents) return;

    await clearStructureOverpaint(plugin, allComponents);
    await this.applyBaseColors();

    for (const h of highlights) {
      const color = h.isCurrent ? Color(0xFFFF00) : Color(0x00BCD4);
      await setStructureOverpaint(
        plugin, allComponents, color,
        async (structureData: Structure) => {
          if (h.dataRef && !this._structureDescendsFromRef(plugin, structureData, h.dataRef))
            return StructureElement.Loci.none(structureData);

          const atomSet = MolScriptBuilder.set(...h.serials);
          const query = MolScriptBuilder.struct.generator.atomGroups({
            'atom-test': MolScriptBuilder.core.set.has([
              atomSet, MolScriptBuilder.ammp('id'),
            ]),
          });
          const sel = Script.getStructureSelection(query, structureData);
          return StructureSelection.toLociWithSourceUnits(sel);
        },
      );
    }
  }

  // -- Internals -------------------------------------------------------------

  private _collectLoadedLigands(): LigandInfo[] {
    const loaded: LigandInfo[] = [];
    const seen = new Set<number>();
    const lig = this.host.getLigands();
    if (!lig) return loaded;

    const push = (item: LigandInfo | null | undefined): void => {
      if (!item) return;
      if (item.rowIdx < 0 || seen.has(item.rowIdx)) return;
      loaded.push(item);
      seen.add(item.rowIdx);
    };

    push(lig.current);
    push(lig.hovered);
    if (lig.selected)
      for (const sel of lig.selected) push(sel);

    return loaded;
  }

  private _getStructureForLigand(ligand: LigandInfo): Structure | undefined {
    const plugin = this.host.getPlugin();
    if (!plugin) return undefined;
    if (!ligand.structureRefs || ligand.structureRefs.length < 4) return undefined;
    const cell = plugin.state.data.cells.get(ligand.structureRefs[3]);
    return cell?.obj?.data;
  }

  private _getAllComponents(): StructureComponentRef[] | null {
    const plugin = this.host.getPlugin();
    if (!plugin) return null;
    const structures = plugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) return null;
    const comps = structures.flatMap((s) => s.components ?? []);
    return comps.length > 0 ? comps : null;
  }

  /** Walks up to 6 levels of Molstar's state cell tree to check lineage.
   *  Fail-open (returns true on error) to preserve prior overpaint behavior. */
  private _structureDescendsFromRef(
    plugin: PluginContext, structureData: Structure, dataRef: string,
  ): boolean {
    try {
      for (const [, cell] of plugin.state.data.cells) {
        if (cell.obj?.data !== structureData) continue;
        let cur: any = cell;
        for (let d = 0; cur && d < 6; d++) {
          if (cur.transform?.ref === dataRef) return true;
          const pRef = cur.transform?.parent;
          cur = pRef ? plugin.state.data.cells.get(pRef) : null;
        }
        return false;
      }
      return false;
    } catch {
      return true; // fail-open: preserve prior behavior on walk error
    }
  }

  private _getHighlightedAtomsFromColumn(rowIdx: number, molColName?: string): number[] {
    const df = this.host.getDataFrame();
    if (!df) return [];

    let molCol: DG.Column | null = null;
    if (molColName)
      molCol = df.col(molColName);
    if (!molCol) {
      molCol = df.columns.toList().find(
        (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE &&
          ((c.temp as Record<string, unknown>)?.[ChemTemps.SUBSTRUCT_PROVIDERS] as
            AtomPickerProvider[] ?? []).some((p) => p.__atomPicker)) ?? null;
    }
    if (!molCol) return [];

    const providers = ((molCol.temp as Record<string, unknown>)?.[ChemTemps.SUBSTRUCT_PROVIDERS] ??
      []) as AtomPickerProvider[];
    const picker = providers.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    if (!picker?.__atoms || picker.__atoms.size === 0) return [];
    return [...picker.__atoms];
  }
}
