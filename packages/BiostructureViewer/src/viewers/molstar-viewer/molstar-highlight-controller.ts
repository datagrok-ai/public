/**
 * Controller for the interactive 2D↔3D atom-highlighting bridge.
 *
 * Extracted from `molstar-viewer.ts` (which was approaching 1800 lines) so
 * the viewer class can stay focused on Molstar lifecycle while this module
 * owns:
 *  - overpaint application state (concurrency guard, last-hover dedup);
 *  - the 3D→2D `CHEM_MOL3D_HOVER_EVENT` emitter;
 *  - per-ligand base-color assignment;
 *  - the `highlightAllLigandAtoms` replay pipeline used when Molstar
 *    rebuilds its scene;
 *  - `highlightLigandAtoms` / `clearLigandAtomHighlight` — the single-pose
 *    overpaint API.
 *
 * The controller reads viewer state through a narrow {@link HighlightHost}
 * interface so it does NOT take a circular compile-time dependency on the
 * full `MolstarViewer` class. The viewer builds a host facade in its
 * constructor and passes it in.
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

import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {
  AtomMapping3D, CHEM_MOL3D_HOVER_EVENT, Mol3DHoverEventArgs,
  computeSerials, getSelectionCache, selectionCacheKey,
} from './molstar-highlight-utils';
import type {LigandMap} from './molstar-viewer';

/** Narrow view of `MolstarViewer` that the highlight controller needs. Only
 *  getter functions for fields that can change at runtime, plus direct
 *  references to the always-available viewSyncer/logger. */
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
  /** Guard against concurrent async highlight updates on the single-pose
   *  `highlightLigandAtoms` path. The multi-pose `highlightAllLigandAtoms`
   *  intentionally does NOT check this flag — see its docstring. */
  private _highlightInProgress = false;

  /** Last-fired 3D→2D hover tuple, used to dedup Molstar's per-mousemove
   *  pixel firing. `null` means no atom is currently hovered. */
  private _lastHoverFired: {
    rowIdx: number, atomSerial: number | null, mode: string,
  } | null = null;

  constructor(private readonly host: HighlightHost) {}

  // -- 3D→2D hover bridge ----------------------------------------------------

  /** Clears the last-hover dedup state. The viewer calls this before
   *  re-subscribing Molstar's hover behavior so the first post-rebuild
   *  hover event is not suppressed as a duplicate. */
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
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, args);
  }

  // -- Replay after viewer rebuilds ------------------------------------------

  /** Queues a highlight replay after the viewer finishes rebuilding its
   *  scene. Reads directly from the molecule column's temp providers (the
   *  source of truth) so highlights survive cell switches even when no
   *  fresh `CHEM_SELECTION_EVENT` arrives. */
  replayHighlightIfCached(): void {
    this.host.logger.debug('[molstar-picker] scheduling replay after buildViewLigands');
    // Run through the viewSyncer so the replay lands AFTER every pending
    // structure rebuild — otherwise our overpaint gets applied to objects
    // that are then immediately replaced.
    this.host.viewSyncer.sync('replayHighlight', async () => {
      this.host.logger.debug('[molstar-picker] replay firing now (synced)');
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

  /** Applies distinct base colors to loaded ligand structures so the user
   *  can tell comparison poses apart from the current row. Runs
   *  unconditionally when 2+ poses are loaded; never requires a selection. */
  async applyBaseColors(): Promise<void> {
    const plugin = this.host.getPlugin();
    if (!plugin) return;

    const loadedLigands = this._collectLoadedLigands();
    if (loadedLigands.length < 2) return;

    const allComponents = this._getAllComponents();
    if (!allComponents) return;

    // Leave the currently-selected row in default colors; tint everyone
    // else's carbons grey so heteroatom colors still read clearly.
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

  /** Applies overpaint for ALL loaded ligands. Atom indices are resolved per
   *  row from (in order): the live event (`liveEvent`), the persistent
   *  selection cache, then the picker providers on the molecule column.
   *
   *  NOTE: no `_highlightInProgress` guard — this method is the single
   *  source of truth for highlight state and must always run to completion,
   *  even when a previous overpaint is still in flight. */
  async highlightAllLigandAtoms(
    liveEvent?: {rowIdx: number, atoms: number[], mapping3D: AtomMapping3D | null} | null,
  ): Promise<void> {
    const plugin = this.host.getPlugin();
    if (!plugin) return;
    this.host.logger.debug('[molstar-picker] highlightAllLigandAtoms called');

    const loadedLigands = this._collectLoadedLigands();
    this.host.logger.debug('[molstar-picker] loaded ligands:', () =>
      loadedLigands.map((l) => ({row: l.rowIdx, hasRefs: !!l.structureRefs})));

    const df = this.host.getDataFrame();
    const currentRowIdx = df?.currentRowIdx ?? -1;
    const molCol = df?.columns.toList().find(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE);
    const dfId = df?.id ?? '';
    const dfName = df?.name ?? '';
    const colName = molCol?.name ?? '';
    const cache = getSelectionCache();

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

      // The first ref from `addLigandOnStage` is the ligand's data root;
      // we scope overpaint to its subtree so we don't paint the protein.
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
      // Paint every component but use the loci getter to filter the
      // structure down to atoms belonging to this ligand's data subtree.
      // This dodges the stale-Structure-identity issue while still keeping
      // per-molecule isolation.
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

  // -- Single-ligand highlight (external API) --------------------------------

  /** Highlights specific atoms in the 3D Molstar view using overpaint.
   *  @param atomIndices 0-based atom indices from the 2D picker.
   *  @param mapping3D Optional pre-computed 2D→3D mapping.
   *  @param precomputedSerials When `true`, `atomIndices` are already
   *    1-based Molstar serials and the 2D→3D mapping step is skipped. */
  async highlightLigandAtoms(
    atomIndices: number[],
    mapping3D?: AtomMapping3D | null,
    precomputedSerials?: boolean,
  ): Promise<void> {
    const plugin = this.host.getPlugin();
    if (!plugin) return;
    if (this._highlightInProgress) return;
    this._highlightInProgress = true;

    try {
      const structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures || structures.length === 0) return;

      const allComponents = structures.flatMap((s: any) => s.components ?? []);
      if (allComponents.length === 0) return;

      await clearStructureOverpaint(plugin, allComponents);

      if (atomIndices.length === 0) return;

      let mol3DSerials: number[];
      if (precomputedSerials)
        mol3DSerials = atomIndices;
      else {
        let structure: Structure | undefined;
        const currentLigand = this.host.getLigands()?.current;
        if (currentLigand?.structureRefs && currentLigand.structureRefs.length >= 4) {
          const cell = plugin.state.data.cells.get(currentLigand.structureRefs[3]);
          if (cell?.obj?.data) structure = cell.obj.data;
        }
        if (!structure && structures.length > 0)
          structure = structures[0]?.cell?.obj?.data;

        mol3DSerials = computeSerials(atomIndices, mapping3D, structure);
      }
      if (mol3DSerials.length === 0) return;

      this.host.logger.debug(`[molstar-picker] overpaint serials [${mol3DSerials.slice(0, 10)}]`);

      const serialSet = mol3DSerials;
      await setStructureOverpaint(
        plugin, allComponents, Color(0xFFFF00),
        async (structureData: Structure) => {
          const atomSet = MolScriptBuilder.set(...serialSet);
          const query = MolScriptBuilder.struct.generator.atomGroups({
            'atom-test': MolScriptBuilder.core.set.has([
              atomSet, MolScriptBuilder.ammp('id'),
            ]),
          });
          const sel = Script.getStructureSelection(query, structureData);
          return StructureSelection.toLociWithSourceUnits(sel);
        },
      );
      this.host.logger.debug('[molstar-picker] overpaint applied');
    } catch (err: unknown) {
      this.host.logger.error(
        `highlightLigandAtoms failed: ${err instanceof Error ? err.message : String(err)}`);
    } finally {
      this._highlightInProgress = false;
    }
  }

  /** Clears any atom-level overpaint applied by {@link highlightLigandAtoms}. */
  async clearLigandAtomHighlight(): Promise<void> {
    const plugin = this.host.getPlugin();
    if (!plugin) return;
    try {
      const structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures || structures.length === 0) return;
      const allComponents = structures.flatMap((s: any) => s.components ?? []);
      if (allComponents.length === 0) return;
      await clearStructureOverpaint(plugin, allComponents);
    } catch {
      /* viewer may be disposed */
    }
  }

  // -- Internals -------------------------------------------------------------

  /** Collects loaded ligands (current + hovered + selected) deduplicated by
   *  `rowIdx`. Order matters: `current` wins over `hovered` wins over
   *  `selected` when the same row appears in more than one slot. */
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

  /** Resolves a ligand's Molstar `Structure` via its 4th structureRef, or
   *  returns `undefined` when the state cell isn't available yet. */
  private _getStructureForLigand(ligand: LigandInfo): Structure | undefined {
    const plugin = this.host.getPlugin();
    if (!plugin) return undefined;
    if (!ligand.structureRefs || ligand.structureRefs.length < 4) return undefined;
    const cell = plugin.state.data.cells.get(ligand.structureRefs[3]);
    return cell?.obj?.data;
  }

  /** Returns all Molstar structure components across all loaded structures,
   *  or `null` when nothing is loaded. */
  private _getAllComponents(): any[] | null {
    const plugin = this.host.getPlugin();
    if (!plugin) return null;
    const structures = plugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) return null;
    const comps = structures.flatMap((s: any) => s.components ?? []);
    return comps.length > 0 ? comps : null;
  }

  /** Walks up Molstar's state cell tree from the cell owning `structureData`
   *  to determine whether it descends from the given `dataRef` (up to 6
   *  levels — the depth of `addLigandOnStage`'s transform chain). Errors
   *  during the walk fail-open (`true`) to preserve the original overpaint
   *  behavior. */
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

  /** Reads highlighted atom indices for `rowIdx` from the atom-picker
   *  provider on `molColName` (or the first molecule column that carries
   *  one, as a fallback). Returns `[]` when no picker is set up. */
  private _getHighlightedAtomsFromColumn(rowIdx: number, molColName?: string): number[] {
    const df = this.host.getDataFrame();
    if (!df) return [];

    let molCol: DG.Column | null = null;
    if (molColName)
      molCol = df.col(molColName);
    if (!molCol) {
      molCol = df.columns.toList().find(
        (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE &&
          ((c.temp as Record<string, unknown>)?.['chem-substruct-providers'] as unknown[] ?? [])
            .some((p: any) => p.__atomPicker)) ?? null;
    }
    if (!molCol) return [];

    const providers = ((molCol.temp as Record<string, unknown>)?.['chem-substruct-providers'] ?? []) as
      Array<{__atomPicker?: boolean; __rowIdx?: number; __atoms?: Set<number>}>;
    const picker = providers.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    if (!picker?.__atoms || picker.__atoms.size === 0) return [];
    return [...picker.__atoms];
  }
}
