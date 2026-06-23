import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as chemSearches from '../chem-searches';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import $ from 'cash-dom';
import {Fingerprint} from '../utils/chem-common';
import {renderMolecule} from '../rendering/render-molecule';
import {ChemSearchBaseViewer, SIMILARITY} from './chem-search-base-viewer';
import {malformedDataWarning} from '../utils/malformed-data-utils';
import '../../css/chem.css';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Subscription} from 'rxjs';

export class ChemSimilarityViewer extends ChemSearchBaseViewer {
  followCurrentRow: boolean;
  sketchButton: HTMLElement;
  sketchedMolecule: string = '';
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  cutoff: number;
  searchAsYouSketch: boolean;
  sketcherDebounceMs: number;
  targetMoleculeIdx: number = 0;
  /** Bitset of all rows similar to the reference (passing the cutoff) — drives select & filter. */
  simBitset: DG.BitSet | null = null;
  /** Cached per-row similarity scores, reused for cheap re-slicing when only `limit` changes. */
  simScores: Float32Array | null = null;
  /** Signature of the inputs the bitset was computed from; recompute only when it changes. */
  computeKey: string = '';
  /** The bitset/limit the display columns were last built from — skip the rebuild when unchanged. */
  displayBitset: DG.BitSet | null = null;
  displayLimit: number = -1;
  /** Whether the collaborative filter (filter the table down to similar molecules) is on. */
  applyFilter: boolean = false;
  filterIcon!: HTMLElement;

  override get maxLimit(): number { return 200; }

  get targetMolecule(): string {
    return this.isEditedFromSketcher ?
      this.sketchedMolecule :
      this.moleculeColumn?.get(this.targetMoleculeIdx);
  }

  constructor() {
    super(SIMILARITY);
    this.cutoff = this.float('cutoff', 0.01, {min: 0, max: 1});
    this.followCurrentRow = this.bool('followCurrentRow', true,
      {description: 'Re-compute similarity search when changing current row'});
    this.searchAsYouSketch = this.bool('searchAsYouSketch', true, {
      description: 'Re-run similarity search dynamically while drawing in the sketcher (debounced)',
    });
    this.sketcherDebounceMs = this.int('sketcherDebounceMs', 100, {
      min: 50,
      max: 2000,
      description: 'Debounce delay (ms) for search-as-you-sketch updates',
    });
    this.sketchButton = ui.icons.edit(() => {
      const sketcher = new grok.chem.Sketcher();
      const savedMolecule = this.targetMolecule;
      const savedIsEditedFromSketcher = this.isEditedFromSketcher;
      const savedSketchedMolecule = this.sketchedMolecule;
      sketcher.setMolecule(this.targetMolecule);

      const applySketcherState = (isEdited: boolean, mol: string) => {
        this.isEditedFromSketcher = isEdited;
        this.sketchedMolecule = mol;
        this.gridSelect = false;
        this.render();
      };

      let liveSub: Subscription | null = null;
      if (this.searchAsYouSketch) {
        liveSub = DG.debounce(sketcher.onChanged, this.sketcherDebounceMs).subscribe(() => {
          const mol = sketcher.getMolFile();
          if (DG.chem.Sketcher.isEmptyMolfile(mol))
            return;
          applySketcherState(true, mol);
        });
      }

      const dialog = ui.dialog()
        .add(sketcher.root)
        .onOK(() => {
          const editedMolecule = sketcher.getMolFile();
          if (DG.chem.Sketcher.isEmptyMolfile(editedMolecule)) {
            grok.shell.error(`Empty molecule cannot be used for similarity search`);
            this.isEditedFromSketcher = true;
            this.sketchedMolecule = savedMolecule;
          } else
            applySketcherState(true, editedMolecule);
        })
        .onCancel(() => applySketcherState(savedIsEditedFromSketcher, savedSketchedMolecule));
      dialog.show();
      dialog.onClose.subscribe(() => liveSub?.unsubscribe());
    }, 'Edit');
    this.sketchButton.classList.add('chem-similarity-search-edit');
    this.sketchButton.classList.add('chem-mol-view-icon');
    const selectIcon = ui.iconSvg('select-all', () => this.selectSimilar(), 'Select all similar');
    this.filterIcon = ui.icons.filter(() => this.toggleFilter(), 'Filter table to similar set');
    this.metricsDiv!.appendChild(ui.divH([selectIcon, this.filterIcon],
      'chem-similarity-action-group'));
    this.updateMetricsLink(this, {});
  }

  init(): void {
    this.isEditedFromSketcher = false;
    this.followCurrentRow = true;
    this.initialized = true;
  }

  async onTableAttached(): Promise<void> {
    await super.onTableAttached();
    if (!this.dataFrame)
      return;
    // collaborative filter: when enabled, narrow the table's filter to the similar set
    this.subs.push(this.dataFrame.onRowsFiltering.subscribe(() => {
      if (this.applyFilter && this.simBitset && this.simBitset.length === this.dataFrame.filter.length)
        this.dataFrame.filter.and(this.simBitset, false);
    }));
  }

  detach(): void {
    const wasFiltering = this.applyFilter;
    super.detach();
    if (wasFiltering)
      this.dataFrame?.rows.requestFilter();
  }

  /** Selects every row passing the cutoff — the same set the filter narrows to. One-shot, like any
   * platform selection: press Esc to clear, Ctrl/Shift-click rows to refine. */
  selectSimilar(): void {
    if (this.dataFrame && this.simBitset)
      this.dataFrame.selection.copyFrom(this.simBitset);
  }

  /** Toggles the collaborative filter on/off and re-runs filtering. */
  toggleFilter(): void {
    this.applyFilter = !this.applyFilter;
    this.filterIcon.classList.toggle('active', this.applyFilter);
    this.dataFrame?.rows.requestFilter();
  }

  isReferenceMolecule(idx: number): boolean {
    return idx === this.targetMoleculeIdx && !this.isEditedFromSketcher;
  }

  async renderInternal(computeData: boolean): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.moleculeColumn && this.dataFrame) {
      if (this.moleculeColumn.type !== DG.TYPE.STRING) {
        this.closeWithError('Incorrect target column type');
        return;
      }
      let progressBar: DG.TaskBarProgressIndicator | null = null;
      this.curIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;
      if (computeData && (!this.gridSelect && this.followCurrentRow || this.isEditedFromSketcher)) {
        this.isComputing = true;
        this.error = '';
        this.root.classList.remove(`chem-malformed-molecule-error`);
        this.targetMoleculeIdx = this.dataFrame.currentRowIdx === -1 ? 0 : this.dataFrame.currentRowIdx;
        if (DG.chem.Sketcher.isEmptyMolfile(this.targetMolecule)) {
          this.closeWithError(`Empty molecule cannot be used for similarity search`);
          return;
        }
        try {
          const moleculeColumn = this.moleculeColumn!;
          const computeKey = `${moleculeColumn.name}|${moleculeColumn.version}|` +
            `${this.targetMolecule}|${this.distanceMetric}|${this.fingerprint}|${this.cutoff}`;
          // Recompute fingerprints/bitset only when the data, reference molecule, metric,
          // fingerprint or threshold change — a pure `limit` change just re-slices the cache.
          if (computeKey !== this.computeKey) {
            // claim the key before the await so a concurrent render doesn't start a duplicate search
            this.computeKey = computeKey;
            progressBar = DG.TaskBarProgressIndicator.create(`Similarity search running...`);
            const res = await chemSimilarityBitset(moleculeColumn, this.targetMolecule,
              this.distanceMetric as BitArrayMetrics, this.cutoff, this.fingerprint as Fingerprint);
            if (!res) {
              this.closeWithError(`Malformed molecule cannot be used for similarity search`, progressBar);
              return;
            }
            this.simScores = res.scores;
            this.simBitset = res.bitset;
            if (this.applyFilter) // re-apply the collaborative filter against the new bitset
              this.dataFrame.rows.requestFilter();
          }
        } catch (e: unknown) {
          grok.shell.error(e instanceof Error ? e.message : String(e));
          return;
        } finally {
          progressBar?.close();
        }
      } else if (this.gridSelect)
        this.gridSelect = false;
      if (this.error) {
        this.closeWithError(this.error, progressBar);
        return;
      }
      // Rebuild the top-`limit` display only when the similar set or `limit` changed; a pure redraw
      // (selection, resize) reuses the cached columns, and a `limit` change re-slices without refingerprinting.
      if (this.simBitset && this.simScores && this.moleculeColumn &&
        (this.simBitset !== this.displayBitset || this.limit !== this.displayLimit)) {
        const df = similarityResultDf(this.moleculeColumn, this.simScores, this.simBitset, this.limit);
        this.molCol = df.getCol('smiles');
        this.idxs = df.getCol('indexes');
        this.scores = df.getCol('score');
        this.displayBitset = this.simBitset;
        this.displayLimit = this.limit;
      }
      this.clearResults();
      const panel = [];
      const grids = [];
      let cnt = 0; let cnt2 = 0;
      panel[cnt++] = this.metricsDiv;
      if (this.molCol && this.idxs && this.scores) {
        if (this.isEditedFromSketcher) {
          const label = this.sketchButton;
          const grid = ui.div([
            renderMolecule(
              this.targetMolecule, { width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height }),
            label], { style: { position: 'relative' } });
          let divClass = 'd4-flex-col';
          divClass += ' d4-current';
          grid.style.boxShadow = '0px 0px 1px var(--grey-6)';
          $(grid).addClass(divClass);
          grids[cnt2++] = grid;
        }
        for (let i = 0; i < this.molCol.length; ++i) {
          const idx = this.idxs.get(i);
          const similarity = this.scores.get(i).toPrecision(2);
          const refMolecule = this.isReferenceMolecule(idx);
          const label = refMolecule ? this.sketchButton : ui.div();
          const molProps = this.createMoleculePropertiesDiv(idx, refMolecule, similarity);
          const grid = ui.div([
            renderMolecule(
              this.molCol?.get(i), { width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height }),
            label,
            molProps], { style: { position: 'relative' } });
          let divClass = 'd4-flex-col';
          if (idx === this.curIdx) {
            divClass += ' d4-current';
            grid.style.backgroundColor = '#ddffd9';
          }
          if (idx === this.targetMoleculeIdx && !this.isEditedFromSketcher) {
            divClass += ' d4-current';
            grid.style.boxShadow = '0px 0px 1px var(--grey-6)';
          }
          if (this.dataFrame?.selection.get(idx)) {
            divClass += ' d4-selected';
            if (divClass === 'd4-flex-col d4-selected')
              grid.style.backgroundColor = '#f8f8df';
            else
              grid.style.backgroundColor = '#d3f8bd';
          }
          $(grid).addClass(divClass);
          grid.addEventListener('click', (event: MouseEvent) => {
            if (this.dataFrame && this.idxs) {
              if (event.shiftKey || event.altKey)
                this.dataFrame.selection.set(idx, true);
              else if (event.metaKey) {
                const selected = this.dataFrame.selection;
                this.dataFrame.selection.set(idx, !selected.get(idx));
              } else {
                this.dataFrame.currentRowIdx = idx;
                this.gridSelect = true;
              }
            }
          });
          grids[cnt2++] = grid;
        }
      }
      panel[cnt++] = ui.divH(grids, 'chem-viewer-grid');
      this.root.appendChild(ui.panel([ui.divV(panel)]));
      progressBar?.close();
    }
  }
}

/** Expensive step: computes fingerprints, the reference fingerprint and per-row similarity,
 * returning the raw scores and a bitset of rows passing `minScore`. The fingerprint column is
 * cached by {@link chemSearches.chemGetFingerprints}, so re-running on a threshold change is cheap. */
export async function chemSimilarityBitset(
  smiles: DG.Column,
  molecule: string,
  metricName: BitArrayMetrics,
  minScore: number,
  fingerprint: Fingerprint,
): Promise<{scores: Float32Array, bitset: DG.BitSet} | null> {
  const targetFingerprint = chemSearches.chemGetFingerprint(molecule, fingerprint, () => {return null;});
  if (!targetFingerprint)
    return null; //returning null in case target molecule is malformed
  const fingerprintCol = await chemSearches.chemGetFingerprints(smiles, fingerprint, false);
  malformedDataWarning(fingerprintCol, smiles);

  const fpSim = similarityMetric[metricName];
  const scores = new Float32Array(fingerprintCol.length);
  for (let row = 0; row < fingerprintCol.length; row++) {
    const fp = fingerprintCol[row];
    scores[row] = (!fp || fp.allFalse) ? -1 : fpSim(targetFingerprint, fp);
  }
  const bitset = DG.BitSet.create(scores.length, (i) => scores[i] >= minScore);
  return {scores, bitset};
}

/** Cheap step: builds the display dataframe (top-`limit` rows by similarity) from the scores and
 * bitset produced by {@link chemSimilarityBitset}. Re-run on its own when only `limit` changes. */
export function similarityResultDf(
  smiles: DG.Column,
  scores: Float32Array,
  bitset: DG.BitSet,
  limit: number,
): DG.DataFrame {
  const indexes = Array.from(bitset.getSelectedIndexes());
  indexes.sort((i1, i2) => scores[i2] - scores[i1]);
  limit = Math.min(indexes.length, limit);
  const molsList = new Array<string>(limit);
  const scoresList = new Array<number>(limit);
  const molsIdxs = new Array<number>(limit);
  for (let n = 0; n < limit; n++) {
    const idx = indexes[n];
    molsIdxs[n] = idx;
    molsList[n] = smiles.get(idx);
    scoresList[n] = scores[idx];
  }
  const mols = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'smiles', molsList);
  mols.semType = DG.SEMTYPE.MOLECULE;
  const scoresCol = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'score', scoresList);
  const newIndexes = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'indexes', molsIdxs);
  return DG.DataFrame.fromColumns([mols, scoresCol, newIndexes]);
}

/** Convenience wrapper: runs both steps and returns the top-`limit` display dataframe,
 * or null when the reference molecule is malformed. */
export async function chemSimilaritySearch(
  smiles: DG.Column,
  molecule: string,
  metricName: BitArrayMetrics,
  limit: number,
  minScore: number,
  fingerprint: Fingerprint,
) : Promise<DG.DataFrame | null> {
  const res = await chemSimilarityBitset(smiles, molecule, metricName, minScore, fingerprint);
  return res ? similarityResultDf(smiles, res.scores, res.bitset, limit) : null;
}
