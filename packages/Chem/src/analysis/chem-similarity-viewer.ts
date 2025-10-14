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

export class ChemSimilarityViewer extends ChemSearchBaseViewer {
  followCurrentRow: boolean;
  sketchButton: HTMLElement;
  sketchedMolecule: string = '';
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  cutoff: number;
  targetMoleculeIdx: number = 0;

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
    this.sketchButton = ui.icons.edit(() => {
      const sketcher = new grok.chem.Sketcher();
      const savedMolecule = this.targetMolecule;
      sketcher.setMolecule(this.targetMolecule);
      ui.dialog()
        .add(sketcher.root)
        .onOK(() => {
          this.isEditedFromSketcher = true;
          const editedMolecule = sketcher.getMolFile();
          if (DG.chem.Sketcher.isEmptyMolfile(editedMolecule)) {
            grok.shell.error(`Empty molecule cannot be used for similarity search`);
            this.sketchedMolecule = savedMolecule;
          } else {
            this.sketchedMolecule = sketcher.getMolFile();
            this.gridSelect = false;
            this.render();
          }
        })
        .show();
    });
    this.sketchButton.classList.add('chem-similarity-search-edit');
    this.sketchButton.classList.add('chem-mol-view-icon');
    this.updateMetricsLink(this, {});
  }

  init(): void {
    this.isEditedFromSketcher = false;
    this.followCurrentRow = true;
    this.initialized = true;
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
      this.curIdx = this.dataFrame.currentRowIdx == -1 ? 0 : this.dataFrame.currentRowIdx;
      if (computeData && (!this.gridSelect && this.followCurrentRow || this.isEditedFromSketcher)) {
        progressBar = DG.TaskBarProgressIndicator.create(`Similarity search running...`);
        this.isComputing = true;
        this.error = '';
        this.root.classList.remove(`chem-malformed-molecule-error`);
        this.targetMoleculeIdx = this.dataFrame.currentRowIdx == -1 ? 0 : this.dataFrame.currentRowIdx;
        if (DG.chem.Sketcher.isEmptyMolfile(this.targetMolecule)) {
          this.closeWithError(`Empty molecule cannot be used for similarity search`, progressBar);
          return;
        }
        try {
          const rowSourceIdxs = this.getRowSourceIndexes();
          rowSourceIdxs.set(this.targetMoleculeIdx, true);
          const df = await chemSimilaritySearch(this.dataFrame!, this.moleculeColumn!,
            this.targetMolecule, this.distanceMetric as BitArrayMetrics, this.limit, this.cutoff,
            this.fingerprint as Fingerprint, rowSourceIdxs);
          if (!df) {
            this.closeWithError(`Malformed molecule cannot be used for similarity search`, progressBar);
            return;
          }
          this.molCol = df.getCol('smiles');
          this.idxs = df.getCol('indexes');
          this.scores = df.getCol('score');
        } catch (e: any) {
          grok.shell.error(e.message);
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
          if (idx == this.curIdx) {
            divClass += ' d4-current';
            grid.style.backgroundColor = '#ddffd9';
          }
          if (idx == this.targetMoleculeIdx && !this.isEditedFromSketcher) {
            divClass += ' d4-current';
            grid.style.boxShadow = '0px 0px 1px var(--grey-6)';
          }
          if (this.dataFrame?.selection.get(idx)) {
            divClass += ' d4-selected';
            if (divClass == 'd4-flex-col d4-selected')
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

export async function chemSimilaritySearch(
  table: DG.DataFrame,
  smiles: DG.Column,
  molecule: string,
  metricName: BitArrayMetrics,
  limit: number,
  minScore: number,
  fingerprint: Fingerprint,
  rowSourceIndexes: DG.BitSet,
) : Promise<DG.DataFrame | null> {
  const targetFingerprint = chemSearches.chemGetFingerprint(molecule, fingerprint, () => {return null;});
  if (!targetFingerprint)
    return null; //returning null in case target molecule is malformed
  const fingerprintCol = await chemSearches.chemGetFingerprints(smiles, fingerprint, false);
  malformedDataWarning(fingerprintCol, smiles);
  const distances: number[] = [];

  const fpSim = similarityMetric[metricName];
  for (let row = 0; row < fingerprintCol.length; row++) {
    const fp = fingerprintCol[row];
    distances[row] = (!fp || fp!.allFalse) ? 100.0 : fpSim(targetFingerprint, fp!);
  }

  function range(end: number) {
    return Array(end).fill(0).map((_, idx) => idx);
  }

  function compare(i1: number, i2: number) {
    if (distances[i1] > distances[i2])
      return -1;

    if (distances[i1] < distances[i2])
      return 1;

    return 0;
  }

  const indexes = range(table.rowCount)
    .filter((idx) => fingerprintCol[idx] && !fingerprintCol[idx]!.allFalse && rowSourceIndexes.get(idx))
    .sort(compare);
  const molsList = [];
  const scoresList = [];
  const molsIdxs = [];
  limit = Math.min(indexes.length, limit);
  for (let n = 0; n < limit; n++) {
    const idx = indexes[n];
    const score = distances[idx];
    if (score < minScore)
      break;

    molsIdxs[n] = idx;
    molsList[n] = smiles.get(idx);
    scoresList[n] = score;
  }
  const mols = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'smiles', molsList);
  mols.semType = DG.SEMTYPE.MOLECULE;
  const scores = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'score', scoresList);
  const newIndexes = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'indexes', molsIdxs);
  return DG.DataFrame.fromColumns([mols, scores, newIndexes]);
}
