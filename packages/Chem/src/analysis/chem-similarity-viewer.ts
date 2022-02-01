import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Property} from 'datagrok-api/dg';
import * as chemSearches from '../chem-searches';
import {similarityMetric} from '@datagrok-libraries/utils/src/similarity-metrics';
import $ from 'cash-dom';
import {Fingerprint} from "../utils/chem-common";
import {chem} from "datagrok-api/grok";
import Sketcher = chem.Sketcher;
import {renderMolecule} from "../rendering/render-molecule";

export class ChemSimilarityViewer extends DG.JsViewer {
  moleculeColumn?: DG.Column;
  isEditedFromSketcher: boolean = false;
  hotSearch: boolean;
  sketchButton: HTMLButtonElement;
  sketchedMolecule: string = "";
  distanceMetric: string;
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  limit: number;
  minScore: number;
  fingerprint: string;
  gridSelect: boolean = false;
  targetMoleculeIdx: number = 0;

  get targetMolecule(): string {
    return this.isEditedFromSketcher
      ? this.sketchedMolecule
      : this.moleculeColumn?.get(this.targetMoleculeIdx);
  }

  constructor() {
    super();
    this.limit = this.int('limit', 10);
    this.minScore = this.float('minScore', 0.1);
    this.distanceMetric = this.string('distanceMetric', 'Tanimoto', {choices: Object.keys(similarityMetric)});
    this.fingerprint = this.string('fingerprint', 'Morgan', {choices: ['Morgan', 'RDKit', 'Pattern']});
    this.hotSearch = this.bool('hotSearch', true);
    this.sketchButton = ui.button('Sketch', () => {
      const sketcher = new Sketcher();
      sketcher.setMolecule(this.targetMolecule);
      ui.dialog()
        .add(sketcher.root)
        .onOK(() => {
            this.isEditedFromSketcher = true;
            this.sketchedMolecule = sketcher.getMolFile();
            this.render();
          })
        .show();
    });
    this.sketchButton.id = 'reference';
  }

  onTableAttached() {
    this.isEditedFromSketcher = false;
    this.hotSearch = true;
    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.onRowsRemoved, 50).subscribe(async (_) => await this.render()));
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe(async (_) => await this.render()));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe(async (_) => await this.render(false)));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe(async (_) => await this.render(false)));
      grok.data.detectSemanticTypes(this.dataFrame).then(() => {
        this.moleculeColumn = this.dataFrame!.columns.bySemType(DG.SEMTYPE.MOLECULE);
        this.render();
      });
    }
  }

  onPropertyChanged(property: Property): void {
    super.onPropertyChanged(property);
    this.render();
  }

  async render(computeData = true) {
    if (this.moleculeColumn) {
      if (!this.gridSelect && this.curIdx != this.dataFrame!.currentRowIdx) {
        this.isEditedFromSketcher = false;
      }
      this.curIdx = this.dataFrame!.currentRowIdx;
      if (computeData && !this.gridSelect) {
        this.targetMoleculeIdx = this.dataFrame!.currentRowIdx;
        const df = await chemSimilaritySearch(this.dataFrame!, this.moleculeColumn!,
          this.targetMolecule, this.distanceMetric, this.limit, this.minScore, this.fingerprint as Fingerprint);
        this.molCol = df.getCol('smiles');
        this.idxs = df.getCol('indexes');
        this.scores = df.getCol('score');
      }
      else if (this.gridSelect)
        this.gridSelect = false;
      if (this.root.hasChildNodes())
        this.root.removeChild(this.root.childNodes[0]);
      const panel = [];
      const grids = []; 
      let cnt = 0, cnt2 = 0;
      panel[cnt++] = ui.h1('Reference');
      panel[cnt++] = renderMolecule(this.targetMolecule);
      panel[cnt++] = this.sketchButton;
      if (this.molCol && this.idxs && this.scores) {
        for (let i = 0; i < this.molCol.length; ++i) {
          let grid = ui.div([
            renderMolecule(this.molCol?.get(i)),
            ui.label(`${this.scores.get(i).toPrecision(2)}`)],
            {style: {width: '200px', height: '120px', margin: '5px'}}
          );
          let divClass = 'd4-flex-col';
          if (this.idxs.get(i) == this.curIdx) {
            divClass += ' d4-current';
            grid.style.backgroundColor = '#ddffd9';
          } 
          if (this.dataFrame!.selection.get(this.idxs.get(i))) {
            divClass += ' d4-selected';
            if (divClass == 'd4-flex-col d4-selected')
              grid.style.backgroundColor = '#f8f8df';
            else 
              grid.style.backgroundColor = '#d3f8bd';
          }
          $(grid).addClass(divClass);
          grid.addEventListener('click', (event: MouseEvent) => {
            if (this.dataFrame && this.idxs) {
              if (event.shiftKey || event.altKey) {
                this.dataFrame.selection.set(this.idxs.get(i), true);
              } else if (event.metaKey) {
                let selected = this.dataFrame.selection;
                this.dataFrame.selection.set(this.idxs.get(i), !selected.get(this.idxs.get(i)));
              } else {
                this.dataFrame.currentRowIdx = this.idxs.get(i);
                this.gridSelect = true;
              }
            }
          });
          grids[cnt2++] = grid;
        }
      }
      panel[cnt++] = ui.div(grids, {classes: 'd4-flex-wrap'});
      this.root.appendChild(ui.div(panel, {style: {margin: '5px'}}));
    }
  }
}

export async function chemSimilaritySearch(
  table: DG.DataFrame,
  smiles: DG.Column,
  molecule: string,
  metricName: string,
  limit: number,
  minScore: number,
  fingerprint: Fingerprint
) {
  limit = Math.min(limit, smiles.length);
  const targetFingerprint = chemSearches.chemGetFingerprint(molecule, fingerprint);
  const fingerprintCol = await chemSearches.chemGetFingerprints(smiles, fingerprint);
  const distances: number[] = [];

  let fpSim = similarityMetric[metricName];
  for (let row = 0; row < fingerprintCol.length; row++) {
    const fp = fingerprintCol[row];
    distances[row] = fp == null ? 100.0 : fpSim(targetFingerprint, fp);
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
    .filter((idx) => fingerprintCol[idx] != null)
    .sort(compare);
  const molsList = [];
  const scoresList = [];
  const molsIdxs = [];

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
