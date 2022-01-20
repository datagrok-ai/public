import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as metric from './chem-common'
import { Property, SIMILARITY_METRIC } from 'datagrok-api/dg';
import { randomInt } from '@datagrok-libraries/utils/src/operations'
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/analysis';
import { chemGetMorganFingerprint, chemGetMorganFingerprints } from './chem-searches';
import $ from 'cash-dom'

export class DiversitySearch extends DG.JsViewer {
  private moleculeColumnName: string;
  private initialized: boolean;
  private distanceMetric: string;
  private limit: number;
  private renderMolIds: number[];
  private fpSim;

  constructor() {
    super();

    this.moleculeColumnName = this.string('moleculeColumnName');
    this.limit = this.int('limit', 10);
    this.distanceMetric = this.string('distanceMetric', 'tanimoto', {choices: Object.values(SIMILARITY_METRIC)});
    this.initialized = false;
    this.renderMolIds = [];
    const metrics: {[Key: string]: any} = {
      'tanimoto': metric.tanimotoSimilarity,
      'dice': metric.diceSimilarity,
      'cosine': metric.cosineSimilarity,
      'sokal': metric.sokalSimilarity,
      'kulczynski': metric.kulczynskiSimilarity,
      'mc-connaughey': metric.mcConnaugheySimilarity,
      'asymmetric': metric.asymmetricSimilarity,
      'braun-blanquet': metric.braunBlanquetSimilarity,
      'russel': metric.russelSimilarity,
      'rogot-goldberg': metric.rogotGoldbergSimilarity,
    }
    this.fpSim = metrics[this.distanceMetric]
  }

  init(): void {
    this.initialized = true;
  }

  async onTableAttached(): Promise<void> {
    this.init();

    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.onRowsRemoved, 50).subscribe(async (_) => await this.render()));
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe(async (_) => await this.render(false)));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe(async (_) => await this.render(false)));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe(async (_) => await this.render(false)));
      for (const col of this.dataFrame.columns) {
        if (col.semType == 'Molecule') {
          this.moleculeColumnName = col.name;
        }
      }
    }

    await this.render();
  }

  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onPropertyChanged(property: Property): void {
    super.onPropertyChanged(property);
    if (!this.initialized)
      return;
    this.render();
  }

  async render(computeData = true): Promise<void> {
    if (!this.initialized)
      return;

    if (this.dataFrame) {
      if (computeData) {
        this.renderMolIds = await chemDiversitySearch(this.dataFrame.getCol(this.moleculeColumnName), this.fpSim, this.limit);
      }

      if (this.root.hasChildNodes())
        this.root.removeChild(this.root.childNodes[0]);
      
        const panel = [];
        const grids = []; 
        let cnt = 0, cnt2 = 0;
        panel[cnt++] = ui.h1('Diverse structures');
        for (let i = 0; i < this.limit; ++i) {
          const mol = grok.chem.svgMol(this.dataFrame.getCol(this.moleculeColumnName).get(this.renderMolIds[i]), 200, 100);
          let grid = ui.div([mol], {style: {width: '200px', height: '100px', margin: '5px'}});
          let divClass = 'd4-flex-col';
  
          if (this.renderMolIds[i] == this.dataFrame.currentRowIdx) {
            divClass += ' d4-current';
            grid.style.backgroundColor = '#ddffd9';
          } 
          if (this.dataFrame.selection.get(this.renderMolIds[i])) {
            divClass += ' d4-selected';
            if (divClass == 'd4-flex-col d4-selected')
              grid.style.backgroundColor = '#f8f8df';
            else 
              grid.style.backgroundColor = '#d3f8bd';
          }
  
          $(grid).addClass(divClass);
          grid.addEventListener('click', (event: MouseEvent) => {
            if (this.dataFrame) {
              if (event.shiftKey || event.altKey) {
                this.dataFrame.selection.set(this.renderMolIds[i], true);
              } else if (event.metaKey) {
                let selected = this.dataFrame.selection;
                this.dataFrame.selection.set(this.renderMolIds[i], selected.get(this.renderMolIds[i]) != true);
              } else {
                this.dataFrame.currentRowIdx = this.renderMolIds[i];
              }
            }
          });
          grids[cnt2++] = grid;
        }
        const gridsDiv = ui.div(grids, {classes: 'd4-flex-wrap'});
        panel[cnt++] = gridsDiv;
        this.root.appendChild(ui.div(panel, {style: {margin: '5px'}}));
    }
  }
}

export async function chemDiversitySearch(smiles: DG.Column, fpSim: (a: BitArray, b: BitArray) => number, limit: number){
  limit = Math.min(limit, smiles.length);
  let fingerprintCol = await chemGetMorganFingerprints(smiles);
  let indexes: number[] = [];

  for (let i = 0; i < fingerprintCol.length; ++i) {
    if (fingerprintCol[i] != null) {
      indexes.push(i);
    }
  }

  let diverseIndexes = getDiverseSubset(indexes.length, limit,
          (i1, i2) => 1 - fpSim(fingerprintCol[indexes[i1]], fingerprintCol[indexes[i2]]));

  let molIds: number[] = [];
  for (let i = 0; i < limit; i++)
    molIds[i] = indexes[diverseIndexes[i]];

  return molIds;
}