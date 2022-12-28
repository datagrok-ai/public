import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {chemGetFingerprints} from '../chem-searches';
import $ from 'cash-dom';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import {Fingerprint} from '../utils/chem-common';
import {renderMolecule} from '../rendering/render-molecule';
import {ChemSearchBaseViewer} from './chem-search-base-viewer';

export class ChemDiversityViewer extends ChemSearchBaseViewer {
  renderMolIds: number[];
  columnNames = [];

  constructor() {
    super('diversity');
    this.renderMolIds = [];
    this.updateMetricsLink(this.metricsDiv, this, {fontSize: '10px', fontWeight: 'normal', paddingBottom: '15px'});
  }


  async render(computeData = true): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.dataFrame && this.moleculeColumn) {
      const progressBar = DG.TaskBarProgressIndicator.create(`Diversity search running...`);
      if (computeData) {
        const rowsWithoutEmptyValues = rowsWithoutEmptyValuesCount(this.moleculeColumn);
        if (this.limit > rowsWithoutEmptyValues)
          this.limit = rowsWithoutEmptyValues;
        this.renderMolIds =
          await chemDiversitySearch(
            this.moleculeColumn, similarityMetric[this.distanceMetric], this.limit, this.fingerprint as Fingerprint);
      }
      if (this.root.hasChildNodes())
        this.root.removeChild(this.root.childNodes[0]);

      const panel = [];
      const grids = [];
      let cnt = 0; let cnt2 = 0;

      panel[cnt++] = this.metricsDiv;
      for (let i = 0; i < this.limit; ++i) {
        const grid = ui.div([
          renderMolecule(
            this.moleculeColumn!.get(this.renderMolIds[i]),
            {width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height}),
        ], {style: {margin: '5px'}});

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
            if (event.shiftKey || event.altKey)
              this.dataFrame.selection.set(this.renderMolIds[i], true);
            else if (event.metaKey) {
              const selected = this.dataFrame.selection;
              this.dataFrame.selection.set(this.renderMolIds[i], !selected.get(this.renderMolIds[i]));
            } else
              this.dataFrame.currentRowIdx = this.renderMolIds[i];
          }
        });
        grids[cnt2++] = grid;
      }

      panel[cnt++] = ui.div(grids, {classes: 'd4-flex-wrap'});
      this.root.appendChild(ui.div(panel, {style: {margin: '5px'}}));
      progressBar.close();
    }
  }
}

function rowsWithoutEmptyValuesCount(col: DG.Column): number {
  const categories = col.categories;
  const rawData = col.getRawData();
  return rawData.filter((it) => categories[it]).length;
}

export async function chemDiversitySearch(
  moleculeColumn: DG.Column, similarity: (a: BitArray, b: BitArray) => number,
  limit: number, fingerprint: Fingerprint): Promise<number[]> {
  limit = Math.min(limit, moleculeColumn.length);
  const fingerprintArray = await chemGetFingerprints(moleculeColumn, fingerprint);
  const indexes = ArrayUtils.indexesOf(fingerprintArray, (f) => f != null);

  const diverseIndexes = getDiverseSubset(fingerprintArray, indexes.length, limit,
    (i1, i2) => 1 - similarity(fingerprintArray[indexes[i1]], fingerprintArray[indexes[i2]]));

  const molIds: number[] = [];
  for (let i = 0; i < limit; i++)
    molIds[i] = indexes[diverseIndexes[i]];

  return molIds;
}
