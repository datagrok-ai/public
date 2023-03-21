import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {similarityMetric} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {chemGetFingerprints} from '../chem-searches';
import $ from 'cash-dom';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import {defaultMorganFpLength, defaultMorganFpRadius, Fingerprint, rdKitFingerprintToBitArray} from '../utils/chem-common';
import {renderMolecule} from '../rendering/render-molecule';
import {ChemSearchBaseViewer} from './chem-search-base-viewer';
import { malformedDataWarning } from '../utils/malformed-data-utils';
import { getRdKitModule } from '../package';
import { RDMol } from '@datagrok-libraries/chem-meta/src/rdkit-api';

export class ChemDiversityViewer extends ChemSearchBaseViewer {
  renderMolIds: number[];
  columnNames = [];
  tooltipUse: boolean;

  constructor(tooltipUse = false, col?: DG.Column) {
    super('diversity', col);
    this.renderMolIds = [];
    this.updateMetricsLink(this.metricsDiv, this, {fontSize: '10px', fontWeight: 'normal', paddingBottom: '15px'});
    this.tooltipUse = tooltipUse;
  }


  async render(computeData = true): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.dataFrame && this.moleculeColumn) {
      let progressBar: DG.TaskBarProgressIndicator;
      if (!this.tooltipUse)
        progressBar = DG.TaskBarProgressIndicator.create(`Diversity search running...`);

      if (computeData) {
        this.renderMolIds =
          await chemDiversitySearch(
            this.moleculeColumn, similarityMetric[this.distanceMetric], this.limit, this.fingerprint as Fingerprint, this.tooltipUse);
      }
      if (this.root.hasChildNodes())
        this.root.removeChild(this.root.childNodes[0]);

      const panel = [];
      const grids = [];
      let cnt = 0; let cnt2 = 0;

      panel[cnt++] = this.metricsDiv;
      for (let i = 0; i < this.renderMolIds.length; ++i) {
        const molProps = this.createMoleculePropertiesDiv(this.renderMolIds[i]);
        const grid = ui.div([
          renderMolecule(
            this.moleculeColumn!.get(this.renderMolIds[i]),
            {width: this.sizesMap[this.size].width, height: this.sizesMap[this.size].height, popupMenu: !this.tooltipUse}),
          molProps],
        {style: {margin: '5px', padding: '3px', position: 'relative'}},
        );

        let divClass = 'd4-flex-col';

        if (this.renderMolIds[i] == this.dataFrame.currentRowIdx) {
          divClass += ' d4-current';
          grid.style.backgroundColor = '#ddffd9';
        }
        if (!this.tooltipUse && this.dataFrame.selection.get(this.renderMolIds[i])) {
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

      if (this.tooltipUse)
        panel[cnt++] = ui.divText('Most diverse structures');
      panel[cnt++] = ui.div(grids, {classes: 'd4-flex-wrap'});
      this.root.appendChild(ui.div(panel, {style: {margin: '5px'}}));
      if (!this.tooltipUse)
        progressBar!.close();
    }
  }
}

export async function chemDiversitySearch(
  moleculeColumn: DG.Column, similarity: (a: BitArray, b: BitArray) => number,
  limit: number, fingerprint: Fingerprint, tooltipUse: boolean = false): Promise<number[]> {

  let fingerprintArray: BitArray[] = [];
  if (tooltipUse) {
    const size = Math.min(moleculeColumn.length, 1000);
    const randomIndexes = Array.from({length: size}, () => Math.floor(Math.random() * moleculeColumn.length));
    for (let i = 0; i < randomIndexes.length; ++i) {
      let mol: RDMol | null = null;
      try {
        mol = getRdKitModule().get_mol(moleculeColumn.get(randomIndexes[i]));
        let fp = mol.get_morgan_fp_as_uint8array(JSON.stringify({
          radius: defaultMorganFpRadius,
          nBits: defaultMorganFpLength,
        }));
        fingerprintArray[i] = rdKitFingerprintToBitArray(fp);
      } catch (e) {
      } finally {
        mol?.delete();
      }
    }
  } else {
    fingerprintArray = await chemGetFingerprints(moleculeColumn, fingerprint);
  }
  if (!tooltipUse)
    malformedDataWarning(fingerprintArray, moleculeColumn.dataFrame);
  const indexes = ArrayUtils.indexesOf(fingerprintArray, (f) => !f.allFalse);
  limit = Math.min(limit, indexes.length);

  const diverseIndexes = getDiverseSubset(indexes.length, limit,
    (i1, i2) => 1 - similarity(fingerprintArray[indexes[i1]], fingerprintArray[indexes[i2]]));

  const molIds: number[] = [];
  for (let i = 0; i < limit; i++)
    molIds[i] = indexes[diverseIndexes[i]];

  return molIds;
}
