import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {SequenceSearchBaseViewer} from './sequence-search-base-viewer';
import {getMonomericMols} from '../calculations/monomerLevelMols';
import {updateDivInnerHTML} from '../utils/ui-utils';
import {Subject} from 'rxjs';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getEncodedSeqSpaceCol} from './sequence-space';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {DistanceMatrixService, dmLinearIndex} from '@datagrok-libraries/ml/src/distance-matrix';

export class SequenceDiversityViewer extends SequenceSearchBaseViewer {
  diverseColumnLabel: string | null; // Use postfix Label to prevent activating table column selection editor

  renderMolIds: number[] | null = null;
  columnNames = [];
  computeCompleted = new Subject<boolean>();

  constructor(
    private readonly seqHelper: ISeqHelper,
  ) {
    super('diversity');
    this.diverseColumnLabel = this.string('diverseColumnLabel', null);
  }

  override async renderInt(computeData: boolean): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.dataFrame) {
      if (computeData && this.moleculeColumn) {
        const sh = this.seqHelper.getSeqHandler(this.moleculeColumn);
        await (sh.isFasta() ? this.computeByMM() : this.computeByChem());

        const diverseColumnName: string = this.diverseColumnLabel != null ? this.diverseColumnLabel :
          `diverse (${this.moleculeColumnName})`;
        const resCol = DG.Column.string(diverseColumnName, this.renderMolIds!.length)
          .init((i) => this.moleculeColumn?.get(this.renderMolIds![i]));
        resCol.semType = DG.SEMTYPE.MACROMOLECULE;
        this.tags.forEach((tag) => resCol.setTag(tag, this.moleculeColumn!.getTag(tag)));
        const resDf = DG.DataFrame.fromColumns([resCol]);
        resDf.onCurrentRowChanged.subscribe(
          (_: any) => { this.dataFrame.currentRowIdx = this.renderMolIds![resDf.currentRowIdx]; });
        updateDivInnerHTML(this.root, resDf.plot.grid().root);
        this.computeCompleted.next(true);
      }
    }
  }

  private async computeByChem() {
    const monomericMols = await getMonomericMols(this.moleculeColumn!, this.seqHelper);
    //need to create df to calculate fingerprints
    const _monomericMolsDf = DG.DataFrame.fromColumns([monomericMols]);
    this.renderMolIds = await grok.functions.call('Chem:callChemDiversitySearch', {
      col: monomericMols,
      metricName: this.distanceMetric,
      limit: this.limit,
      fingerprint: this.fingerprint,
    });
  }

  private async computeByMM() {
    const encodedSequences =
      (await getEncodedSeqSpaceCol(this.moleculeColumn!, MmDistanceFunctionsNames.LEVENSHTEIN)).seqList;
    const distanceMatrixService = new DistanceMatrixService(true, false);
    const distanceMatrixData = await distanceMatrixService.calc(encodedSequences, MmDistanceFunctionsNames.LEVENSHTEIN);
    distanceMatrixService.terminate();
    const len = this.moleculeColumn!.length;
    const linearizeFunc = dmLinearIndex(len);
    this.renderMolIds = getDiverseSubset(len, Math.min(len, this.limit),
      (i1: number, i2: number) => {
        return this.moleculeColumn!.isNone(i1) || this.moleculeColumn!.isNone(i2) ? 0 :
          distanceMatrixData[linearizeFunc(i1, i2)];
      });
  }
}
