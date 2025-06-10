
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {SequenceSearchBaseViewer} from './sequence-search-base-viewer';
import {getMonomericMols} from '../calculations/monomerLevelMols';
import {adjustGridcolAfterRender, updateDivInnerHTML} from '../utils/ui-utils';
import {Subject} from 'rxjs';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getEncodedSeqSpaceCol} from './sequence-space';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {DistanceMatrixService, dmLinearIndex} from '@datagrok-libraries/ml/src/distance-matrix';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {_package} from '../package';
import {getMonomerSubstitutionMatrix} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

const MAX_SAMPLE_SIZE = 10000;


export class SequenceDiversityViewer extends SequenceSearchBaseViewer {
  diverseColumnLabel: string | null;

  renderMolIds: number[] | null = null;
  columnNames = [];
  computeCompleted = new Subject<boolean>();

  private sampledIndices: number[] | null = null;

  constructor(
    private readonly seqHelper: ISeqHelper,
  ) {
    super('diversity', DG.SEMTYPE.MACROMOLECULE);
    this.diverseColumnLabel = this.string('diverseColumnLabel', null);
  }

  override async renderInt(computeData: boolean): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.dataFrame) {
      if (computeData && this.targetColumn) {
        await this.computeByMM();

        const diverseColumnName: string = this.diverseColumnLabel != null ? this.diverseColumnLabel :
          `diverse (${this.targetColumnName})`;
        const resCol = DG.Column.string(diverseColumnName, this.renderMolIds!.length)
          .init((i) => this.targetColumn?.get(this.renderMolIds![i]));
        resCol.semType = DG.SEMTYPE.MACROMOLECULE;
        this.tags.forEach((tag) => resCol.setTag(tag, this.targetColumn!.getTag(tag)));
        const resDf = DG.DataFrame.fromColumns([resCol]);
        resCol.temp[MmcrTemps.maxMonomerLength] = 4;

        const _ = resDf.onCurrentRowChanged.subscribe((_: any) => {
          this.dataFrame.currentRowIdx = this.renderMolIds![resDf.currentRowIdx];
        });

        const grid = resDf.plot.grid();
        adjustGridcolAfterRender(grid, resCol.name, 450, 30);

        updateDivInnerHTML(this.root, grid.root);
        this.computeCompleted.next(true);
      }
    }
  }

  private async computeByMM() {
    const totalLength = this.targetColumn!.length;
    let workingIndices: number[];

    // Determine if we need to sample the data
    if (this.requiresSampling && totalLength > MAX_SAMPLE_SIZE) {
      workingIndices = this.createRandomSample(totalLength, MAX_SAMPLE_SIZE);
      this.sampledIndices = workingIndices;
      grok.shell.info(`Sampled ${workingIndices.length} sequences from ${totalLength} for diversity calculation.`);
    } else {
      // Use the full dataset
      workingIndices = Array.from({length: totalLength}, (_, i) => i);
      this.sampledIndices = null;
    }

    const distanceFunction = this.distanceMetric as MmDistanceFunctionsNames;
    const params = this.getDistanceFunctionParams();

    const encodedResult = await getEncodedSeqSpaceCol(this.targetColumn!, distanceFunction, params);
    const fullEncodedSequences = encodedResult.seqList;
    const options = encodedResult.options;

    // Extract only the sequences we need for the working set
    const workingEncodedSequences = workingIndices.map((idx) => fullEncodedSequences[idx]);

    const distanceMatrixService = new DistanceMatrixService(true, false);
    const distanceMatrixData = await distanceMatrixService.calc(
      workingEncodedSequences,
      distanceFunction,
      true, // normalize
      options
    );
    distanceMatrixService.terminate();

    const workingLength = workingIndices.length;
    const linearizeFunc = dmLinearIndex(workingLength);

    const diverseIndicesInWorkingSet = getDiverseSubset(
      workingLength,
      Math.min(workingLength, this.limit),
      (i1: number, i2: number) => {
        return this.targetColumn!.isNone(workingIndices[i1]) || this.targetColumn!.isNone(workingIndices[i2]) ? 0 :
          distanceMatrixData[linearizeFunc(i1, i2)];
      }
    );

    // Map back to original indices
    this.renderMolIds = diverseIndicesInWorkingSet.map((workingIndex) => workingIndices[workingIndex]);
  }

  private createRandomSample(totalLength: number, sampleSize: number): number[] {
    const validIndices: number[] = [];
    for (let i = 0; i < totalLength; i++) {
      if (!this.targetColumn!.isNone(i))
        validIndices.push(i);
    }

    if (validIndices.length <= sampleSize)
      return validIndices;

    for (let i = validIndices.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      const temp = validIndices[i];
      validIndices[i] = validIndices[j];
      validIndices[j] = temp;
    }

    // Return first sampleSize elements, sorted for better cache performance
    const result = validIndices.slice(0, sampleSize);
    result.sort((a, b) => a - b);
    return result;
  }

  // // Helper method to get information about sampling (useful for debugging/info)
  // getSamplingInfo(): {isSampled: boolean, originalSize?: number, sampleSize?: number, sampledIndices?: number[]} {
  //   if (this.sampledIndices) {
  //     return {
  //       isSampled: true,
  //       originalSize: this.targetColumn?.length,
  //       sampleSize: this.sampledIndices.length,
  //       sampledIndices: this.sampledIndices
  //     };
  //   }
  //   return {isSampled: false};
  // }
}
