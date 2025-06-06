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

const MAX_SAMPLE_SIZE = 10000; // Maximum size for sampling

export class SequenceDiversityViewer extends SequenceSearchBaseViewer {
  diverseColumnLabel: string | null; // Use postfix Label to prevent activating table column selection editor

  renderMolIds: number[] | null = null;
  columnNames = [];
  computeCompleted = new Subject<boolean>();

  // Store the original indices when sampling is used
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
    let workingColumn: DG.Column;

    // Determine if we need to sample the data
    if (this.requiresSampling && totalLength > MAX_SAMPLE_SIZE) {
      // Create a random sample of indices
      workingIndices = this.createRandomSample(totalLength, MAX_SAMPLE_SIZE);
      this.sampledIndices = workingIndices;

      // Create a new column with only the sampled sequences
      workingColumn = DG.Column.string(`${this.targetColumn!.name}_sampled`, workingIndices.length)
        .init((i) => this.targetColumn!.get(workingIndices[i]));
      workingColumn.semType = this.targetColumn!.semType;

      // Copy all relevant tags
      this.tags.forEach((tag) => {
        if (this.targetColumn!.getTag(tag))
          workingColumn.setTag(tag, this.targetColumn!.getTag(tag));
      });

      grok.shell.info(`Sampled ${workingIndices.length} sequences from ${totalLength} for diversity calculation.`);
    } else {
      // Use the full dataset
      workingIndices = Array.from({length: totalLength}, (_, i) => i);
      workingColumn = this.targetColumn!;
      this.sampledIndices = null;
    }

    // Get the selected distance function and its options
    const distanceFunction = this.distanceMetric as MmDistanceFunctionsNames;
    const distanceFunctionOptions = this.getDistanceFunctionOptions();

    // Encode sequences using the selected distance function
    // Note: getEncodedSeqSpaceCol only needs the fingerprint type, not full options
    const fingerprintType = distanceFunctionOptions.fingerprintType || 'Morgan';
    const encodedSequences = (await getEncodedSeqSpaceCol(
      workingColumn,
      distanceFunction,
      fingerprintType
    )).seqList;

    // Calculate distance matrix for the working dataset
    const distanceMatrixService = new DistanceMatrixService(true, false);
    const distanceMatrixData = await distanceMatrixService.calc(
      encodedSequences,
      distanceFunction,
      true, // normalize
      distanceFunctionOptions
    );
    distanceMatrixService.terminate();

    const workingLength = workingColumn.length;
    const linearizeFunc = dmLinearIndex(workingLength);

    // Apply diversity selection on the working dataset
    const diverseIndicesInWorkingSet = getDiverseSubset(
      workingLength,
      Math.min(workingLength, this.limit),
      (i1: number, i2: number) => {
        return workingColumn.isNone(i1) || workingColumn.isNone(i2) ? 0 :
          distanceMatrixData[linearizeFunc(i1, i2)];
      }
    );

    // Map back to original indices
    this.renderMolIds = diverseIndicesInWorkingSet.map((workingIndex) => workingIndices[workingIndex]);
  }

  private createRandomSample(totalLength: number, sampleSize: number): number[] {
    // Create array of all indices
    const allIndices = Array.from({length: totalLength}, (_, i) => i);

    // Filter out indices with missing values to avoid wasting sample space
    const validIndices = allIndices.filter((i) => !this.targetColumn!.isNone(i));

    if (validIndices.length <= sampleSize)
      return validIndices;


    // Fisher-Yates shuffle to get random sample
    const shuffled = [...validIndices];
    for (let i = shuffled.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }

    return shuffled.slice(0, sampleSize).sort((a, b) => a - b);
  }

  // Helper method to get information about sampling (useful for debugging/info)
  getSamplingInfo(): {isSampled: boolean, originalSize?: number, sampleSize?: number, sampledIndices?: number[]} {
    if (this.sampledIndices) {
      return {
        isSampled: true,
        originalSize: this.targetColumn?.length,
        sampleSize: this.sampledIndices.length,
        sampledIndices: this.sampledIndices
      };
    }
    return {isSampled: false};
  }
}
