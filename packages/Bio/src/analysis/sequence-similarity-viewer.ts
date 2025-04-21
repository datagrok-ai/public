import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SequenceSearchBaseViewer} from './sequence-search-base-viewer';
import {createDifferenceCanvas, createDifferencesWithPositions} from './sequence-activity-cliffs';
import {adjustGridcolAfterRender, updateDivInnerHTML} from '../utils/ui-utils';
import {Subject} from 'rxjs';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {alignSequencePair} from '@datagrok-libraries/bio/src/utils/macromolecule/alignment';
import {KnnResult, SparseMatrixService} from '@datagrok-libraries/ml/src/distance-matrix/sparse-matrix-service';
import {getEncodedSeqSpaceCol} from './sequence-space';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {MmcrTemps, tempTAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';

export class SequenceSimilarityViewer extends SequenceSearchBaseViewer {
  cutoff: number;
  hotSearch: boolean;
  similarColumnLabel: string | null; // Use postfix Label to prevent activating table column selection editor

  sketchedMolecule: string = '';
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  gridSelect: boolean = false;
  targetMoleculeIdx: number = 0;
  computeCompleted = new Subject<boolean>();
  distanceMatrixComputed: boolean = false;
  mmDistanceMatrix: Float32Array;
  knn?: KnnResult;
  kPrevNeighbors: number = 0;
  demo?: boolean;
  analysisGrid?: DG.Grid;
  subInited: boolean = false;

  constructor(
    private readonly seqHelper: ISeqHelper,
    demo?: boolean,
  ) {
    super('similarity', DG.SEMTYPE.MACROMOLECULE);
    this.cutoff = this.float('cutoff', 0.01, {min: 0, max: 1});
    this.hotSearch = this.bool('hotSearch', true);
    this.similarColumnLabel = this.string('similarColumnLabel', null);
    this.demo = demo;
  }

  init(): void {
    this.hotSearch = true;
    this.initialized = true;
  }

  override async renderInt(computeData: boolean): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.targetColumn) {
      this.curIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;
      if (computeData && !this.gridSelect) {
        this.targetMoleculeIdx = (this.dataFrame!.currentRowIdx ?? -1) < 0 ? 0 : this.dataFrame!.currentRowIdx;

        await this.computeByMM();
        const similarColumnName: string = this.similarColumnLabel != null ? this.similarColumnLabel :
          `similar (${this.targetColumn})`;
        this.molCol = DG.Column.string(similarColumnName,
          this.idxs!.length).init((i) => this.targetColumn?.get(this.idxs?.get(i)));
        this.molCol.semType = DG.SEMTYPE.MACROMOLECULE;
        this.tags.forEach((tag) => this.molCol!.setTag(tag, this.targetColumn!.getTag(tag)));
        const resDf = DG.DataFrame.fromColumns([this.idxs!, this.molCol!, this.scores!]);
        await resDf.meta.detectSemanticTypes();
        await grok.data.detectSemanticTypes(resDf);
        this.molCol.temp[tempTAGS.referenceSequence] = this.targetColumn!.get(this.targetMoleculeIdx);
        this.molCol.temp[MmcrTemps.maxMonomerLength] = 4;
        let prevTimer: any = null;
        const _ = resDf.onCurrentRowChanged.subscribe((_: any) => {
          prevTimer && clearTimeout(prevTimer);
          this.dataFrame.currentRowIdx = resDf.col('indexes')!.get(resDf.currentRowIdx);
          prevTimer = setTimeout(() => { this.createPropertyPanel(resDf); }, 300);
          this.gridSelect = true;
        });
        if (!this.analysisGrid) {
          this.analysisGrid = resDf.plot.grid();
          updateDivInnerHTML(this.root, this.analysisGrid.root);
        } else {
          this.analysisGrid.dataFrame = resDf;
          this.analysisGrid.invalidate();
        }
        this.analysisGrid.col('indexes')!.visible = false;
        adjustGridcolAfterRender(this.analysisGrid, this.molCol!.name, 450, 30, true);
        const targetMolRow = this.idxs?.getRawData().findIndex((it) => it == this.targetMoleculeIdx);
        const targetScoreCell = this.analysisGrid.cell('score', targetMolRow!);
        targetScoreCell.cell.value = null;
        const view = grok.shell.tv;
        if (!this.subInited) {
          view.grid.root.addEventListener('click', (_event: MouseEvent) => {
            this.gridSelect = false;
          });
          this.subInited = true;
        }
        this.computeCompleted.next(true);
      }
    }
  }

  private async computeByMM() {
    const len = this.targetColumn!.length;
    const actualLimit = Math.min(this.limit, len - 1);
    if (!this.knn || this.kPrevNeighbors !== actualLimit) {
      const encodedSequences =
        (await getEncodedSeqSpaceCol(this.targetColumn!, MmDistanceFunctionsNames.LEVENSHTEIN)).seqList;

      this.kPrevNeighbors = actualLimit;
      this.knn = await (new SparseMatrixService()
        .getKNN(encodedSequences, MmDistanceFunctionsNames.LEVENSHTEIN, Math.min(this.limit, len - 1)));
    }
    const indexWScore = new Array(actualLimit).fill(0).map((_, i) => ({
      idx: this.knn!.knnIndexes[this.targetMoleculeIdx][i],
      score: 1 - this.knn!.knnDistances[this.targetMoleculeIdx][i],
    }));
    indexWScore.sort((a, b) => b.score - a.score);
    indexWScore.unshift({idx: this.targetMoleculeIdx, score: DG.FLOAT_NULL});
    this.idxs = DG.Column.int('indexes', actualLimit + 1).init((i) => indexWScore[i].idx);
    this.scores = DG.Column.float('score', actualLimit + 1).init((i) => indexWScore[i].score);
  }

  createPropertyPanel(resDf: DG.DataFrame) {
    const propPanel = ui.div();
    const molDifferences: { [key: number]: HTMLCanvasElement } = {};
    const molColName = this.molCol?.name!;
    const resCol: DG.Column<string> = resDf.col(molColName)!;
    const molColSh = this.seqHelper.getSeqHandler(this.targetColumn!);
    const resSh = this.seqHelper.getSeqHandler(resCol);
    const subParts1 = molColSh.getSplitted(this.targetMoleculeIdx);
    const subParts2 = resSh.getSplitted(resDf.currentRowIdx);
    const alignment = alignSequencePair(subParts1, subParts2);
    const canvas =
      createDifferenceCanvas(alignment.seq1Splitted, alignment.seq2Splitted, resSh.defaultBiotype, molDifferences);
    propPanel.append(ui.div(canvas, {style: {width: '300px', overflow: 'scroll'}}));
    if (subParts1.length !== subParts2.length) {
      propPanel.append(ui.divV([
        ui.divText(`Different sequence length:`, {style: {fontWeight: 'bold'}}),
        ui.divText(`target: ${subParts1.length} monomers`),
        ui.divText(`selected: ${subParts2.length} monomers`),
      ], {style: {paddingBottom: '10px'}}));
    }
    propPanel.append(createDifferencesWithPositions(molDifferences));
    const acc = ui.accordion();
    const accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(`Similarity search`)]));
    acc.addPane('Differences', () => propPanel, true);
    grok.shell.o = acc.root;
  }
}
