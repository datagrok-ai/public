import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SequenceSearchBaseViewer} from './sequence-search-base-viewer';
import {getMonomericMols} from '../calculations/monomerLevelMols';
import {createDifferenceCanvas, createDifferencesWithPositions} from './sequence-activity-cliffs';
import {updateDivInnerHTML} from '../utils/ui-utils';
import {Subject} from 'rxjs';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {alignSequencePair} from '@datagrok-libraries/bio/src/utils/macromolecule/alignment';
import {KnnResult, SparseMatrixService} from '@datagrok-libraries/ml/src/distance-matrix/sparse-matrix-service';
import {getEncodedSeqSpaceCol} from './sequence-space';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

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

  constructor(demo?: boolean) {
    super('similarity');
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
    if (this.moleculeColumn) {
      this.curIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;
      if (computeData && !this.gridSelect) {
        this.targetMoleculeIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;
        const sh = SeqHandler.forColumn(this.moleculeColumn!);

        await (!sh.isHelm() ? this.computeByMM() : this.computeByChem());
        const similarColumnName: string = this.similarColumnLabel != null ? this.similarColumnLabel :
          `similar (${this.moleculeColumnName})`;
        this.molCol = DG.Column.string(similarColumnName,
          this.idxs!.length).init((i) => this.moleculeColumn?.get(this.idxs?.get(i)));
        this.molCol.semType = DG.SEMTYPE.MACROMOLECULE;
        this.tags.forEach((tag) => this.molCol!.setTag(tag, this.moleculeColumn!.getTag(tag)));
        const resDf = DG.DataFrame.fromColumns([this.idxs!, this.molCol!, this.scores!]);
        resDf.onCurrentRowChanged.subscribe((_: any) => {
          this.dataFrame.currentRowIdx = resDf.col('indexes')!.get(resDf.currentRowIdx);
          setTimeout(() => { this.createPropertyPanel(resDf); }, 1000);
          this.gridSelect = true;
        });
        const grid = resDf.plot.grid();
        grid.col('indexes')!.visible = false;
        const targetMolRow = this.idxs?.getRawData().findIndex((it) => it == this.targetMoleculeIdx);
        const targetScoreCell = grid.cell('score', targetMolRow!);
        targetScoreCell.cell.value = null;
        const view = this.demo ? (grok.shell.view('Browse')! as DG.BrowseView)!.preview! as DG.TableView : grok.shell.v as DG.TableView;
        view.grid.root.addEventListener('click', (_event: MouseEvent) => {
          this.gridSelect = false;
        });
        updateDivInnerHTML(this.root, grid.root);
        this.computeCompleted.next(true);
      }
    }
  }

  private async computeByChem() {
    const monomericMols = await getMonomericMols(this.moleculeColumn!);
    //need to create df to calculate fingerprints
    const _monomericMolsDf = DG.DataFrame.fromColumns([monomericMols]);
    const df = await grok.functions.call('Chem:callChemSimilaritySearch', {
      df: this.dataFrame,
      col: monomericMols,
      molecule: monomericMols.get(this.targetMoleculeIdx),
      metricName: this.distanceMetric,
      limit: this.limit,
      minScore: this.cutoff,
      fingerprint: this.fingerprint,
    });
    this.idxs = df.getCol('indexes');
    this.scores = df.getCol('score');
  }

  private async computeByMM() {
    const len = this.moleculeColumn!.length;
    const actualLimit = Math.min(this.limit, len - 1);
    if (!this.knn || this.kPrevNeighbors !== actualLimit) {
      const encodedSequences =
        (await getEncodedSeqSpaceCol(this.moleculeColumn!, MmDistanceFunctionsNames.LEVENSHTEIN)).seqList;

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
    const molColSh = SeqHandler.forColumn(this.moleculeColumn!);
    const resSh = SeqHandler.forColumn(resCol);
    const subParts1 = molColSh.getSplitted(this.targetMoleculeIdx);
    const subParts2 = resSh.getSplitted(resDf.currentRowIdx);
    const alignment = alignSequencePair(subParts1, subParts2);
    const canvas = createDifferenceCanvas(alignment.seq1Splitted, alignment.seq2Splitted, resSh.units, molDifferences);
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
