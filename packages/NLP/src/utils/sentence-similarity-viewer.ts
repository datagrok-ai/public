import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SearchBaseViewer} from '@datagrok-libraries/ml/src/viewers/search-base-viewer';
import {KnnResult, SparseMatrixService} from '@datagrok-libraries/ml/src/distance-matrix/sparse-matrix-service';
import {Subject} from 'rxjs';
import {VectorMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {multiColWebGPUKNN} from '@datagrok-libraries/math/src/webGPU/multi-col-knn/multiCol-KNN';
import {WEBGPUDISTANCE} from '@datagrok-libraries/math/src/webGPU/multi-col-distances/webGPU-multicol-distances';
import {WEBGSLAGGREGATION} from '@datagrok-libraries/math/src/webGPU/multi-col-distances/webGPU-aggregation';
import {getEmbeddings} from '../package';
import {DistanceAggregationMethods} from '@datagrok-libraries/ml/src/distance-matrix/types';

export class SentenceSimilarityViewer extends SearchBaseViewer {
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  targetMoleculeIdx: number = 0;
  computeCompleted = new Subject<boolean>();
  splashScreen: HTMLElement | null = null;
  knn?: KnnResult;
  sentenceCol: DG.Column<string> | null = null;
  additionalColumnNames: string[];
  indexWScore: { idx: number; score: number; }[] | null = null;
  prevTargetColumnName: string | null = null;
  embeddings: Float32Array[] | null = null;
  similarDf: DG.DataFrame | null = null;
  similarGrid: DG.Grid | null = null;

  constructor() {
    super('similarity', DG.SEMTYPE.TEXT);
    this.additionalColumnNames = this.addProperty('additionalColumnNames', DG.TYPE.COLUMN_LIST, null,
      {description: 'Additional columns to display with similar sentences'});
    this.recomputeOnCurrentRowChange = false;
    this.skipRecomputingProperies.push('additionalColumnNames');
  }

  init(): void {
    this.initialized = true;
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();
    this.sub(
      grok.events.onContextMenu.subscribe((e) => {
        if (e.causedBy && e.causedBy.target && this.similarGrid?.root?.contains(e.causedBy.target) &&
        this.similarDf && e.args?.menu && e.args.menu.item) {
          e.causedBy.preventDefault();
          e.causedBy.stopPropagation();
          e.causedBy.stopImmediatePropagation();
          e.args.menu.item('Download Results',
            () => DG.Utils.download('similar_sentences.csv', this.similarDf!.toCsv()));
        }
      }),
    );
  }

  override async renderInt(computeData: boolean): Promise<void> {
    if (!this.beforeRender()) return;

    if (!this.targetColumn) {
      this.showErrorMessage(`There is no ${this.semType} column available.`);
      return;
    }

    const currentRowIdx = this.dataFrame!.currentRowIdx;
    this.curIdx = currentRowIdx === -1 ? 0 : currentRowIdx;

    if (computeData) {
      this.targetMoleculeIdx = currentRowIdx === -1 ? 0 : currentRowIdx;
      ui.empty(this.root);
      this.showSplashScreen('Getting embeddings ...');

      try {
        this.updateSplashScreen('Finding similar sentences...');
        await this.compute();

        const actualLimit = Math.min(Math.min(this.limit, this.targetColumn!.length - 1), this.indexWScore!.length - 1);
        const idxsData = new Int32Array(actualLimit + 1);
        const scoresData = new Float32Array(actualLimit + 1);
        const sentencesData: string[] = new Array(actualLimit + 1);

        for (let i = 0; i <= actualLimit; i++) {
          idxsData[i] = this.indexWScore![i].idx;
          scoresData[i] = this.indexWScore![i].score;
          sentencesData[i] = this.targetColumn?.get(idxsData[i]) || '';
        }

        this.idxs = DG.Column.fromInt32Array('indexes', idxsData);
        this.scores = DG.Column.fromFloat32Array('score', scoresData);
        this.sentenceCol = DG.Column.fromStrings(this.targetColumnName, sentencesData);
        this.sentenceCol.semType = DG.SEMTYPE.TEXT;

        const additionalColumns = this.additionalColumnNames ?
          this.additionalColumnNames.reduce((acc: DG.Column[], colName: string) => {
            if (colName !== this.targetColumnName) {
              const col = this.dataFrame.columns.byName(colName);
              // eslint-disable-next-line max-len
              acc.push(DG.Column.fromType(col.type, col.name, this.idxs!.length).init((i) => col.get(this.idxs!.get(i))));
            }
            return acc;
          }, []) :
          [];

        this.similarDf = DG.DataFrame.fromColumns([this.idxs!, this.sentenceCol!, this.scores!, ...additionalColumns]);
        this.similarDf.onCurrentRowChanged.subscribe((_: any) => {
          this.dataFrame.currentRowIdx = this.similarDf!.col('indexes')!.get(this.similarDf!.currentRowIdx);
          this.gridSelect = true;
        });

        this.similarGrid = this.similarDf.plot.grid();
        this.similarGrid.props.rowHeight = 80;
        this.similarGrid.col('indexes')!.visible = false;

        const view = grok.shell.v as DG.TableView;
        view.grid.root.addEventListener('click', (_event: MouseEvent) => {
          this.gridSelect = false;
        });

        this.updateDivInnerHTML(this.root, this.similarGrid.root);
        this.computeCompleted.next(true);
      } catch (error) {
        this.showErrorMessage('An error occurred during the computation');
        console.error(error);
      }
      this.removeSplashScreen();
    }
  }

  async getEmbeddings(): Promise<Float32Array[]> {
    if (!this.embeddings || this.targetColumnName !== this.prevTargetColumnName)
      this.embeddings = await getEmbeddings(this.targetColumn!.toList());
    return this.embeddings;
  }

  private async compute() {
    if (!this.targetColumn)
      return;
    const len = Math.min(this.maxLimit, this.targetColumn!.length);

    const needsKnnRecompute = this.targetColumn && (!this.knn || this.targetColumnName !== this.prevTargetColumnName);
    // if data is more than 10000 rows, KNN computation can take a long time, so we will compute on fly
    const calulateOnTheFly = this.targetColumn.length > 5000;
    if (needsKnnRecompute && !calulateOnTheFly) {
      const embeddings = await this.getEmbeddings();
      try {
        this.knn = await multiColWebGPUKNN(
          // eslint-disable-next-line max-len
          [embeddings], len - 1, [WEBGPUDISTANCE.VECTOR_COSINE], WEBGSLAGGREGATION.MANHATTAN,
          [1], [{'preprocessingFuncArgs': {}}],
        ) as unknown as KnnResult;
      } catch (e) {
        console.error(e);
      }

      if (!this.knn) {
        this.knn = await (new SparseMatrixService()
          .getKNN(embeddings, VectorMetricsNames.Cosine, len - 1));
      }
    } else if (calulateOnTheFly) {
      const embeddings = await this.getEmbeddings();
      const localNN = await (new SparseMatrixService()
        .multiColumnSingleValueKNN([embeddings], this.targetMoleculeIdx,
          [VectorMetricsNames.Cosine], this.maxLimit, [{}], [1], DistanceAggregationMethods.MANHATTAN,
        ));
      this.indexWScore = Array.from({length: localNN.knnDistances.length}, (_, i) => ({
        idx: localNN!.knnIndexes[i],
        score: 1 - localNN!.knnDistances[i],
      }));
    } else if (this.knn) {
      this.indexWScore = Array.from({length: len - 1}, (_, i) => ({
        idx: this.knn!.knnIndexes[this.targetMoleculeIdx][i],
        score: 1 - this.knn!.knnDistances[this.targetMoleculeIdx][i],
      }));
    }
    this.indexWScore ??= [];
    this.indexWScore.sort((a, b) => b.score - a.score);
    // knn can return non-meaningful indexes like -1 or 99999, so we need to filter them out
    this.indexWScore = this.indexWScore
      .filter((el) => el.idx != null && !isNaN(el.idx) && el.idx >= 0 &&
        el.idx < this.targetColumn!.length && !isNaN(el.score) && el.score != null &&
        el.idx !== this.targetMoleculeIdx);
    this.indexWScore.unshift({idx: this.targetMoleculeIdx, score: DG.FLOAT_NULL});
    this.prevTargetColumnName = this.targetColumnName;
  }

  private showSplashScreen(message: string): void {
    if (!this.splashScreen) {
      this.splashScreen = ui.div(
        [ui.loader(), ui.divText(message, {style: {'margin-top': '50px'}})],
        'nlp-splash-container',
      );

      this.root.appendChild(this.splashScreen);
    } else
      this.splashScreen.style.display = 'flex';
  }

  private updateSplashScreen(message: string): void {
    if (this.splashScreen && this.splashScreen.children[1]!.textContent !== message)
      this.splashScreen.children[1]!.textContent = message;
  }

  private removeSplashScreen(): void {
    if (this.splashScreen) {
      this.splashScreen.style.display = 'none';
      this.splashScreen = null;
    }
  }

  showErrorMessage(message: string): void {
    const errorMessageDiv = ui.div([
      ui.divText(message, {style: {'color': 'red'}}),
    ], 'nlp-splash-container');

    this.root.appendChild(errorMessageDiv);
  }

  updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
    div.innerHTML = '';
    div.append(content);
  }
}
