import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { SearchBaseViewer } from '@datagrok-libraries/ml/src/viewers/search-base-viewer';
import { KnnResult, SparseMatrixService } from '@datagrok-libraries/ml/src/distance-matrix/sparse-matrix-service';
import { Subject } from 'rxjs';
import { VectorMetricsNames } from '@datagrok-libraries/ml/src/typed-metrics';
import { multiColWebGPUKNN } from '@datagrok-libraries/math/src/webGPU/multi-col-knn/multiCol-KNN'
import { WEBGPUDISTANCE } from '@datagrok-libraries/math/src/webGPU/multi-col-distances/webGPU-multicol-distances';
import { WEBGSLAGGREGATION } from '@datagrok-libraries/math/src/webGPU/multi-col-distances/webGPU-aggregation';

export class SentenceSimilarityViewer extends SearchBaseViewer {
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  gridSelect: boolean = false;
  targetMoleculeIdx: number = 0;
  computeCompleted = new Subject<boolean>();
  splashScreen: HTMLElement | null = null;
  knn?: KnnResult;
  kPrevNeighbors: number = 0;
  sentenceCol: DG.Column<string> | null = null;

  constructor() {
    super('similarity', DG.SEMTYPE.TEXT);
  }

  init(): void {
    this.initialized = true;
  }

  override async renderInt(computeData: boolean): Promise<void> {
    if (!this.beforeRender()) return;

    if (!this.targetColumn) {
      const noTargetColumnMessage = ui.divV([
        ui.divText(`There is no ${this.semType} column available.`, { style: { 'color': 'red', 'margin-top': '20px' } })],
        'nlp-splash-container'
      );
      this.updateDivInnerHTML(this.root, noTargetColumnMessage);
      return;
    }

    this.curIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;

    if (computeData && !this.gridSelect) {
      this.targetMoleculeIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;

      ui.empty(this.root);

      this.showSplashScreen("Getting embeddings ...");

      try {
        const embeddings = await grok.functions.call('NLP: getEmbeddings', { sentences: this.targetColumn.toList() });

        this.updateSplashScreen("Finding similar sentences...");
        await this.compute(JSON.parse(embeddings));

        this.sentenceCol = DG.Column.string(this.targetColumnName,
          this.idxs!.length).init((i) => this.targetColumn?.get(this.idxs?.get(i)));
        this.sentenceCol.semType = DG.SEMTYPE.TEXT;
        const resDf = DG.DataFrame.fromColumns([this.idxs!, this.sentenceCol!, this.scores!]);
        resDf.onCurrentRowChanged.subscribe((_: any) => {
          this.dataFrame.currentRowIdx = resDf.col('indexes')!.get(resDf.currentRowIdx);
          this.gridSelect = true;
        });

        const grid = resDf.plot.grid();
        grid.col('indexes')!.visible = false;

        const view = grok.shell.v as DG.TableView;
        view.grid.root.addEventListener('click', (_event: MouseEvent) => {
          this.gridSelect = false;
        });

        this.updateDivInnerHTML(this.root, grid.root);
        this.computeCompleted.next(true);
      } catch (error) {
        console.error("Error during computation:", error);
      } finally {
        this.removeSplashScreen();
      }
    }
  }

  private async compute(embeddings: Array<Array<number>>) {
    const len = this.targetColumn!.length;
    const actualLimit = Math.min(this.limit, len - 1);

    if (!this.knn || this.kPrevNeighbors !== actualLimit) {
      try {
        this.knn = await multiColWebGPUKNN(
          [embeddings.map((c) => new Float32Array(c))], Math.min(this.limit, len - 1), [WEBGPUDISTANCE.VECTOR_COSINE], WEBGSLAGGREGATION.MANHATTAN,
          [1], [{'preprocessingFuncArgs': {}}]
        ) as unknown as KnnResult;
      } catch (e) {
        console.error(e);
      }
      
      if (!this.knn) {
        this.kPrevNeighbors = actualLimit;
        this.knn = await (new SparseMatrixService()
          .getKNN(embeddings, VectorMetricsNames.Cosine, Math.min(this.limit, len - 1)));
      }
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

  private showSplashScreen(message: string): void {
    if (!this.splashScreen) {
      this.splashScreen = ui.div(
        [ui.loader(), ui.divText(message, {style: {'margin-top': '50px'}})],
        'nlp-splash-container'
      );

      this.root.appendChild(this.splashScreen);
    } else {
      this.splashScreen.style.display = 'flex';
    }
  }

  private updateSplashScreen(message: string): void {
    if (this.splashScreen) {
      this.splashScreen.children[1]!.textContent = message;
    }
  }

  private removeSplashScreen(): void {
    if (this.splashScreen) {
      this.splashScreen.style.display = 'none';
    }
  }

  updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
    div.innerHTML = '';
    div.append(content);
  }
}