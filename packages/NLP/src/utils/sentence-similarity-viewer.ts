import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { SearchBaseViewer } from '@datagrok-libraries/ml/src/viewers/search-base-viewer';
import { Subject } from 'rxjs';

export class SentenceSimilarityViewer extends SearchBaseViewer {
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  gridSelect: boolean = false;
  targetMoleculeIdx: number = 0;
  computeCompleted = new Subject<boolean>();
  splashScreen: HTMLElement | null = null;
  //threshold: number;

  constructor() {
    super('similarity', DG.SEMTYPE.TEXT);
    //this.threshold = this.float('threshold', 0.8, {min: 0, max: 1});
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

        const resDf = DG.DataFrame.fromCsv(await grok.functions.call('NLP: getSimilar', {embeddingsStr: embeddings, sentences: this.targetColumn.toList(), currentRow: this.targetMoleculeIdx}));
        //resDf.rows.removeWhere((row) => row.get('score') < this.threshold);
        resDf.rows.removeWhereIdx((idx) => idx > this.limit);

        resDf.onCurrentRowChanged.subscribe((_: any) => {
          this.dataFrame.currentRowIdx = resDf.col('idx')!.get(resDf.currentRowIdx);
          this.gridSelect = true;
        });

        const grid = resDf.plot.grid();
        grid.col('idx')!.visible = false;

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