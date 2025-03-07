import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { SearchBaseViewer } from '@datagrok-libraries/ml/src/viewers/search-base-viewer';
import { Subject } from 'rxjs';
import { getEmbeddings, getSimilar } from './package';

export class SentenceSimilarityViewer extends SearchBaseViewer {
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  gridSelect: boolean = false;
  targetMoleculeIdx: number = 0;
  computeCompleted = new Subject<boolean>();

  constructor() {
    super('similarity', DG.SEMTYPE.TEXT);
  }

  init(): void {
    this.initialized = true;
  }

  override async renderInt(computeData: boolean): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.targetColumn) {
      this.curIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;
      if (computeData && !this.gridSelect) {
        this.targetMoleculeIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;

        const embeddings = await getEmbeddings(this.targetColumn);
        const resDf = DG.DataFrame.fromCsv(await getSimilar(embeddings, this.targetColumn, this.targetMoleculeIdx));
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
      }
    }
  }

  updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
    div.innerHTML = '';
    div.append(content);
  }
}