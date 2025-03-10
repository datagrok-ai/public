import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {ReactionData, Tree} from './aizynth-api';
import {SAMPLE_TREE} from './mock-data';
import {RotatePath} from './const';
import {createPathsTreeTabs} from './utils';

export class AiZynthFinderViewer extends DG.JsViewer {
  scoreCutoff: number;
  moleculeColumn?: DG.Column|null;
  currentMolecule = '';
  rotate = 'Clockwise';
  opts = {suppressChiralText: true};
  // padding = 4;
  // molSize = 120;
  // verticalSpacing = 20;
  // horizontalSpacing = 30;
  paths: Tree[] = SAMPLE_TREE;

  constructor() {
    super();
    this.scoreCutoff = this.float('scoreCutoff', 0.99, {min: 0, max: 1});
    this.rotate = this.string('rotate', RotatePath.Clockwise, {choices: Object.keys(RotatePath)});
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    if (this.dataFrame) {
      this.moleculeColumn = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      if (!this.moleculeColumn) {
        grok.shell.error(`Dataframe doesn't contain molecule columns`);
        return;
      }
      if (this.dataFrame.currentRowIdx !== -1)
        this.currentMolecule = this.moleculeColumn.get(this.dataFrame.currentRowIdx);
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50)
        .subscribe(async (_: any) => {
          this.currentMolecule = this.moleculeColumn?.get(this.dataFrame.currentRowIdx);
          this.renderPaths();
        }));

      this.renderPaths();
    }
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
  }


  async renderPaths(computeData = true) {
    ui.empty(this.root);
    if (computeData) {
      ui.setUpdateIndicator(this.root, true, `Generating retrosynthesis paths`);
      const result = await grok.functions.call('Aizynthfinder:calculateRetroSynthesisPaths',
        {molecule: this.currentMolecule});
      const reactionData = JSON.parse(result) as ReactionData;
      if (reactionData.data?.length) {
        const paths = reactionData.data[0].trees as Tree[];
        if (paths?.length) {
          const tabControl = createPathsTreeTabs(this.paths);
          this.root.append(tabControl.root);
        } else
          this.root.append(ui.divText('No paths found for the molecule'));
      } else
        this.root.append(ui.divText('No paths found for the molecule'));
    }
    ui.setUpdateIndicator(this.root, false);
  }
}
