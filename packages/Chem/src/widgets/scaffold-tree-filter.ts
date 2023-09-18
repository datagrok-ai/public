import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ScaffoldTreeViewer, BitwiseOp} from './scaffold-tree';
import {filter, debounce} from 'rxjs/operators';
import {interval} from 'rxjs';
import { FILTER_SCAFFOLD_TAG, HIGHLIGHT_BY_SCAFFOLD_TAG } from '../constants';
import { IColoredScaffold } from '../rendering/rdkit-cell-renderer';

const COLUMN_NAME_CHANGED = 'column-name-changed';

interface IScaffoldFilterState {
  colName: string;
}

function clearNotIcon(viewer: ScaffoldTreeViewer, tree: DG.TreeViewNode[]) {
  tree.map((group) => {
    const castedGroup = group as DG.TreeViewGroup;
    viewer.setNotBitOperation(castedGroup, false);
    if (castedGroup.children) {
      clearNotIcon(viewer, castedGroup.children)
    }
  });
}

export class ScaffoldTreeFilter extends DG.Filter {
  viewer: ScaffoldTreeViewer = new ScaffoldTreeViewer();
  savedTree: string = '';
  scaffolds: IColoredScaffold[] = [];

  constructor() {
    super();
    this.root = ui.divV([]);
    this.subs = this.viewer.subs;
    this.subs.push(grok.events.onResetFilterRequest.subscribe((_) => {
      this.viewer.clearFilters();
      clearNotIcon(this.viewer, this.viewer.tree.children);
      if (this.viewer._bitOpInput)
        this.viewer._bitOpInput.value = BitwiseOp.OR;
    }));
  }

  get isFiltering(): boolean {
    return super.isFiltering && this.viewer.tree.items.filter((v) => v.checked).length > 0;
  }

  get filterSummary(): string {
    return ScaffoldTreeViewer.TYPE;
  }

  attach(dataFrame: DG.DataFrame): void {
    super.attach(dataFrame);
    this.column ??= dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.columnName ??= this.column?.name;
    this.subs.push(this.dataFrame!.onRowsFiltering
      .pipe(filter((_) => !this.isFiltering), debounce((_) => interval(100)))
      .subscribe((_) => {
        //delete this.column!.temp[FILTER_SCAFFOLD_TAG];
        this.column!.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, '');
        this.viewer.updateFilters(this.isFiltering);
      }),
    );
    this.subs.push(grok.events.onCustomEvent(COLUMN_NAME_CHANGED).subscribe((state: IScaffoldFilterState) => {
      if (state.colName === this.columnName)
        this.createViewer(dataFrame);
    }));
  }

  saveState(): any {
    const state = super.saveState();
    state.savedTree = JSON.stringify(ScaffoldTreeViewer.serializeTrees(this.viewer.tree));
    state.scaffolds = this.viewer.scaffolds;
    return state;
  }

  applyState(state: any): void {
    super.applyState(state);
    grok.events.fireCustomEvent(COLUMN_NAME_CHANGED, {
      colName: state.columnName,
    });
    if (state.savedTree) {
      this.viewer.loadTreeStr(state.savedTree);
    }
    if (state.scaffolds) {
      this.viewer.molCol!.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, JSON.stringify(state.scaffolds));
    }
  }

  applyFilter(): void {
  }

  detach(): void {
    super.detach();
    this.viewer.clearFilters();
    this.viewer.molCol!.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, '');
     grok.shell.tv.dataFrame.fireValuesChanged();
  }

  createViewer(dataFrame: DG.DataFrame) {
    this.viewer.allowGenerate = false;
    this.viewer.size = 'extra';
    this.viewer.resizable = true;
    this.viewer.molCol = this.column;
    this.viewer.MoleculeColumn = this.columnName!;
    this.viewer.dataFrame = dataFrame;
    this.root.appendChild(this.viewer.root);
  }
}
