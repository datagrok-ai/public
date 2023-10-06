import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ScaffoldTreeViewer, BitwiseOp, isOrphans} from './scaffold-tree';
import {filter, debounce} from 'rxjs/operators';
import {interval} from 'rxjs';
import { SCAFFOLD_TREE_HIGHLIGHT } from '../constants';

const COLUMN_NAME_CHANGED = 'column-name-changed';

interface IScaffoldFilterState {
  colName: string;
}

function clearNotIcon(viewer: ScaffoldTreeViewer, tree: DG.TreeViewNode[]) {
  tree.map((group) => {
    const castedGroup = group as DG.TreeViewGroup;
    if (!isOrphans(castedGroup)) {
      viewer.setNotBitOperation(castedGroup, false);
      if (castedGroup.children)
        clearNotIcon(viewer, castedGroup.children)
    }
  });
}

export class ScaffoldTreeFilter extends DG.Filter {
  viewer: ScaffoldTreeViewer = new ScaffoldTreeViewer();
  savedTree: string = '';
  colorCodedScaffolds: string = '';

  constructor() {
    super();
    this.root = ui.divV([]);
    this.subs = this.viewer.subs;
    this.subs.push(grok.events.onResetFilterRequest.subscribe((_) => {
      this.viewer.clearFilters();
      clearNotIcon(this.viewer, this.viewer.tree.children);
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
      .pipe(debounce((_) => interval(100)))
      .subscribe((_) => { 
        this.column!.setTag(SCAFFOLD_TREE_HIGHLIGHT, JSON.stringify(this.viewer.colorCodedScaffolds));
        if (!this.isFiltering) {
          this.viewer.bitset = null;
          this.dataFrame!.rows.requestFilter();
        } else {
          this.viewer.updateFilters();
        }
      }),
    );
    this.subs.push(grok.events.onCustomEvent(COLUMN_NAME_CHANGED).subscribe((state: IScaffoldFilterState) => {
      if (state.colName === this.columnName)
        this.createViewer(dataFrame);
    }));
  }

  saveState(): any {
    const state = super.saveState();
    state.savedTree = JSON.stringify(this.viewer.serializeTrees(this.viewer.tree));
    state.colorCodedScaffolds = JSON.stringify(this.viewer.colorCodedScaffolds);
    return state;
  }

  applyState(state: any): void {
    super.applyState(state);
    grok.events.fireCustomEvent(COLUMN_NAME_CHANGED, {
      colName: state.columnName,
    });

    if (state.savedTree)
      this.viewer.loadTreeStr(state.savedTree);

    if (state.colorCodedScaffolds)
      this.viewer.colorCodedScaffolds = JSON.parse(state.colorCodedScaffolds);
  }

  applyFilter(): void {
  }

  detach(): void {
    super.detach();
    this.viewer.clearFilters();
    this.viewer.molCol!.setTag(SCAFFOLD_TREE_HIGHLIGHT, '');
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
