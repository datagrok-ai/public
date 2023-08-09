import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ScaffoldTreeViewer, BitwiseOp} from './scaffold-tree';
import {filter, debounce} from 'rxjs/operators';
import {interval} from 'rxjs';

const COLUMN_NAME_CHANGED = 'column-name-changed';

interface IScaffoldFilterState {
  colName: string;
}

export class ScaffoldTreeFilter extends DG.Filter {
  viewer: ScaffoldTreeViewer = new ScaffoldTreeViewer();
  savedTree: string = '';

  constructor() {
    super();
    this.root = ui.divV([]);
    this.subs = this.viewer.subs;
    this.subs.push(grok.events.onResetFilterRequest.subscribe((_) => {
      this.viewer.clearFilters();
      this.viewer.tree.children.map((group) => {
        this.viewer.setNotBitOperation(group as DG.TreeViewGroup, false);
      });
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
        delete this.column!.temp['chem-scaffold-filter'];
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
    return state;
  }

  applyState(state: any): void {
    super.applyState(state);
    grok.events.fireCustomEvent(COLUMN_NAME_CHANGED, {
      colName: state.columnName,
    });
    if (state.savedTree) this.viewer.loadTreeStr(state.savedTree);
  }

  applyFilter(): void {
  }

  detach(): void {
    super.detach();
    this.viewer.clearFilters();
  }

  createViewer(dataFrame: DG.DataFrame) {
    this.viewer.allowGenerate = false;
    this.viewer.size = 'large';
    this.viewer.resizable = true;
    this.viewer.molCol = this.column;
    this.viewer.MoleculeColumn = this.columnName!;
    this.viewer.dataFrame = dataFrame;
    this.root.appendChild(this.viewer.root);
  }
}
