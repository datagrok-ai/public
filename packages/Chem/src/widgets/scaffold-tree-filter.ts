import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ScaffoldTreeViewer} from "./scaffold-tree";
import {filter, debounce} from 'rxjs/operators';
import {interval} from 'rxjs';

export class ScaffoldTreeFilter extends DG.Filter {
  viewer: ScaffoldTreeViewer = new ScaffoldTreeViewer();
  savedTree: string = '';
  
  constructor() {
    super();
    this.root = ui.divV([]);
    this.subs = this.viewer.subs;
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
    this.createViewer(dataFrame);
    this.subs.push(this.dataFrame!.onRowsFiltering
      .pipe(filter((_) => !this.isFiltering), debounce(_ => interval(100)))
      .subscribe((_) => {
        delete this.column!.temp['chem-scaffold-filter'];
        this.viewer.updateFilters(this.isFiltering);
      })
    );
  }

  saveState(): any {
    const state = super.saveState();
    state.savedTree = JSON.stringify(ScaffoldTreeViewer.serializeTrees(this.viewer.tree));
    return state;
  }

  applyState(state: any): void {
    super.applyState(state);
    this.viewer.loadTreeStr(state.savedTree);
  }

  applyFilter(): void {
  }

  detach(): void {
    super.detach();
    this.viewer.clearFilters();
  }

  createViewer(dataFrame: DG.DataFrame) {
    this.viewer.allowGenerate = false;
    this.viewer.size = 'small';
    this.viewer.dataFrame = dataFrame;
    this.root.appendChild(this.viewer.root);
  }
}