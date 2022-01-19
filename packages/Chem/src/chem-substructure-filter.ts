/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {chem} from 'datagrok-api/grok';
import {chemSubstructureSearchLibrary} from "./chem-searches";
import {initRdKitService} from "./chem-common-rdkit";
import {Observable, Subscription} from "rxjs";
import {debounceTime} from 'rxjs/operators';
import Sketcher = chem.Sketcher;

export class SubstructureFilter extends DG.Filter {
  sketcher: Sketcher = new Sketcher();
  bitset: DG.BitSet | null = null;
  loader = ui.loader();
  readonly WHITE_MOL = '\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n';
  onSketcherChangedSubs?: Subscription;

  _indicateProgress(on = true) {
    this.loader.style.display = on ? 'block' : 'none';
  }

  get filterSummary(): string {
    return this.sketcher.getSmiles();
  }

  get isFiltering(): boolean {
    return !this.sketcher.getMolFile().endsWith(this.WHITE_MOL);
  }

  constructor() {
    super();
    initRdKitService();  // No await
    this.root = ui.divV([]);
    this._indicateProgress(false);
    this.loader.style.position = 'absolute';
    this.loader.style.right = '60px';
    this.loader.style.top = '4px';
    this.root.appendChild(this.sketcher.root);
    this.root.firstChild!.firstChild!.firstChild!.appendChild(this.loader);
  }

  get _debounceTime(): number {
    const sz = this.column!.length;
    const szMin = 500, szMax = 10000, msecMax = 1000;
    if (sz < szMin) { return 0; }
    if (sz > szMax) { return msecMax; }
    return Math.floor(msecMax * ((sz - szMin) / (szMax - szMin)));
  }

  attach(dataFrame: DG.DataFrame): void {
    if (this.onSketcherChangedSubs)
      this.onSketcherChangedSubs.unsubscribe();
    this.column = dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    chemSubstructureSearchLibrary(this.column!, '', ''); // No await
    this.columnName = this.column?.name;
    let onChangedEvent: Observable<any> = this.sketcher.onChanged;
    const t = this._debounceTime;
    if (t > 0)
      onChangedEvent = onChangedEvent.pipe(debounceTime(t));
    this.onSketcherChangedSubs = onChangedEvent
      .subscribe(async _ => await this._onSketchChanged());
    super.attach(dataFrame);
    if (this.column?.temp['chem-scaffold-filter'])
      this.sketcher.setMolFile(this.column?.temp['chem-scaffold-filter']);
  }

  applyFilter(): void {
    if (this.bitset) {
      this.dataFrame?.filter.and(this.bitset);
      this.column!.temp['chem-scaffold-filter'] = this.sketcher.getMolFile();
    }
  }

  async _onSketchChanged(): Promise<void> {
    if (!this.isFiltering) {
      if (this.column?.temp['chem-scaffold-filter'])
        delete this.column.temp['chem-scaffold-filter'];
      this.bitset = null;
    } else {
      this._indicateProgress();
      try {
        this.bitset = await chemSubstructureSearchLibrary(
          this.column!, this.sketcher.getMolFile(), await this.sketcher.getSmarts());
      }
      finally {
        this._indicateProgress(false);
      }
    }
    this.dataFrame?.rows.requestFilter();
  }
}
