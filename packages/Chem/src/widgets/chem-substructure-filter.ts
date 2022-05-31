/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {chem} from 'datagrok-api/grok';
import {chemSubstructureSearchLibrary} from '../chem-searches';
import {initRdKitService} from '../utils/chem-common-rdkit';
import {Subscription} from 'rxjs';
import {debounceTime, filter} from 'rxjs/operators';
import Sketcher = chem.Sketcher;
import wu from 'wu';

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
    return super.isFiltering && !this.sketcher?.getMolFile()?.endsWith(this.WHITE_MOL);
  }

  constructor() {
    super();
    initRdKitService(); // No await
    this.root = ui.divV([]);
    this._indicateProgress(false);
    this.loader.style.position = 'absolute';
    this.loader.style.right = '60px';
    this.loader.style.top = '4px';
    this.root.appendChild(this.sketcher.root);
    this.root.firstChild!.firstChild!.firstChild!.appendChild(this.loader);
  }

  get _debounceTime(): number {
    if (this.column == null)
      return 1000;
    const sz = this.column.length;
    const szMin = 500;
    const szMax = 10000;
    const msecMax = 1000;
    if (sz < szMin) return 0;
    if (sz > szMax) return msecMax;
    return Math.floor(msecMax * ((sz - szMin) / (szMax - szMin)));
  }

  attach(dataFrame: DG.DataFrame): void {
    super.attach(dataFrame);

    this.column = dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.columnName = this.column?.name;
    this.onSketcherChangedSubs?.unsubscribe();

    // hide the scaffold when user deactivates the filter
    this.subs.push(this.dataFrame!.onRowsFiltering
      .pipe(filter((_) => this.column != null && !this.isFiltering))
      .subscribe((_) => delete this.column!.temp['chem-scaffold-filter']));

    chemSubstructureSearchLibrary(this.column!, '', '')
      .then((_) => {}); // Nothing, just a warmup

    let onChangedEvent: any = this.sketcher.onChanged;
    onChangedEvent = onChangedEvent.pipe(debounceTime(this._debounceTime));
    this.onSketcherChangedSubs = onChangedEvent.subscribe(async (_: any) => await this._onSketchChanged());

    if (this.column?.temp['chem-scaffold-filter'])
      this.sketcher.setMolFile(this.column?.temp['chem-scaffold-filter']);
  }

  detach() {
    super.detach();
    if (this.column?.temp['chem-scaffold-filter'])
      this.column.temp['chem-scaffold-filter'] = null;
  }

  applyFilter(): void {
    if (this.bitset && !this.isDetached) {
      this.dataFrame?.filter.and(this.bitset);
      this.column!.temp['chem-scaffold-filter'] = this.sketcher.getMolFile();
    }
  }

  /** Override to save filter state. */
  saveState(): any {
    const state = super.saveState();
    state.molBlock = this.sketcher.getMolFile();
    return state;
  }

  /** Override to load filter state. */
  applyState(state: any): void {
    super.applyState(state);
    if (state.molBlock)
      this.sketcher.setMolFile(state.molBlock);
    //
    // if (state.molBlock)
    //   setTimeout(() => this._onSketchChanged(), 1000);
  }

  async _onSketchChanged(): Promise<void> {
    if (!this.isFiltering) {
      if (this.column?.temp['chem-scaffold-filter'])
        delete this.column.temp['chem-scaffold-filter'];
      this.bitset = null;
    }
    else if (wu(this.dataFrame!.rows.filters).has(`${this.columnName}: ${this.filterSummary}`)) {
      // some other filter is already filtering for the exact same thing
      return;
    }
    else {
      this._indicateProgress();
      try {
        this.bitset = await chemSubstructureSearchLibrary(
          this.column!, this.sketcher.getMolFile(), await this.sketcher.getSmarts());
        this.dataFrame?.rows.requestFilter();
      } finally {
        this._indicateProgress(false);
      }
    }
  }
}
