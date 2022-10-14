/**
 * Macromolecules substructure filter that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import wu from 'wu';
import {linearSubstructureSearch} from '../substructure-search/substructure-search';
import {Subject, Subscription} from 'rxjs';
import * as C from '../utils/constants';

export class BioSubstructureFilter extends DG.Filter {
  bioFilter: FastaFilter | SeparatorFilter | null = null;
  bitset: DG.BitSet | null = null;
  loader: HTMLDivElement = ui.loader();
  onBioFilterChangedSubs?: Subscription;

  get calculating(): boolean { return this.loader.style.display == 'initial'; }

  set calculating(value: boolean) { this.loader.style.display = value ? 'initial' : 'none'; }

  get filterSummary(): string {
    return this.bioFilter!.substructure;
  }

  get isFiltering(): boolean {
    return super.isFiltering && this.bioFilter!.substructure !== '';
  }

  get isReadyToApplyFilter(): boolean {
    return !this.calculating && this.bitset != null;
  }

  constructor() {
    super();
    this.root = ui.divV([]);
    this.calculating = false;
  }

  attach(dataFrame: DG.DataFrame): void {
    super.attach(dataFrame);
    this.column = dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    this.columnName = this.column?.name;
    const notation = this.column?.getTag(DG.TAGS.UNITS);
    this.bioFilter = notation === bio.NOTATION.FASTA ?
      new FastaFilter() : new SeparatorFilter(this.column!.getTag(C.TAGS.SEPARATOR));
    this.root.appendChild(this.bioFilter!.filterPanel);
    this.root.appendChild(this.loader);

    this.onBioFilterChangedSubs?.unsubscribe();
    const onChangedEvent: any = this.bioFilter.onChanged;
    this.onBioFilterChangedSubs = onChangedEvent.subscribe(async (_: any) => await this._onInputChanged());
  }

  detach() {
    super.detach();
  }

  applyFilter(): void {
    if (this.bitset && !this.isDetached)
      this.dataFrame?.filter.and(this.bitset);
  }

  /** Override to save filter state. */
  saveState(): any {
    const state = super.saveState();
    state.bioSubstructure = this.bioFilter?.substructure;
    return state;
  }

  /** Override to load filter state. */
  applyState(state: any): void {
    super.applyState(state);
    if (state.bioSubstructure)
      this.bioFilter!.substructureInput.value = state.bioSubstructure;

    const that = this;
    if (state.bioSubstructure)
      setTimeout(function() { that._onInputChanged(); }, 1000);
  }

  /**
   * Performs the actual filtering
   * When the results are ready, triggers `rows.requestFilter`, which in turn triggers `applyFilter`
   * that would simply apply the bitset synchronously.
   */
  async _onInputChanged(): Promise<void> {
    if (!this.isFiltering) {
      this.bitset = null;
      this.dataFrame?.rows.requestFilter();
    } else if (wu(this.dataFrame!.rows.filters).has(`${this.columnName}: ${this.filterSummary}`)) {
      // some other filter is already filtering for the exact same thing
      return;
    } else {
      this.calculating = true;
      try {
        this.bitset = linearSubstructureSearch(this.bioFilter!.substructure, this.column!);
        this.calculating = false;
        this.dataFrame?.rows.requestFilter();
      } finally {
        this.calculating = false;
      }
    }
  }
}

class FastaFilter {
  onChanged: Subject<any> = new Subject<any>();
  substructureInput: DG.InputBase<string> = ui.stringInput('', '', () => {
    this.onChanged.next();
  }, {placeholder: 'Substructure'});

  constructor() {
  }

  get filterPanel() {
    return this.substructureInput.root;
  }

  get substructure() {
    return this.substructureInput.value;
  }
}

class SeparatorFilter extends FastaFilter {
  separatorInput: DG.InputBase<string> = ui.stringInput('', '', () => {
    this.onChanged.next();
  }, {placeholder: 'Separator'});
  colSeparator = '';

  constructor(separator: string) {
    super();
    this.colSeparator = separator;
    this.separatorInput.value = separator;
  }

  get filterPanel() {
    return ui.divV([
      this.substructureInput.root,
      this.separatorInput.root
    ]);
  }

  get substructure() {
    return this.separatorInput.value && this.separatorInput.value !== this.colSeparator ?
      this.substructureInput.value.replaceAll(this.separatorInput.value, this.colSeparator) :
      this.substructureInput.value;
  }
}
