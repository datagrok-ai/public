/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {chemSubstructureSearchLibrary} from '../chem-searches';
import {initRdKitService} from '../utils/chem-common-rdkit';
import {Subscription} from 'rxjs';
import {debounceTime, filter} from 'rxjs/operators';
import wu from 'wu';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';
import { chem } from 'datagrok-api/dg';
import { _convertMolNotation } from '../utils/convert-notation-utils';
import { getRdKitModule } from '../package';

const FILTER_SYNC_EVENT = 'chem-substructure-filter';
const SKETCHER_TYPE_CHANGED = 'chem-sketcher-type-changed';
let id = 0;

interface ISubstructureFilterState {
  bitset?: DG.BitSet;
  molblock?: string;
  colName: string;
  filterId: number;
  tableName: string
}

export class SubstructureFilter extends DG.Filter {
  // @ts-ignore
  sketcher: DG.chem.Sketcher = new DG.chem.Sketcher();
  bitset: DG.BitSet | null = null;
  loader: HTMLDivElement = ui.loader();
  onSketcherChangedSubs?: Subscription;
  active: boolean = true;
  syncEvent = false;
  filterId: number;
  tableName: string = '';

  get calculating(): boolean {return this.loader.style.display == 'initial';}
  set calculating(value: boolean) {this.loader.style.display = value ? 'initial' : 'none';}

  get filterSummary(): string {
    return _convertMolNotation(this.sketcher.getMolFile(), DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts, getRdKitModule());
  }

  get isFiltering(): boolean {
    const molFile = this.sketcher?.getMolFile();
    return super.isFiltering && (!!molFile && !chem.Sketcher.isEmptyMolfile(molFile));
  }

  get isReadyToApplyFilter(): boolean {
    return !this.calculating && this.bitset != null;
  }

  constructor() {
    super();
    initRdKitService(); // No await
    this.filterId = id++;
    this.root = ui.divV([]);
    this.calculating = false;
    this.root.appendChild(this.sketcher.root);
    this.root.appendChild(this.loader);
    this.subs.push(grok.events.onResetFilterRequest.subscribe((_) => { {
      this.sketcher.setMolFile(DG.WHITE_MOLBLOCK);
      this.updateExternalSketcher();
    } }));
    this.subs.push(grok.events.onCustomEvent(FILTER_SYNC_EVENT).subscribe((state: ISubstructureFilterState) => {
      if (state.colName === this.columnName && this.tableName == state.tableName && this.filterId !== state.filterId) {
        if (this.sketcher.sketcher?.isInitialized) //setting syncEvent to true only if base sketcher is initialized. If base sketcher is initialized, it will fire onChange event
          this.syncEvent = true;
        this.bitset = state.bitset!;
        this.sketcher.setMolFile(state.molblock!);
        this.updateExternalSketcher();
      }
    }));
    this.subs.push(grok.events.onCustomEvent(SKETCHER_TYPE_CHANGED).subscribe((state: ISubstructureFilterState) => {
      if (state.colName === this.columnName && this.tableName == state.tableName && this.filterId !== state.filterId) {
        if (this.sketcher.sketcher?.isInitialized) {
          if (DG.chem.currentSketcherType !== this.sketcher.sketcher!.name) {
            this.sketcher.sketcherType = DG.chem.currentSketcherType;
          }
        }
      }
    }))
  }

  get _debounceTime(): number {
    if (this.column == null)
      return 1000;
    const length = this.column.length;
    const minLength = 500;
    const maxLength = 10000;
    const msecMax = 1000;
    if (length < minLength) return 0;
    if (length > maxLength) return msecMax;
    return Math.floor(msecMax * ((length - minLength) / (maxLength - minLength)));
  }

  attach(dataFrame: DG.DataFrame): void {
    super.attach(dataFrame);

    this.column ??= dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.columnName ??= this.column?.name;
    this.tableName = dataFrame.name;
    this.onSketcherChangedSubs?.unsubscribe();

    // hide the scaffold when user deactivates the filter
    this.subs.push(this.dataFrame!.onRowsFiltering
      .pipe(filter((_) => this.column != null && !this.isFiltering))
      .subscribe((_: any) => delete this.column!.temp['chem-scaffold-filter']));

    chemSubstructureSearchLibrary(this.column!, '', '')
      .then((_) => {}); // Nothing, just a warmup

    let onChangedEvent: any = this.sketcher.onChanged;
    onChangedEvent = onChangedEvent.pipe(debounceTime(this._debounceTime));
    this.onSketcherChangedSubs = onChangedEvent.subscribe(async (_: any) => {
      this.syncEvent === true ? this.syncEvent = false : await this._onSketchChanged();
    });

  }

  refresh() {
    if (!this.sketcher.sketcherTypeChanged)
      this.sketcher.sketcher?.refresh();
  }

  detach() {
    super.detach();
    if (this.column?.temp['chem-scaffold-filter'])
      this.column.temp['chem-scaffold-filter'] = null;
  }

  applyFilter(): void {
    if (this.bitset && !this.isDetached) {
      this.dataFrame?.filter.and(this.bitset);
      this.dataFrame?.rows.addFilterState(this.saveState());
      this.column!.temp['chem-scaffold-filter'] = this.sketcher.getMolFile();
      this.active = true;
    }
  }

  /** Override to save filter state. */
  saveState(): any {
    const state = super.saveState();
    state.type = 'Chem:substructureFilter';
    state.molBlock = this.sketcher.getMolFile();
    return state;
  }

  /** Override to load filter state. */
  applyState(state: any): void {
      super.applyState(state);
      this.active = state.active ?? true;
      if (this.column?.temp['chem-scaffold-filter'])
        state.molBlock = this.column?.temp['chem-scaffold-filter'];
      if (state.molBlock) {
        this.sketcher.setMolFile(state.molBlock);
        this.updateExternalSketcher();
      }

      const that = this;
      if (state.molBlock)
        setTimeout(function() {that._onSketchChanged();}, 1000);
  }

  /**
   * Performs the actual filtering
   * When the results are ready, triggers `rows.requestFilter`, which in turn triggers `applyFilter`
   * that would simply apply the bitset synchronously.
   */
  async _onSketchChanged(): Promise<void> {
    grok.events.fireCustomEvent(SKETCHER_TYPE_CHANGED, {colName: this.columnName,
      filterId: this.filterId, tableName: this.tableName});
    if (!this.isFiltering) {
      this.bitset = !this.active ? await this.getFilterBitset() : null;
      if (this.column?.temp['chem-scaffold-filter'])
        delete this.column.temp['chem-scaffold-filter'];
      grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset, colName: this.columnName,
        molblock: this.sketcher.getMolFile(), filterId: this.filterId, tableName: this.tableName});
      this.dataFrame?.rows.requestFilter();
    } else if (wu(this.dataFrame!.rows.filters).has(`${this.columnName}: ${this.filterSummary}`)) {
      // some other filter is already filtering for the exact same thing
      return;
    } else {
      this.calculating = true;
      try {
        const bitset = await this.getFilterBitset();
        if (!bitset)
          return;
        this.bitset = bitset;
        this.calculating = false;
        grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
          molblock: this.sketcher.getMolFile(), colName: this.columnName, filterId: this.filterId, tableName: this.tableName});
        this.dataFrame?.rows.requestFilter();
      } finally {
        this.calculating = false;
      }
    }
  }

  async getFilterBitset(): Promise<DG.BitSet | null> {
    const smarts = await this.sketcher.getSmarts();
        if (StringUtils.isEmpty(smarts) && StringUtils.isEmpty(this.sketcher.getMolFile()))
          return null;
    return await chemSubstructureSearchLibrary(this.column!, this.sketcher.getMolFile(), smarts!);
  }

  updateExternalSketcher() {
    if (this.sketcher._mode === DG.chem.SKETCHER_MODE.EXTERNAL)
      this.sketcher.updateExtSketcherContent(); //updating image in minimized sketcher panel
  }

}
