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
import {TaskBarProgressIndicator, chem} from 'datagrok-api/dg';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../package';
import {MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT, getSearchProgressEventName, getTerminateEventName} from '../constants';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

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
  errorDiv = ui.divText(`Too many rows, maximum for substructure search is ${MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT}`,
    'chem-substructure-limit');
  sketcherType = DG.chem.currentSketcherType;
  currentSearches = new Set<string>();
  progressBar: TaskBarProgressIndicator | null = null;
  batchResultObservable: Subscription | null = null;
  terminateEventName: string = '';
  progressEventName: string = '';
  currentMolfile: string = '';
  initListeners = false;

  get calculating(): boolean {return this.loader.style.display == 'initial';}
  set calculating(value: boolean) {this.loader.style.display = value ? 'initial' : 'none';}

  get filterSummary(): string {
    return _convertMolNotation(this.currentMolfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts,
      getRdKitModule());
  }

  get isFiltering(): boolean {
    const molFile = this.sketcher?.getMolFile();
    return super.isFiltering && (!!molFile && !chem.Sketcher.isEmptyMolfile(molFile));
  }

  get isReadyToApplyFilter(): boolean {
    return this.bitset != null;
  }

  constructor() {
    super();
    initRdKitService(); // No await
    this.filterId = id++;
    this.root = ui.divV([]);
    this.calculating = false;
    this.root.appendChild(this.sketcher.root);
    this.root.appendChild(this.loader);
    this.subs.push(grok.events.onResetFilterRequest.subscribe((_) => {
      {
        this.sketcher.setMolFile(DG.WHITE_MOLBLOCK);
        this.updateExternalSketcher();
      }
    }));
    this.subs.push(grok.events.onCustomEvent(FILTER_SYNC_EVENT).subscribe((state: ISubstructureFilterState) => {
      if (state.colName === this.columnName && this.tableName == state.tableName && this.filterId !== state.filterId) {
        /* setting syncEvent to true only if base sketcher is initialized.
        If base sketcher is initialized, it will fire onChange event */
        if (this.currentSearches.size > 0) 
            grok.events.fireCustomEvent(this.terminateEventName, this.currentSearches.values().next().value);
        if (this.sketcher.sketcher?.isInitialized)
          this.syncEvent = true;
        this.bitset = state.bitset!;
        this.sketcher.setMolFile(state.molblock!);
        this.updateExternalSketcher();
      }
    }));
    this.subs.push(grok.events.onCustomEvent(SKETCHER_TYPE_CHANGED).subscribe((state: ISubstructureFilterState) => {
      if (state.colName === this.columnName && this.tableName == state.tableName && this.filterId !== state.filterId) {
        if (this.sketcher.sketcher?.isInitialized) {
          if (DG.chem.currentSketcherType !== this.sketcherType &&
              this.sketcher._mode !== DG.chem.SKETCHER_MODE.EXTERNAL) {
            this.sketcherType = DG.chem.currentSketcherType;
            this.sketcher.sketcherType = DG.chem.currentSketcherType;
          }
        }
      }
    }));
  }

  get _debounceTime(): number {
    return 0;
    // if (this.column == null)
    //   return 1000;
    // const length = this.column.length;
    // const minLength = 500;
    // const maxLength = 10000;
    // const msecMax = 1000;
    // if (length < minLength) return 0;
    // if (length > maxLength) return msecMax;
    // return Math.floor(msecMax * ((length - minLength) / (maxLength - minLength)));
  }

  attach(dataFrame: DG.DataFrame): void {
    if (dataFrame.rowCount > MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT) {
      ui.tools.waitForElementInDom(this.sketcher.root).then(() => {
        this.sketcher.root.children[0].classList.add('chem-hide-filter');
        this.sketcher.root.append(this.errorDiv);
        return;
      });
    }
    super.attach(dataFrame);
    this.column ??=dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.columnName ??= this.column?.name ?? '';
    this.tableName = dataFrame.name ?? '';
    this.onSketcherChangedSubs?.unsubscribe();

    // hide the scaffold when user deactivates the filter
    this.subs.push(this.dataFrame!.onRowsFiltering
      .pipe(filter((_) => this.column != null && !this.isFiltering))
      .subscribe((_: any) => delete this.column!.temp['chem-scaffold-filter']));

    this.currentSearches.add('');
    chemSubstructureSearchLibrary(this.column!, '', '', false, false)
      .then((_) => {}); // Precalculating fingerprints

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
    this.batchResultObservable?.unsubscribe();
  }

  applyFilter(): void {
   // console.log(`in apply filter ${this.sketcher.getMolFile()}`)
      if (this.dataFrame && this.bitset && !this.isDetached) {
        this.dataFrame.filter.and(this.bitset);
        this.dataFrame.rows.addFilterState(this.saveState());
        this.column!.temp['chem-scaffold-filter'] = this.currentMolfile;
        this.active = true;
    }
  }

  /** Override to save filter state. */
  saveState(): any {
    const state = super.saveState();
    state.type = 'Chem:substructureFilter';
    state.molBlock = this.currentMolfile;
    return state;
  }

  /** Override to load filter state. */
  applyState(state: any): void {
    super.applyState(state);
    if (!this.initListeners) {
      this.initListeners = true;
      this.terminateEventName = getTerminateEventName(this.tableName, this.columnName!);
      this.progressEventName = getSearchProgressEventName(this.tableName, this.columnName!);

      this.subs.push(grok.events.onCustomEvent(this.terminateEventName).subscribe((queryMol: string) => {
        this.finishSearch(queryMol);
    }));
    }
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
    const newMolFile = this.sketcher.getMolFile();
    const newSmarts = await this.sketcher.getSmarts();
    grok.events.fireCustomEvent(SKETCHER_TYPE_CHANGED, {colName: this.columnName,
      filterId: this.filterId, tableName: this.tableName});
    if (!this.isFiltering) {
      this.currentMolfile = newMolFile; 
      this.bitset = !this.active ? DG.BitSet.fromBytes((await this.getFilterBitset())!.buffer.buffer, this.column!.length) : null;//TODO
      if (this.column?.temp['chem-scaffold-filter'])
        delete this.column.temp['chem-scaffold-filter'];
      this.dataFrame?.rows.requestFilter();
      this.terminatePreviousSearch();
      this.finishSearch(newSmarts ?? '');
      grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
        molblock: this.currentMolfile, colName: this.columnName,
        filterId: this.filterId, tableName: this.tableName}); 
    } else if (wu(this.dataFrame!.rows.filters)
        .has(`${this.columnName}: ${_convertMolNotation(newMolFile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts,
      getRdKitModule())}`)) {
      // some other filter is already filtering for the exact same thing
      return;
    } else {
        this.terminatePreviousSearch();
        this.currentMolfile = newMolFile; 
        this.currentSearches.add(newSmarts ?? '');
        this.calculating = true;
        this.progressBar ??= DG.TaskBarProgressIndicator.create(`Starting substructure search...`);
        try {
          grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
            molblock: this.currentMolfile, colName: this.columnName,
            filterId: this.filterId, tableName: this.tableName});  

          const bitArray = await this.getFilterBitset();
          this.bitset = DG.BitSet.fromBytes(bitArray.buffer.buffer, this.column!.length);
          this.batchResultObservable?.unsubscribe();
          this.batchResultObservable = grok.events.onCustomEvent(this.progressEventName).subscribe((progress: number) => {

            this.bitset = DG.BitSet.fromBytes(bitArray.buffer.buffer, this.column!.length);
            this.dataFrame?.rows.requestFilter();
            this.progressBar!.update(progress, `${progress?.toFixed(2)}% of search completed`);
            
        })
        } catch {
          this.finishSearch(newSmarts ?? '');
        }
    }
  }

  async getFilterBitset(): Promise<BitArray> {
    console.log(`getFilterBitset currentSearches: ${this.currentSearches}`);
    const smarts = await this.sketcher.getSmarts();
    return await chemSubstructureSearchLibrary(this.column!, this.currentMolfile, smarts!, false, false);
  }

  updateExternalSketcher() {
    if (this.sketcher._mode === DG.chem.SKETCHER_MODE.EXTERNAL)
      this.sketcher.updateExtSketcherContent(); //updating image in minimized sketcher panel
  }

  terminatePreviousSearch() {
    if (this.currentSearches.size > 0) 
      grok.events.fireCustomEvent(this.terminateEventName, this.currentSearches.values().next().value);
  }

  finishSearch(queryMol: string) {
    const finish = () => {
      if (this.currentSearches.size === 0) {
        this.calculating = false;
        this.progressBar?.close();
        this.progressBar = null;
        this.batchResultObservable?.unsubscribe();
        console.log(`Unsubscribed from batchResultObservable  Filter ${this.filterId}`);
      }
    }
    if (this.currentSearches.has(queryMol)) {
      this.currentSearches.delete(queryMol);
      finish();
    }
    if (queryMol == null && this.currentSearches.size === 1) {
      const v = this.currentSearches.values().next().value;
      if (v != null)
        this.currentSearches.delete(v);
      finish();
    }
  }

}
