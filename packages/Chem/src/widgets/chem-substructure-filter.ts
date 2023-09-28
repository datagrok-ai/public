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
import {Subject, Subscription} from 'rxjs';
import {debounceTime, filter} from 'rxjs/operators';
import wu from 'wu';
import {TaskBarProgressIndicator, chem} from 'datagrok-api/dg';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../package';
import {AVAILABLE_FPS, FILTER_SCAFFOLD_TAG, MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT, SubstructureSearchType, getSearchProgressEventName, getTerminateEventName} from '../constants';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { IColoredScaffold } from '../rendering/rdkit-cell-renderer';
import { Fingerprint } from '../utils/chem-common';

const FILTER_SYNC_EVENT = 'chem-substructure-filter';
const SKETCHER_TYPE_CHANGED = 'chem-sketcher-type-changed';
let id = 0;

interface ISubstructureFilterState {
  bitset?: DG.BitSet;
  molblock?: string;
  colName: string;
  filterId: number;
  tableName: string;
  searchType: SubstructureSearchType;
  simCutOff: number;
  fp: Fingerprint
}

export class SubstructureFilter extends DG.Filter {
  // @ts-ignore
  sketcher: DG.chem.Sketcher = new DG.chem.Sketcher();
  bitset: DG.BitSet | null = null;
  loader: HTMLDivElement = ui.loader();
  onSketcherChangedSubs?: Subscription[] = [];
  active: boolean = true;
  syncEvent = false;
  requestedFilter = false; // required for filters sync
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
  //searchTypeLink: HTMLButtonElement;
  searchType: SubstructureSearchType = SubstructureSearchType.CONTAINS;
  similarityCutOff = 0.8;
  fp: Fingerprint = Fingerprint.Morgan;
  searchTypes = [SubstructureSearchType.CONTAINS, SubstructureSearchType.INCLUDED_IN,
    SubstructureSearchType.EXACT_MATCH, SubstructureSearchType.IS_SIMILAR];
  fpsTypes = AVAILABLE_FPS;
  searchTypeInput: DG.InputBase;
  similarityCutOffInput: DG.InputBase;
  fpInput: DG.InputBase;
  simOptionsDiv = ui.divH([], {style: {marginBottom: '3px', gap: '5px'}});
  searchTypeChanged = new Subject();
  
  get calculating(): boolean {return this.loader.style.display == 'initial';}
  set calculating(value: boolean) {this.loader.style.display = value ? 'initial' : 'none';}

  get filterSummary(): string {
    const smarts = _convertMolNotation(this.currentMolfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts,
      getRdKitModule());
    return `${smarts}_${this.searchType}_${this.similarityCutOff}_${this.fp}`;
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

    this.searchTypeInput = ui.choiceInput('', this.searchType, this.searchTypes, () => { 
      this.searchType = this.searchTypeInput.value;
      if (this.searchType === SubstructureSearchType.IS_SIMILAR) {
        this.simOptionsDiv.append(this.fpInput.root);
        this.simOptionsDiv.append(this.similarityCutOffInput.root);
      } else {
        ui.empty(this.simOptionsDiv);
      }
      if (!this.syncEvent)
        this.searchTypeChanged.next();
    });

    this.fpInput = ui.choiceInput('FP', this.fp, this.fpsTypes, () => {
      this.fp = this.fpInput.value;
        if (!this.syncEvent)
          this.searchTypeChanged.next();      
    });
    this.fpInput.input.classList.add('chem-filter-fp-editor');
    this.fpInput.captionLabel.classList.add('chem-filter-fp-label');

    const property =
    {
      "name": "lim",
      "type": DG.TYPE.FLOAT,
      "showSlider": true,
      "min": 0,
      "max": 1,
      "nullable": false,
    };
    const slider = DG.Property.fromOptions(property);
    const initialCutOff = { lim: 0.8 };
    this.similarityCutOffInput = ui.input.forProperty(slider, initialCutOff);
    this.similarityCutOffInput.input.classList.add('chem-filter-sim-cutoff-editor');
    this.similarityCutOffInput.captionLabel.classList.add('chem-filter-sim-cutoff-label');
    this.similarityCutOffInput.onChanged(() => {
      this.similarityCutOff = this.similarityCutOffInput.value;
      if (!this.syncEvent)
        this.searchTypeChanged.next();  
    });

    ui.tooltip.bind((this.similarityCutOffInput.root.getElementsByClassName('ui-input-slider')[0] as HTMLElement)!, 'Similarity cutoff');

    this.root.appendChild(ui.div(this.searchTypeInput.root));
    this.root.appendChild(this.simOptionsDiv);
    this.root.appendChild(ui.div(this.sketcher.root, {style: {position: 'relative'}}));
    this.root.appendChild(this.loader);
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
    if (dataFrame.rowCount > MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT) {
      ui.tools.waitForElementInDom(this.sketcher.root).then(() => {
        this.sketcher.root.children[0].classList.add('chem-hide-filter');
        this.sketcher.root.append(this.errorDiv);
        return;
      });
    }
    super.attach(dataFrame);
    this.column ??= dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.columnName ??= this.column?.name ?? '';
    this.tableName = dataFrame.name ?? '';
    this.onSketcherChangedSubs?.forEach((it) => it.unsubscribe());

    // hide the scaffold when user deactivates the filter
    this.subs.push(this.dataFrame!.onRowsFiltering
      .pipe(filter((_) => this.column != null && !this.isFiltering))
      .subscribe((_: any) => delete this.column!.temp[FILTER_SCAFFOLD_TAG]));

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
        this.currentMolfile = state.molblock!;
        this.bitset = state.bitset!;
        this.searchTypeInput.value = state.searchType;
        this.similarityCutOffInput.value = state.simCutOff;
        this.fpInput.value = state.fp;
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

    this.currentSearches.add('');
    chemSubstructureSearchLibrary(this.column!, '', '', false, false)
      .then((_) => { }); // Precalculating fingerprints

    let onChangedEvent: any = this.sketcher.onChanged;
    //onChangedEvent = onChangedEvent.pipe(debounceTime(this._debounceTime));
    this.onSketcherChangedSubs?.push(onChangedEvent.subscribe(async (_: any) => {
      this.syncEvent === true ? this.syncEvent = false : await this._onSketchChanged();
    }));
    let searchTypeChanged = this.searchTypeChanged.pipe(debounceTime(this._debounceTime));
    this.onSketcherChangedSubs?.push(searchTypeChanged.subscribe(async (_: any) => {
      await this._onSketchChanged();
    }));
  }

  refresh() {
    if (!this.sketcher.sketcherTypeChanged)
      this.sketcher.sketcher?.refresh();
  }

  detach() {
    this.sketcher.getSmarts().then((smarts) => {
      this.terminatePreviousSearch();
      this.finishSearch(smarts ?? '');
      this.onSketcherChangedSubs?.forEach((it) => it.unsubscribe());
    });
    super.detach();
    if (this.column?.temp[FILTER_SCAFFOLD_TAG])
      this.column.temp[FILTER_SCAFFOLD_TAG] = null;
    this.batchResultObservable?.unsubscribe();
  }

  applyFilter(): void {
    // console.log(`in apply filter ${this.sketcher.getMolFile()}`)
    if (this.dataFrame && this.bitset && !this.isDetached && this.requestedFilter) {
      this.dataFrame.filter.and(this.bitset);
      this.dataFrame.rows.addFilterState(this.saveState());
      this.column!.temp[FILTER_SCAFFOLD_TAG] = JSON.stringify([{
        molecule: this.currentMolfile,
        isSuperstructure: this.searchType === SubstructureSearchType.INCLUDED_IN
      }]);
      this.active = true;
      this.requestedFilter = false;
    }
  }

  /** Override to save filter state. */
  saveState(): any {
    const state = super.saveState();
    state.type = 'Chem:substructureFilter';
    state.molBlock = this.currentMolfile;
    state.searchType = this.searchType;
    state.simCutOff = this.similarityCutOff;
    state.fp = this.fp;
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
    if (this.column?.temp[FILTER_SCAFFOLD_TAG])
      state.molBlock = (JSON.parse(this.column?.temp[FILTER_SCAFFOLD_TAG]) as IColoredScaffold[])[0].molecule;
    if (state.molBlock) {
      this.currentMolfile = state.molBlock;
      this.sketcher.setMolFile(state.molBlock);
      this.updateExternalSketcher();
    }
    if (state.searchType)
      this.searchTypeInput.value = state.searchType;
    if (state.simCutOff)
      this.similarityCutOffInput.value = state.simCutOff;
    if (state.fp)
      this.fpInput.value = state.fp;

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
      this.bitset = !this.active ?
        DG.BitSet.fromBytes((await this.getFilterBitset())!.buffer.buffer, this.column!.length) : null;//TODO
      if (this.column?.temp[FILTER_SCAFFOLD_TAG])
        delete this.column.temp[FILTER_SCAFFOLD_TAG];
      this.terminatePreviousSearch();
      this.finishSearch(newSmarts ?? '');
      grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
        molblock: this.currentMolfile, colName: this.columnName, filterId: this.filterId, 
        tableName: this.tableName, searchType: this.searchType, simCutOff: this.similarityCutOff, fp: this.fp});
      this.requestedFilter = true;
      this.dataFrame?.rows.requestFilter();
    } else if (wu(this.dataFrame!.rows.filters)
      .has(`${this.columnName}: ${_convertMolNotation(newMolFile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts,
        getRdKitModule())}_${this.searchType}_${this.similarityCutOff}_${this.fp}`)) {
      // some other filter is already filtering for the exact same thing
      // value to pass into has() is created similarly to filterSummary property 
      return;
    } else {
      this.terminatePreviousSearch();
      this.currentMolfile = newMolFile;
      this.currentSearches.add(newSmarts ?? '');
      this.calculating = true;
      this.progressBar ??= DG.TaskBarProgressIndicator.create(`Starting substructure search...`);
      try {
        grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
          molblock: this.currentMolfile, colName: this.columnName, filterId: this.filterId, 
          tableName: this.tableName, searchType: this.searchType, simCutOff: this.similarityCutOff, fp: this.fp});

        const bitArray = await this.getFilterBitset();
        this.bitset = DG.BitSet.fromBytes(bitArray.buffer.buffer, this.column!.length);
        this.batchResultObservable?.unsubscribe();
        this.batchResultObservable = grok.events.onCustomEvent(this.progressEventName).subscribe((progress: number) => {
          this.bitset = DG.BitSet.fromBytes(bitArray.buffer.buffer, this.column!.length);
          this.requestedFilter = true;
          this.dataFrame?.rows.requestFilter();
            this.progressBar?.update(progress, `${progress?.toFixed(2)}% of search completed`);
        });
      } catch {
        this.finishSearch(newSmarts ?? '');
      }
    }
  }

  async getFilterBitset(): Promise<BitArray> {
    console.log(`getFilterBitset currentSearches: ${this.currentSearches}`);
    const smarts = await this.sketcher.getSmarts();
    return await chemSubstructureSearchLibrary(this.column!, this.currentMolfile, smarts!, false, false,
      this.searchType, this.similarityCutOff, this.fp);
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
    };
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
