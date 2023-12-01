/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {FILTER_TYPES, chemSubstructureSearchLibrary} from '../chem-searches';
import {initRdKitService} from '../utils/chem-common-rdkit';
import {Subject, Subscription} from 'rxjs';
import {debounceTime, filter} from 'rxjs/operators';
import wu from 'wu';
import {TaskBarProgressIndicator, chem} from 'datagrok-api/dg';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {_package, getRdKitModule} from '../package';
import {AVAILABLE_FPS, CHEM_APPLY_FILTER_SYNC, FILTER_SCAFFOLD_TAG, MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT,
  SubstructureSearchType, getSearchProgressEventName, getSearchQueryAndType, getTerminateEventName} from '../constants';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { IColoredScaffold } from '../rendering/rdkit-cell-renderer';
import { Fingerprint } from '../utils/chem-common';
import $ from 'cash-dom';

const FILTER_SYNC_EVENT = 'chem-substructure-filter';
const SKETCHER_TYPE_CHANGED = 'chem-sketcher-type-changed';
let id = 0;

const searchTypeHints  = {
  [SubstructureSearchType.CONTAINS]: 'search structures which contain sketched pattern as a substructure',
  [SubstructureSearchType.INCLUDED_IN]: 'search structures for which sketched pattern is a superstructure',
  [SubstructureSearchType.NOT_CONTAINS]: 'search structures which DO NOT contain sketched pattern as a substructure',
  [SubstructureSearchType.NOT_INCLUDED_IN]: 'search structures for which sketched pattern is NOT a superstructure',
  [SubstructureSearchType.EXACT_MATCH]: 'search structures which exactly match sketched pattern',
  [SubstructureSearchType.IS_SIMILAR]: 'search structures similar to sketched pattern',
}

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

interface IApplyFilterSync{
  totalSubstrFiltersOnCol: number[];
  activeFilterId?: number;
  applyFilterCallsCounter: number;
}

export class SubstructureFilter extends DG.Filter {
  // @ts-ignore
  sketcher: DG.chem.Sketcher = new DG.chem.Sketcher();
  bitset: DG.BitSet | null = null;
  loader: HTMLDivElement = ui.loader();
  onSketcherChangedSubs?: Subscription[] = [];
  active: boolean = true;
  syncEvent = false;
  searchTypeSync = false;
  similarityCutOffSync = false;
  fpSync = false;
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
    SubstructureSearchType.EXACT_MATCH, SubstructureSearchType.IS_SIMILAR,
    SubstructureSearchType.NOT_CONTAINS, SubstructureSearchType.NOT_INCLUDED_IN];
  fpsTypes = AVAILABLE_FPS;
  searchTypeInput: DG.InputBase;
  similarityCutOffInput: DG.InputBase;
  fpInput: DG.InputBase;
  similarityOptionsDiv = ui.divH([], 'chem-filter-similarity-options');
  sketcherDiv = ui.div('', {style: {position: 'relative', width: '100%'}})
  emptySketcherDiv = ui.divH([], 'empty-filter');
  optionsIcon: HTMLElement;
  searchTypeChanged = new Subject();
  searchOptionsDiv = ui.div('', 'filter-search-options');
  showOptions = false;
  searchNotCompleted = false;
  
  get calculating(): boolean {return this.loader.style.display == 'initial';}
  set calculating(value: boolean) {this.loader.style.display = value ? 'initial' : 'none';}

  get filterSummary(): string {
    return this.getFilterSummary(this.currentMolfile);
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
      this.onSearchTypeChanged();
    });
    ui.tooltip.bind(this.searchTypeInput.input, () => {
      console.log(`in tooltip`);
      return searchTypeHints[this.searchTypeInput.value as SubstructureSearchType]
    });

    this.fpInput = ui.choiceInput('FP', this.fp, this.fpsTypes, () => {
      this.fp = this.fpInput.value;
      !this.fpSync ? this.searchTypeChanged.next() : this.fpSync = false;
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
      !this.similarityCutOffSync ? this.searchTypeChanged.next() : this.similarityCutOffSync = false;
    });

    ui.tooltip.bind((this.similarityCutOffInput.root.getElementsByClassName('ui-input-slider')[0] as HTMLElement)!, 'Similarity cutoff');


    this.optionsIcon = ui.icons.settings(() => {
      this.onShowOptionsChanged();
    });
    $(this.optionsIcon).addClass('chem-search-options-icon');

    this.sketcherDiv.append(this.sketcher.root);
    this.emptySketcherDiv.append(this.searchTypeInput.root);
    this.emptySketcherDiv.append(this.sketcherDiv);
    this.similarityOptionsDiv.append(this.fpInput.root);
    this.similarityOptionsDiv.append(this.similarityCutOffInput.root);
    this.root.appendChild(this.searchOptionsDiv);
    this.root.appendChild(this.emptySketcherDiv);
    this.root.appendChild(this.loader);
    this.root.classList.add('chem-filter');
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
    };
    //this fix is required for Marvin JS sync between filter panel and hamburger menu 
    ui.tools.waitForElementInDom(this.sketcher.root).then(() => {
      if (this.sketcher._mode === DG.chem.SKETCHER_MODE.INPLACE) {
        this.root.append(this.sketcher.root);
        this.searchOptionsDiv.append(this.searchTypeInput.root);
      } else {
        this.updateFilterUiOnSketcherChanged(this.currentMolfile);
      }
    });
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
        this.searchTypeInput.value = SubstructureSearchType.CONTAINS;
        this.sketcher.setMolFile(DG.WHITE_MOLBLOCK);
      }
    }));
    this.subs.push(grok.events.onCustomEvent(FILTER_SYNC_EVENT).subscribe((state: ISubstructureFilterState) => {
      if (state.colName === this.columnName && this.tableName == state.tableName && this.filterId !== state.filterId) {
        /* setting syncEvent to true only if base sketcher is initialized.
        If base sketcher is initialized, it will fire onChange event */      
        if (this.currentSearches.size > 0)
          grok.events.fireCustomEvent(this.terminateEventName, this.currentSearches.values().next().value);
        if (this.sketcher.sketcher?.isInitialized || this.sketcher._mode == DG.chem.SKETCHER_MODE.INPLACE) 
          this.syncEvent = true;
        this.currentMolfile = state.molblock!;
        this.bitset = state.bitset!;
        this.searchTypeSync = true;
        this.searchTypeInput.value = state.searchType;
        this.similarityCutOffInput.value = state.simCutOff;
        this.similarityCutOffSync = true;
        this.fpInput.value = state.fp;
        this.fpSync = true;
        this.sketcher.setMolFile(state.molblock!);
        this.updateFilterUiOnSketcherChanged(state.molblock!);
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

    this.subs.push(this.dataFrame!.onEvent('d4-filter-control-active-changed').subscribe(() => {
      //in case filter is swithed off via checkbox on filter panel, we finish current search
      if (!this.isFiltering && this.bitset) {
        if (!this.batchResultObservable?.closed) {
          this.searchNotCompleted = true;
          this.sketcher.getSmarts().then((smarts) => {
            this.terminatePreviousSearch();
            this.finishSearch(getSearchQueryAndType(smarts, this.searchType, this.fp, this.similarityCutOff));
          });
        }
      }
    }));

    this.currentSearches.add('');
    chemSubstructureSearchLibrary(this.column!, '', '', FILTER_TYPES.substructure, false, false)
      .then((_) => { }); // Precalculating fingerprints

    let onChangedEvent: any = this.sketcher.onChanged;
    //onChangedEvent = onChangedEvent.pipe(debounceTime(this._debounceTime));
    this.onSketcherChangedSubs?.push(onChangedEvent.subscribe(async (_: any) => {
      _package.logger.debug(`in filter onChangedEvent, sync event: ${this.syncEvent}`);
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
      this.finishSearch(getSearchQueryAndType(smarts, this.searchType, this.fp, this.similarityCutOff));
      this.onSketcherChangedSubs?.forEach((it) => it.unsubscribe());
    });
    grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: DG.BitSet.create(this.dataFrame!.rowCount).setAll(true),
      molblock: DG.WHITE_MOLBLOCK, colName: this.columnName, filterId: this.filterId, 
      tableName: this.tableName, searchType: this.searchType, simCutOff: this.similarityCutOff, fp: this.fp});
    super.detach();
    if (this.column?.temp[FILTER_SCAFFOLD_TAG])
      this.column.temp[FILTER_SCAFFOLD_TAG] = null;
    this.batchResultObservable?.unsubscribe();
  }

  applyFilter(): void {
    if (!this.updateApplyFilterSyncTag(false))
      return;
    if (this.dataFrame && this.bitset && !this.isDetached) {
      this.dataFrame.filter.and(this.bitset);
        this.dataFrame.rows.addFilterState(this.saveState());
        this.column!.temp[FILTER_SCAFFOLD_TAG] = JSON.stringify([{
          molecule: this.currentMolfile,
          isSuperstructure: this.searchType === SubstructureSearchType.INCLUDED_IN
        }]);
        this.active = true;
        if (this.searchNotCompleted)
          this._onSketchChanged();
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

    const applyFilterSyncTag: IApplyFilterSync = this.column?.getTag(CHEM_APPLY_FILTER_SYNC) ?
    JSON.parse(this.column!.getTag(CHEM_APPLY_FILTER_SYNC)!) : {totalSubstrFiltersOnCol: [], applyFilterCallsCounter: 0};
    if (!applyFilterSyncTag.totalSubstrFiltersOnCol.includes(this.filterId))
      applyFilterSyncTag.totalSubstrFiltersOnCol.push(this.filterId);
    this.column!.setTag(CHEM_APPLY_FILTER_SYNC, JSON.stringify(applyFilterSyncTag));

    if (!this.initListeners) {
      this.initListeners = true;
      this.terminateEventName = getTerminateEventName(this.tableName, this.columnName!);
      this.progressEventName = getSearchProgressEventName(this.tableName, this.columnName!);

      this.subs.push(grok.events.onCustomEvent(this.terminateEventName).subscribe((queryMol: string) => {
        _package.logger.debug(`in ${this.terminateEventName} handler, querymol: ${queryMol}`);
        this.finishSearch(queryMol);
      }));
    }
    this.active = state.active ?? true;
    if (this.column?.temp[FILTER_SCAFFOLD_TAG])
      state.molBlock ??= (JSON.parse(this.column?.temp[FILTER_SCAFFOLD_TAG]) as IColoredScaffold[])[0].molecule;
    if (state.molBlock) {
      this.currentMolfile = state.molBlock;
      this.sketcher.setMolFile(state.molBlock);
      this.updateFilterUiOnSketcherChanged(this.currentMolfile);
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
    _package.logger.debug(`newMolfile ${newMolFile}`);
    const newSmarts = await this.getSmartsToFilter();
    _package.logger.debug(`newSmarts ${newSmarts}`);
    if (this.currentMolfile !== newMolFile)
      this.updateFilterUiOnSketcherChanged(newMolFile);
    grok.events.fireCustomEvent(SKETCHER_TYPE_CHANGED, {colName: this.columnName,
      filterId: this.filterId, tableName: this.tableName});
    if (!this.isFiltering) {
      _package.logger.debug(`not filtering ${newSmarts}`);
      this.currentMolfile = newMolFile;
      this.bitset = !this.active ?
        DG.BitSet.fromBytes((await this.getFilterBitset())!.buffer.buffer, this.column!.length) : null; //TODO
      if (this.column?.temp[FILTER_SCAFFOLD_TAG])
        delete this.column.temp[FILTER_SCAFFOLD_TAG];
      this.terminatePreviousSearch();
      this.finishSearch(getSearchQueryAndType(newSmarts, this.searchType, this.fp, this.similarityCutOff));
      grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
        molblock: this.currentMolfile, colName: this.columnName, filterId: this.filterId, 
        tableName: this.tableName, searchType: this.searchType, simCutOff: this.similarityCutOff, fp: this.fp});
      this.dataFrame?.rows.requestFilter();
    } else if (wu(this.dataFrame!.rows.filters)
      .has(`${this.columnName}: ${this.getFilterSummary(newMolFile)}`) && !this.searchNotCompleted) {
      // some other filter is already filtering for the exact same thing
      // value to pass into has() is created similarly to filterSummary property
      _package.logger.debug(`already filter by the same structure ${this.getFilterSummary(newMolFile)}`);
      return;
    } else {
      this.searchNotCompleted = false;
      this.terminatePreviousSearch();
      this.currentMolfile = newMolFile;
      this.currentSearches.add(getSearchQueryAndType(newSmarts, this.searchType, this.fp, this.similarityCutOff));
      this.calculating = true;
      this.progressBar ??= DG.TaskBarProgressIndicator.create(`Starting substructure search...`);
      _package.logger.debug(`starting filter by ${this.currentMolfile}`);
      try {
        grok.events.fireCustomEvent(FILTER_SYNC_EVENT, {bitset: this.bitset,
          molblock: this.currentMolfile, colName: this.columnName, filterId: this.filterId, 
          tableName: this.tableName, searchType: this.searchType, simCutOff: this.similarityCutOff, fp: this.fp});
        const bitArray = await this.getFilterBitset();
        this.bitset = DG.BitSet.fromBytes(bitArray.buffer.buffer, this.column!.length);
        this.batchResultObservable?.unsubscribe();
        this.batchResultObservable = grok.events.onCustomEvent(this.progressEventName).subscribe((progress: number) => {
          _package.logger.debug(`progress: ${progress}, molfile: ${this.currentMolfile}}`);
          this.bitset = DG.BitSet.fromBytes(bitArray.buffer.buffer, this.column!.length);
          this.updateApplyFilterSyncTag();
          this.dataFrame?.rows.requestFilter();
            this.progressBar?.update(progress, `${progress?.toFixed(2)}% of search completed`);
        });
      } catch {
        this.finishSearch(getSearchQueryAndType(newSmarts, this.searchType, this.fp, this.similarityCutOff));
      }
    }
  }

  async getSmartsToFilter() {
    return this.sketcher.sketcher?.isInitialized ? await this.sketcher.getSmarts() :
      _convertMolNotation(this.currentMolfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts, getRdKitModule());
  }

  updateApplyFilterSyncTag(add = true): boolean {
    const applyFilterSyncTag: IApplyFilterSync = JSON.parse(this.column!.getTag(CHEM_APPLY_FILTER_SYNC)!);
    if (add) {
      applyFilterSyncTag.applyFilterCallsCounter = applyFilterSyncTag.totalSubstrFiltersOnCol.length;
      applyFilterSyncTag.activeFilterId = this.filterId;
    } else {
      if (applyFilterSyncTag.applyFilterCallsCounter === 0)
        return true;
      applyFilterSyncTag.applyFilterCallsCounter -=1;
      if (this.filterId === applyFilterSyncTag.activeFilterId)
        return true;
    }
    this.column!.setTag(CHEM_APPLY_FILTER_SYNC, JSON.stringify(applyFilterSyncTag));
    return false;
  }

  async getFilterBitset(): Promise<BitArray> {
    const smarts = await this.getSmartsToFilter();
    return await chemSubstructureSearchLibrary(this.column!, this.currentMolfile, smarts!, FILTER_TYPES.substructure, false, false,
      this.searchType, this.similarityCutOff, this.fp);
  }

  updateExternalSketcher() {
    if (this.sketcher._mode === DG.chem.SKETCHER_MODE.EXTERNAL)
      this.sketcher.updateExtSketcherContent(); //updating image in minimized sketcher panel
  }

  updateFilterUiOnSketcherChanged(newMolFile: string){
    if (this.sketcher._mode !== DG.chem.SKETCHER_MODE.INPLACE) {
      if (!!newMolFile && !chem.Sketcher.isEmptyMolfile(newMolFile)){
        this.removeChildIfExists(this.root, this.emptySketcherDiv, 'empty-filter');
        this.root.appendChild(this.sketcherDiv);
        if (this.searchType === SubstructureSearchType.CONTAINS) {
          this.sketcher.root.appendChild(this.optionsIcon);
        } else {
          this.root.prepend(this.searchOptionsDiv);
          this.searchOptionsDiv.append(this.searchTypeInput.root);
          if (this.searchType === SubstructureSearchType.IS_SIMILAR)
            this.searchOptionsDiv.append(this.similarityOptionsDiv);
          else
            this.removeChildIfExists(this.searchOptionsDiv, this.similarityOptionsDiv, 'chem-filter-similarity-options');
        }
      }
      else {
        this.emptySketcherDiv.append(this.searchTypeInput.root);
        this.emptySketcherDiv.append(this.sketcherDiv);
        this.root.append(this.emptySketcherDiv);
        this.removeChildIfExists(this.sketcher.root, this.optionsIcon, 'chem-search-options-icon');
        if (this.searchType !== SubstructureSearchType.IS_SIMILAR)
          this.removeChildIfExists(this.searchOptionsDiv, this.similarityOptionsDiv, 'chem-filter-similarity-options');
      }
    }
    this.updateExternalSketcher();
  }

  removeChildIfExists(parent: HTMLElement, child: HTMLElement, className: string) {
    if (parent.getElementsByClassName(className).length)
      parent.removeChild(child);
  }

  onSearchTypeChanged() {
    this.searchType = this.searchTypeInput.value;
    if (this.searchType !== SubstructureSearchType.CONTAINS)
      this.removeChildIfExists(this.sketcher.root, this.optionsIcon, 'chem-search-options-icon');
    else {
      if (!chem.Sketcher.isEmptyMolfile(this.sketcher.getMolFile()) && this.sketcher._mode !== DG.chem.SKETCHER_MODE.INPLACE)
        this.sketcher.root.appendChild(this.optionsIcon);
    }
    if (this.searchType === SubstructureSearchType.IS_SIMILAR)
      this.searchOptionsDiv.append(this.similarityOptionsDiv);
    else
      this.removeChildIfExists(this.searchOptionsDiv, this.similarityOptionsDiv, 'chem-filter-similarity-options');
      !this.searchTypeSync ? this.searchTypeChanged.next() : this.searchTypeSync = false;
  }

  onShowOptionsChanged() {
    this.showOptions = !this.showOptions;
    if (this.showOptions) {
      this.root.prepend(this.searchOptionsDiv);
      this.searchOptionsDiv.append(this.searchTypeInput.root);
    } else
      this.removeChildIfExists(this.root, this.searchOptionsDiv, 'filter-search-options');
  }

  terminatePreviousSearch() {
    if (this.currentSearches.size > 0)
      grok.events.fireCustomEvent(this.terminateEventName, this.currentSearches.values().next().value);
  }

  finishSearch(queryMolAndType: string) {
    _package.logger.debug(`in finishSearch ${this.currentSearches}`);
    const finish = () => {
      _package.logger.debug(`in finish function ${queryMolAndType}`);
      if (this.currentSearches.size === 0) {
        this.calculating = false;
        this.progressBar?.close();
        this.progressBar = null;
        this.batchResultObservable?.unsubscribe();
        console.log(`Unsubscribed from batchResultObservable  Filter ${this.filterId}`);
      }
    };
    if (this.currentSearches.has(queryMolAndType)) {
      _package.logger.debug(`in finishSearch first if ${queryMolAndType}`);
      this.currentSearches.delete(queryMolAndType);
      finish();
    }
    if (queryMolAndType == null && this.currentSearches.size === 1) {
      _package.logger.debug(`in finishSearch second if ${queryMolAndType}`);
      const v = this.currentSearches.values().next().value;
      if (v != null)
        this.currentSearches.delete(v);
      finish();
    }
  }

  getFilterSummary(molfile: string): string {
    const smarts = _convertMolNotation(molfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts,
      getRdKitModule());
    return this.searchType === SubstructureSearchType.IS_SIMILAR ?
      `${smarts}_${this.searchType}_${this.similarityCutOff}_${this.fp}` : `${smarts}_${this.searchType}`;
  }
}
