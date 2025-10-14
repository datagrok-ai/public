/* eslint-disable max-lines */
/**
 * Macromolecules substructure filter that uses Datagrok's collaborative filtering.
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import wu from 'wu';

import {Observable, Subject, Unsubscribable} from 'rxjs';

import {TAGS as bioTAGS, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {IRenderer} from '@datagrok-libraries/bio/src/types/renderer';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {linearSubstructureSearch} from '../substructure-search/substructure-search';
import {BioFilterBase, BioFilterProps, IBioFilter, IFilterProps}
  from '@datagrok-libraries/bio/src/substructure-filter/bio-substructure-filter-types';
import {HelmBioFilter} from './bio-substructure-filter-helm';
import {_package} from '../package';

const FILTER_SYNC_EVENT: string = 'bio-substructure-filter';

class FilterState {
  constructor(
    public readonly props: IFilterProps,
    public readonly filterId: number,
    public readonly dataFrameId: string,
    public readonly columnName: string,
    public readonly bitset: DG.BitSet | null,
  ) {}
}

export class SeparatorFilterProps extends BioFilterProps {
  constructor(
    substructure: string,
    public readonly separator?: string,
    logger?: DG.PackageLogger
  ) {
    super(substructure, false, logger);
    this.readOnly = true;
  }
}

export class BioSubstructureFilter extends DG.Filter implements IRenderer {
  bioFilter: IBioFilter | null = null;
  bitset: DG.BitSet | null = null;
  readonly loader: HTMLDivElement;
  notation: string | undefined = undefined;

  readonly filterSyncer: PromiseSyncer;

  get calculating(): boolean { return this.loader.style.display == 'initial'; }

  set calculating(value: boolean) { this.loader.style.display = value ? 'initial' : 'none'; }

  get filterSummary(): string {
    return this.bioFilter!.filterSummary;
  }

  get isFiltering(): boolean {
    return super.isFiltering && (this.bioFilter?.isFiltering ?? false);
  }

  get isReadyToApplyFilter(): boolean {
    return !this.calculating && this.bitset != null;
  }

  get debounceTime(): number {
    if (this.column == null)
      return 1000;
    const length = this.column.length;
    const minLength = 500;
    const maxLength = 10000;
    const msecMax = 1000;
    if (length < minLength) return 0;
    if (length > maxLength) return msecMax;
    const res = Math.floor(msecMax * ((length - minLength) / (maxLength - minLength)));
    return res;
  }

  constructor(
    private readonly seqHelper: ISeqHelper,
    private logger: ILogger
  ) {
    super();
    this.root = ui.divV([]);
    this.loader = ui.loader();
    this.calculating = false;
    this.filterSyncer = new PromiseSyncer(this.logger);
  }

  private static filterCounter: number = -1;
  private readonly filterId: number = ++BioSubstructureFilter.filterCounter;

  private filterToLog(): string { return `BioSubstructureFilter<${this.filterId}>`; }

  private viewSubs: Unsubscribable[] = [];

  attach(dataFrame: DG.DataFrame): void {
    const superAttach = super.attach.bind(this);
    const logPrefix = `${this.filterToLog()}.attach()`;
    this.filterSyncer.sync(logPrefix, async () => {
      superAttach(dataFrame);

      if (!this.column) {
        if (this.columnName)
          this.column = this.dataFrame!.getCol(this.columnName);
        else
          this.column = dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
      }
      const sh = this.seqHelper.getSeqHandler(this.column!);
      this.columnName ??= this.column?.name;
      this.notation ??= this.column?.meta.units!;

      this.bioFilter = this.notation === NOTATION.FASTA ?
        new FastaBioFilter() : this.notation === NOTATION.SEPARATOR ?
          new SeparatorBioFilter(this.column!.getTag(bioTAGS.separator)) : new HelmBioFilter(this.seqHelper);
      this.root.appendChild(this.bioFilter!.filterPanel);
      this.root.appendChild(this.loader);
      await this.bioFilter.attach(); // may await waitForElementInDom

      this.viewSubs.push(DG.debounce(this.bioFilter!.onChanged, this.debounceTime)
        .subscribe(this.bioFilterOnChangedDebounced.bind(this)));
      this.viewSubs.push(grok.events.onResetFilterRequest
        .subscribe(this.grokEventsOnResetFilterRequest.bind(this)));
      this.viewSubs.push(grok.events.onCustomEvent(FILTER_SYNC_EVENT)
        .subscribe(this.filterOnSync.bind(this)));
    });
  }

  detach() {
    const superDetach = super.detach.bind(this);
    const logPrefix = `${this.filterToLog()}.detach()`;
    this.filterSyncer.sync(logPrefix, async () => {
      for (const sub of this.viewSubs) sub.unsubscribe();
      this.viewSubs = [];
      superDetach(); // requests this.isFiltering
      if (this.bioFilter) this.bioFilter.detach();
      this.bioFilter = null;
    });
  }

  // -- Sync -

  private filterOnSync(state: FilterState): void {
    if (state.filterId === this.filterId) return;
    if (state.dataFrameId !== this.dataFrame!.id || state.columnName !== this.columnName) return;

    this.bioFilter!.props = state.props;
  }

  // -- Layout --

  applyFilter(): void {
    const logPrefix = `${this.filterToLog()}.applyFilter()`;
    this.logger.debug(`${logPrefix}, IN`);
    if (this.bitset && !this.isDetached)
      this.dataFrame?.filter.and(this.bitset);
  }

  /** Override to save filter state.
   * @return {any} - filter state
   */
  saveState(): any {
    const logPrefix = `${this.filterToLog()}.saveState()`;
    const state = super.saveState();
    this.logger.debug(`${logPrefix}, super.state = ${JSON.stringify(state)}`);
    state.props = this.bioFilter!.saveProps();
    return state;
  }

  /** Override to load filter state.
   * @param {any} state - filter state
   */
  applyState(state: any): void {
    const logPrefix = `${this.filterToLog()}.applyState()`;
    super.applyState(state); //column, columnName

    this.filterSyncer.sync(logPrefix, async () => {
      if (state.props && this.bioFilter)
        this.bioFilter.props = DG.toJs(state.props ?? {});
    });
  }

  private fireFilterSync(): void {
    const logPrefix = `${this.filterToLog()}.fireFilterSync()`;
    this.logger.debug(`${logPrefix}, ` +
      `bioFilter = ${!!this.bioFilter ? this.bioFilter.constructor.name : 'null'}` +
      (!!this.bioFilter ? `, props = ${JSON.stringify(this.bioFilter!.saveProps())}` : ''));

    grok.events.fireCustomEvent(FILTER_SYNC_EVENT, new FilterState(
      this.bioFilter!.props, this.filterId, this.dataFrame!.id, this.columnName!, this.bitset));
  }

  // -- Handle events

  /**
   * Performs the actual filtering
   * When the results are ready, triggers `rows.requestFilter`, which in turn triggers `applyFilter`
   * that would simply apply the bitset synchronously.
   */
  bioFilterOnChangedDebounced(): void {
    if (!this.dataFrame) return; // Debounced event can be handled postponed
    const logPrefix = `${this.filterToLog()}.bioFilterOnChangedDebounced()`;
    this.logger.debug(`${logPrefix}, start, ` +
      `isFiltering = ${this.isFiltering}, ` +
      `props = ${JSON.stringify(this.bioFilter!.saveProps())}`);

    if (!this.isFiltering) {
      this.bitset = null;
      this.dataFrame!.rows.requestFilter();
      return;
    }

    // some other filter is already filtering for the exact same thing
    if (wu(this.dataFrame!.rows.filters).has(`${this.columnName}: ${this.filterSummary}`))
      return;

    this.filterSyncer.sync(logPrefix, async () => {
      this.calculating = true;
      try {
        this.logger.debug(`${logPrefix}, before substructureSearch`);
        this.bitset = await this.bioFilter?.substructureSearch(this.column!)!;
        this.logger.debug(`${logPrefix}, after substructureSearch`);
        this.calculating = false;
        this.fireFilterSync();
        this.dataFrame?.rows.requestFilter();
      } finally {
        this.calculating = false;
        this.logger.debug(`${logPrefix}, end`);
      }
    });
  }

  grokEventsOnResetFilterRequest(): void {
    const logPrefix = `${this.filterToLog()}.grokEventsOnResetFilterRequest()`;
    this.logger.debug(`${logPrefix}`);
    this.bioFilter?.resetFilter();
  }

  // -- IRenderer --

  private _onRendered = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    const logPrefix = `${this.filterToLog()}.invalidate(${caller ? ` <- ${caller} ` : ''})`;
    this.filterSyncer.sync(logPrefix, async () => {
      // TODO: Request re-render
      this._onRendered.next();
    });
  }

  async awaitRendered(timeout: number = 10000): Promise<void> {
    const callLog = `awaitRendered( ${timeout} )`;
    const logPrefix = `${this.filterToLog()}.${callLog}`;
    await delay(10);
    await testEvent(this.onRendered, () => {
      this.logger.debug(`${logPrefix}, ` + '_onRendered event caught');
    }, () => {
      this.invalidate(callLog);
    }, timeout, `${logPrefix} timeout`);

    // Rethrow stored syncer error (for test purposes)
    const viewErrors = this.filterSyncer.resetErrors();
    if (viewErrors.length > 0) throw viewErrors[0];
  }
}

export class FastaBioFilter extends BioFilterBase<BioFilterProps> {
  readonly emptyProps = new BioFilterProps('', undefined, _package.logger);

  readonly substructureInput: DG.InputBase<string>;

  get type(): string { return 'FastaBioFilter'; }

  constructor() {
    super();

    this.substructureInput = ui.input.string('', {
      value: '', onValueChanged: (value) => {
        window.setTimeout(() => {
          this.props = new BioFilterProps(value, undefined, _package.logger);
          if (!this._propsChanging) this.onChanged.next();
        }, 0 /* next event cycle */);
      }, placeholder: 'Substructure'
    });
  }

  public applyProps() {
    if (this.substructureInput.value !== this.props.substructure)
      this.substructureInput.value = this.props.substructure;
  }

  get filterPanel() {
    return this.substructureInput.root;
  }

  get isFiltering(): boolean { return this.substructureInput.value !== ''; }

  async substructureSearch(column: DG.Column): Promise<DG.BitSet | null> {
    return linearSubstructureSearch(this.props.substructure, column);
  }

  async attach(): Promise<void> {}

  async detach(): Promise<void> {
    await super.detach();
  }
}

export class SeparatorBioFilter extends BioFilterBase<SeparatorFilterProps> {
  readonly emptyProps = new SeparatorFilterProps('', undefined, _package.logger);

  readonly substructureInput: DG.InputBase<string>;
  readonly separatorInput: DG.InputBase<string>;
  colSeparator = '';

  get type(): string { return 'SeparatorBioFilter'; }

  constructor(colSeparator: string) {
    super();

    this.substructureInput = ui.input.string('', {
      value: '', onValueChanged: (value) => {
        this.props = new SeparatorFilterProps(value, this.props.separator, _package.logger);
        setTimeout(() => {
          if (!this._propsChanging) this.onChanged.next();
        });
      }, placeholder: 'Substructure'
    });
    this.separatorInput = ui.input.string('', {
      value: this.colSeparator = colSeparator, onValueChanged: (value) => {
        const separator: string | undefined = !!value ? value : '';
        this.props = new SeparatorFilterProps(this.props.substructure, separator, _package.logger);
        setTimeout(() => {
          if (!this._propsChanging) this.onChanged.next();
        });
      }, placeholder: 'Separator'
    });
  }

  applyProps(): void {
    if (this.substructureInput.value !== this.props.substructure)
      this.substructureInput.value = this.props.substructure;

    const separatorValue = this.props.separator ?? this.colSeparator;
    if (this.separatorInput.value !== separatorValue)
      this.separatorInput.value = separatorValue;
  }

  get filterSummary(): string {
    const _sep: string = this.props.separator ? this.props.separator : this.colSeparator;
    return `${this.props.substructure}, {sep}`;
  };

  get isFiltering(): boolean { return this.props.substructure !== ''; };

  resetFilter(): void {
    this.props = new SeparatorFilterProps('', undefined, _package.logger);
  }

  get filterPanel() {
    return ui.divV([
      this.substructureInput.root,
      this.separatorInput.root,
    ]);
  }

  get substructure() {
    return this.separatorInput.value && this.separatorInput.value !== this.colSeparator ?
      this.substructureInput.value.replaceAll(this.separatorInput.value, this.colSeparator) :
      this.substructureInput.value;
  }

  set substructure(s: string) {
    this.substructureInput.value = s;
  }

  async substructureSearch(column: DG.Column): Promise<DG.BitSet | null> {
    return linearSubstructureSearch(this.substructure, column, this.colSeparator);
  }

  async attach(): Promise<void> {}

  async detach(): Promise<void> {
    await super.detach();
  }
}
