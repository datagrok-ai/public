import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {CHEM_SIMILARITY_METRICS} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import '../../css/chem.css';
import {Subject, Subscription} from 'rxjs';
import {AVAILABLE_FPS} from '../constants';
import {pickTextColorBasedOnBgColor} from '../utils/ui-utils';

export const SIMILARITY = 'similarity';
export const DIVERSITY = 'diversity';
// Shared default cap (kept small — the diversity viewer's selection is O(n^2)). The similarity
// viewer raises its own `limit` property max to MAX_LIMIT_SIMILARITY in its constructor.
export const MAX_LIMIT = 50;
export enum RowSourceTypes {
  All = 'All',
  Filtered = 'Filtered',
  Selected = 'Selected',
  FilteredSelected = 'FilteredSelected',
}
export class ChemSearchBaseViewer extends DG.JsViewer {
  isEditedFromSketcher: boolean = false;
  gridSelect: boolean = false;
  name: string = '';
  distanceMetric: string;
  limit: number;
  fingerprint: string;
  metricsProperties = ['distanceMetric', 'fingerprint'];
  fingerprintChoices = AVAILABLE_FPS;
  rowSourceChoices = [RowSourceTypes.All, RowSourceTypes.Filtered, RowSourceTypes.Selected, RowSourceTypes.FilteredSelected];
  sizesMap: {[key: string]: {[key: string]: number}} = {
    'small': {height: 60, width: 120},
    'normal': {height: 100, width: 200},
    'large': {height: 150, width: 300}};
  size: string;
  moleculeColumn?: DG.Column|null;
  moleculeColumnName: string;
  initialized: boolean = false;
  metricsDiv: HTMLElement | null;
  moleculeProperties: string[];
  renderCompleted = new Subject<void>();
  isComputing = false;
  rowSource: string;
  error = '';

  constructor(name: string, col?: DG.Column) {
    super();
    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0],
      {choices: this.fingerprintChoices, category: 'Similarity search'});
    this.limit = this.int('limit', 12, {min: 1, max: MAX_LIMIT, category: 'Similarity search'});
    this.distanceMetric = this.string('distanceMetric', CHEM_SIMILARITY_METRICS[0],
      {choices: CHEM_SIMILARITY_METRICS, category: 'Similarity search'});
    this.size = this.string('size', Object.keys(this.sizesMap)[0], {choices: Object.keys(this.sizesMap)});
    // Per-viewer default (e.g. the similarity viewer raises it to 'Filtered' in its own constructor).
    this.rowSource = this.string('rowSource', this.rowSourceChoices[0], {choices: this.rowSourceChoices});
    this.moleculeColumnName = this.addProperty('moleculeColumnName', DG.TYPE.COLUMN, '', {semType: DG.SEMTYPE.MOLECULE});
    this.name = name;
    this.moleculeProperties = this.columnList('moleculeProperties', [],
      {description: 'Adds selected fields from the grid to similarity search viewer'});
    if (col) {
      this.moleculeColumn = col;
      this.moleculeColumnName = col.name!;
    }
    const header = this.name === DIVERSITY ? `Most diverse structures` : `Most similar structures`;
    this.metricsDiv = ui.divH([ui.divText(header)], 'chem-similarity-header');
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    if (this._debRenderTimeout)
      clearTimeout(this._debRenderTimeout);
    this.subs.forEach((sub) => sub.unsubscribe());
    // Clear the array too: onTableAttached pushes fresh subscriptions, so without this a re-attach would
    // accumulate dead subs (and a second live onResetFilterRequest would toggle the filter twice).
    this.subs = [];
  }

  async onTableAttached(): Promise<void> {
    // Drop any prior subscriptions before re-subscribing (mirrors detach), so a re-attach WITHOUT a
    // preceding detach can't accumulate duplicate subs — e.g. a second onResetFilterRequest that would
    // make "Reset all filters" toggle the similarity filter twice.
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
    // Clamp a restored/over-set limit to this viewer's own property max (project restore writes the
    // field directly, bypassing onPropertyChanged) — 50 for diversity, MAX_LIMIT_SIMILARITY for similarity.
    const limitMax = this.getProperty('limit')?.max;
    if (limitMax != null && this.limit > limitMax)
      this.limit = limitMax;
    this.init();
    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.onRowsRemoved, 50).subscribe(async (_: any) => await this.render()));
      const compute = this.name !== DIVERSITY;
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50)
        .subscribe(async (_: any) => {
          if (this.isEditedFromSketcher && !this.gridSelect)
            this.isEditedFromSketcher = false;
          await this.render(compute);
        }));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
        .subscribe(async (_: any) =>
          await this.render(this.rowSource === RowSourceTypes.Selected || this.rowSource === RowSourceTypes.FilteredSelected)));
      this.subs.push(DG.debounce((this.dataFrame.onFilterChanged), 50)
        .subscribe(async (_: any) => await this.onExternalFilterChanged()));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50)
        .subscribe(async (_: any) => await this.render(false)));
      this.subs.push(DG.debounce((this.dataFrame.onMetadataChanged), 50)
        .subscribe(async (_: any) => await this.render(false)));
      this.moleculeColumn = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      this.moleculeColumnName = this.moleculeColumn?.name ?? '';
    }
    await this.render(true);
  }


  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (!this.initialized)
      return;
    if (this.metricsProperties.includes(property?.name))
      this.updateMetricsLink(this, {});
    if (property?.name === 'moleculeColumnName') {
      const col = this.dataFrame.col(property.get(this));
      this.moleculeColumn = col;
    }
    // Clamp `limit` to the property's own (per-viewer) max — 50 for diversity, MAX_LIMIT_SIMILARITY for
    // similarity — so a value written via setOptions (which bypasses the UI slider's max) stays bounded.
    if (property?.name === 'limit' && this.limit > property.max)
      this.limit = property.max;
    if (property.name === 'moleculeProperties') {
      this.debouncedRender(false);
      return;
    }
    this.debouncedRender();
  }

  updateMetricsLink(object: any, options: {[key: string]: string}): void {
    const metricsButton = ui.link(` ${this.distanceMetric}, ${this.fingerprint}`, () => {
      if (!grok.shell.windows.showProperties)
        grok.shell.windows.showProperties = true;
      grok.shell.o = object;
    }, 'Distance metric and fingerprint', '');
    metricsButton.classList.add('chem-similarity-metrics-link');
    Object.keys(options).forEach((it: any) => metricsButton.style[it] = options[it]);
    // Locate the previous link by class (not by child index) so that any action
    // icons injected into metricsDiv are not clobbered when this is refreshed.
    const existingLink = this.metricsDiv!.querySelector('.chem-similarity-metrics-link');
    if (existingLink)
      this.metricsDiv!.removeChild(existingLink);
    this.metricsDiv!.appendChild(metricsButton);
  }

  private _debRenderTimeout: ReturnType<typeof setTimeout> | null = null;
  private _debComputeFlag = false;
  protected async debouncedRender(computeData = true): Promise<void> {
    if (this._debRenderTimeout)
      clearTimeout(this._debRenderTimeout);
    this._debComputeFlag = this._debComputeFlag || computeData;
    this._debRenderTimeout = setTimeout(async () => {
      const flag = this._debComputeFlag;
      this._debComputeFlag = false;
      await this.render(flag);
    }, 200);
  }

  private _renderEpoch = 0;
  async render(computeData = true): Promise<void> {
    const epoch = ++this._renderEpoch;
    try {
      if (!this.moleculeColumn) {
        ui.empty(this.root);
        this.root.append(ui.div(ui.divText(`No molecule columns found in dataset`, {style: {color: 'red'}}),
          'chem-similarity-diversity-search-error'));
        return;
      }
      await this.renderInternal(computeData);
    } finally {
      // Only the LATEST render emits completion, so renderCompleted reflects the winning render — not a
      // stale one that bailed mid-flight (a fast row-navigation can otherwise fire it before the winner).
      if (this.isComputing && epoch === this._renderEpoch) {
        this.isComputing = false;
        this.renderCompleted.next();
      }
    }
  }

  /** Called (debounced) when the dataframe's filter changes. Subclasses may override to coalesce or
   * defer the re-render; the base simply re-renders, recomputing only in filter-driven row sources. */
  protected async onExternalFilterChanged(): Promise<void> {
    await this.render(this.rowSource === RowSourceTypes.Filtered || this.rowSource === RowSourceTypes.FilteredSelected);
  }

  async renderInternal(compute = true) {

  }

  beforeRender() {
    return this.initialized && this.dataFrame;
  }

  createMoleculePropertiesDiv(idx: number, refMolecule: boolean, similarity?: number): HTMLDivElement {
    const propsDict: {[key: string]: any} = {};
    if (!grok.shell.tv)
      return ui.div();
    if (similarity) {
      if (refMolecule)
        propsDict['Reference'] = {val: ''};
      else
        propsDict[SIMILARITY] = {val: similarity};
    }
    for (const col of this.moleculeProperties) {
      const propCol = this.moleculeColumn!.dataFrame.col(col);
      if (propCol) {
        propsDict[col] = {val: propCol.getString(idx)};
        const colorCoding = propCol.meta.colors.getType();
        if (colorCoding && colorCoding !== DG.COLOR_CODING_TYPE.OFF) {
          propsDict[col].color = this.dataFrame.col(col)?.meta.colors.getColor(idx);
          propsDict[col].isTextColorCoded = grok.shell.tv.grid?.col(col)?.isTextColorCoded;
        }
      }
    }
    //const item = ui.divH([], 'similarity-prop-item');
    const div = ui.divV([], {style: {marginTop: '5px'}});
    for (const key of Object.keys(propsDict)) {
      const labelName = key === SIMILARITY ? '' : key;
      const label = ui.divText(`${labelName}`, 'chem-similarity-prop-label');
      const value = ui.divText(`${propsDict[key].val}`, 'chem-similarity-prop-value');
      ui.tooltip.bind(value, key);
      if (propsDict[key].color) {
        const color = DG.Color.toHtml(propsDict[key].color);
        if (!propsDict[key].isTextColorCoded) {
          value.style.backgroundColor = color;
          value.style.color = DG.Color.toHtml(DG.Color.getContrastColor(propsDict[key].color));
        } else
          value.style.color = color;
      }
      const item = ui.divH([
        label,
        value,
      ], 'chem-similarity-prop-item');
      div.append(item);
    }
    return div;
  }

  getPropsColumnsNames(): string[] {
    let fingerprintTag = '';
    for (const t of this.moleculeColumn!.tags.keys()) {
      if (t.endsWith('.Column')) {
        fingerprintTag = t;
        break;
      }
    }
    return this.moleculeColumn!.dataFrame.columns.names()
      .filter((name) => name !== this.moleculeColumn!.getTag(fingerprintTag) && name !== this.moleculeColumn!.name);
  }

  closeWithError(error: string, progressBar?: DG.TaskBarProgressIndicator | null) {
    this.error = error;
    this.clearResults();
    this.root.append(ui.divText(this.error));
    this.root.classList.add(`chem-malformed-molecule-error`);
    progressBar?.close();
  }

  clearResults() {
    if (this.root.hasChildNodes())
      this.root.removeChild(this.root.childNodes[0]);
  }

  getRowSourceIndexes(): DG.BitSet {
    const bitset = DG.BitSet.create(this.dataFrame.rowCount);
    switch (this.rowSource) {
    case RowSourceTypes.All:
      bitset.setAll(true);
      break;
    case RowSourceTypes.Filtered:
      bitset.copyFrom(this.dataFrame.filter);
      break;
    case RowSourceTypes.Selected:
      bitset.copyFrom(this.dataFrame.selection);
      break;
    case RowSourceTypes.FilteredSelected:
      const filterCopy = DG.BitSet.create(this.dataFrame.rowCount).copyFrom(this.dataFrame.filter);
      bitset.copyFrom(filterCopy.and(this.dataFrame.selection));
      break;
    default:
      break;
    }
    return bitset;
  }
}
