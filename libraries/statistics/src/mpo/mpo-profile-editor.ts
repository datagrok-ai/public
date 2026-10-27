import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subject, Subscription} from 'rxjs';
import {
  DEFAULT_AGGREGATION, WEIGHTED_AGGREGATIONS_LIST, MPO_WEIGHTING_CAPTIONS,
  DesirabilityMode, DesirabilityProfile, MpoWeighting, PropertyDesirability, WeightedAggregation,
  createDefaultCategorical, createDefaultNumerical, isMpoNumericColumn, isNumerical, migrateDesirability,
  rangeNumericalToColumn, applyDesirabilityPatch,
} from './mpo';
import {DesirabilityEditor} from './editors/desirability-editor-factory';
import {createMpoRow, MpoRow} from './editors/mpo-property-row';
import {discoverComputeFunctions, chemFunctionsDialog} from '../compute-functions';

import '../../css/styles.css';

import {MPO_SCORE_CHANGED_EVENT} from './utils';

const MAX_CATEGORICAL_CATEGORIES = 20;
const COLUMN_DROPDOWN_OPEN_SELECTOR = '.d4-column-selector-backdrop';

interface RowCtx {
  mpoRow: MpoRow;
  row: HTMLElement;
  editor: DesirabilityEditor;
  sub: Subscription;
  propertyCell?: HTMLElement;
  nameInput?: DG.InputBase<string>;
  nameInputSub?: Subscription;
  name: string;
  /** Correlation-aware weighting: the pin input, present while the property is being pinned. */
  pinInput?: DG.InputBase<number | null>;
  refreshWeight?: () => void;
}

const PIN_SUM_TOLERANCE = 1e-9;

function isMpoEligibleColumn(c: DG.Column): boolean {
  return c.isCategorical ? c.categories.length <= MAX_CATEGORICAL_CATEGORIES : isMpoNumericColumn(c);
}

export class MpoProfileEditor {
  readonly root = ui.div([]);
  readonly onChanged = new Subject<void>();
  readonly aggregationInput: DG.ChoiceInput<WeightedAggregation | null>;

  readonly onWeightingInputChanged = new Subject<void>();

  profile?: DesirabilityProfile;
  dataFrame?: DG.DataFrame;
  design = false;
  preview = false;
  /** Set by the host that owns the Weighting control (the profile view): shows priorities and pins. Other hosts
   *  show the weights of an Equal or Correlation-aware profile read-only. */
  weightingControl = false;
  /** "Why this weight" tooltip of a property name in Correlation-aware weighting. */
  explainWeight: ((name: string) => HTMLElement | null) | null = null;
  private baseWeights: Record<string, number> | null = null;

  private rowIds: Record<string, string> = {};
  private rowCtx = new Map<string, RowCtx>();
  private propertyOrder: string[] = [];
  columnMapping: Record<string, string | null> = {};

  constructor(dataFrame?: DG.DataFrame, design = false, preview = false) {
    this.dataFrame = dataFrame;
    this.design = design;
    this.preview = preview;
    this.aggregationInput = ui.input.choice('Aggregation', {
      items: WEIGHTED_AGGREGATIONS_LIST,
      value: DEFAULT_AGGREGATION,
      nullable: false,
      onValueChanged: (v) => {
        if (this.profile)
          this.profile.aggregation = v;
        this.emitChange();
      },
    });
    this.aggregationInput.setTooltip('Score aggregation method');
  }

  private newRowId(): string {
    return crypto.randomUUID();
  }

  private uniquePropertyName(): string {
    let n = 1;
    while (this.profile!.properties[`NewProperty${n}`])
      ++n;
    return `NewProperty${n}`;
  }

  setProfile(profile?: DesirabilityProfile): void {
    if (profile !== this.profile) {
      this.columnMapping = {};
      this.rowIds = {};
    }
    this.baseWeights = null;

    if (profile) {
      for (const key of Object.keys(profile.properties)) {
        const prop = migrateDesirability(profile.properties[key]);
        if (isNumerical(prop))
          prop.mode ??= DesirabilityMode.Freeform;
        profile.properties[key] = prop;
      }
    }

    this.resetRows();

    this.profile = profile;
    this.aggregationInput.value = profile?.aggregation ?? DEFAULT_AGGREGATION;
    this.propertyOrder = profile ? Object.keys(profile.properties) : [];

    for (const name of this.propertyOrder)
      this.rowIds[name] ??= this.newRowId();

    this.render();
    this.runAllComputeFunctions();
  }

  private disposeRowCtx(ctx: RowCtx): void {
    ctx.sub.unsubscribe();
    ctx.nameInputSub?.unsubscribe();
  }

  private resetRows(): void {
    for (const ctx of this.rowCtx.values())
      this.disposeRowCtx(ctx);
    this.rowCtx.clear();
  }

  private runAllComputeFunctions(): void {
    if (!this.dataFrame || !this.profile)
      return;
    for (const prop of Object.values(this.profile.properties)) {
      if (prop.computeFunction)
        this.runComputeFunction(prop);
    }
  }

  getProfile(): DesirabilityProfile | undefined {
    return this.profile;
  }

  get weighting(): MpoWeighting {
    return this.profile?.weighting ?? 'manual';
  }

  setWeighting(mode: MpoWeighting): void {
    if (!this.profile)
      return;
    if (mode === 'manual')
      delete this.profile.weighting;
    else
      this.profile.weighting = mode;
    this.baseWeights = null;
    this.resetRows();
    this.render();
    this.emitChange();
  }

  setPriority(name: string, priority: number): void {
    const prop = this.profile?.properties[name];
    if (!prop)
      return;
    if (priority === 1)
      delete prop.priority;
    else
      prop.priority = priority;
    this.rerenderWeightingInputs();
  }

  setPinnedWeight(name: string, weight: number | null): void {
    const prop = this.profile?.properties[name];
    if (!prop)
      return;
    if (weight == null)
      delete prop.pinnedWeight;
    else
      prop.pinnedWeight = weight;
    this.rerenderWeightingInputs();
  }

  private rerenderWeightingInputs(): void {
    this.resetRows();
    this.render();
    this.weightingInputChanged();
  }

  /** Writes computed weights (and the base weights the Priority hint needs) without a change per property;
   *  fires one change if any weight changed at display precision. Returns whether it did. */
  setWeights(weights: Record<string, number>, base: Record<string, number> | null): boolean {
    if (!this.profile)
      return false;
    this.baseWeights = base;
    let changed = false;
    for (const [name, w] of Object.entries(weights)) {
      const prop = this.profile.properties[name];
      if (!prop || prop.weight.toFixed(3) === w.toFixed(3))
        continue;
      prop.weight = w;
      changed = true;
    }
    for (const ctx of this.rowCtx.values())
      ctx.refreshWeight?.();
    if (changed)
      this.emitChange();
    return changed;
  }

  setDesignMode(on: boolean): void {
    if (this.design === on)
      return;
    this.design = on;

    this.resetRows();
    this.render();
  }

  private render(): void {
    ui.empty(this.root);

    if (!this.profile)
      return this.renderEmpty('No profile specified.');

    const rows = this.propertyOrder.map((name) => {
      const rowId = this.rowIds[name];
      return this.rowCtx.get(rowId)?.row ?? this.buildRow(name, rowId, this.profile!.properties[name]);
    });

    if (!rows.length) {
      if (!this.preview)
        this.root.append(this.buildHeader());
      if (this.design)
        return this.renderDesignEmpty();
      return this.renderEmpty('No properties defined.');
    }

    if (!this.preview)
      this.root.append(this.buildHeader());
    this.root.append(ui.divV(rows));

    this.revalidateNameInputs();
  }

  private renderEmpty(text: string): void {
    this.root.append(ui.divText(text));
  }

  private renderDesignEmpty(): void {
    const icon = ui.iconFA('chart-line');
    const heading = ui.h3('No properties yet');
    const msg = ui.p('Select a dataset to auto-populate properties from its numerical columns, or add them manually.');
    const addBtn = ui.link('+ Add Property', () => this.addProperty());
    const container = ui.divV([icon, heading, msg, addBtn], 'statistics-mpo-empty-state');
    this.root.append(container);
  }

  addProperty(name?: string): void {
    if (!this.profile)
      return;

    const newName = name ?? this.uniquePropertyName();
    const newRowId = this.newRowId();

    this.profile.properties[newName] = createDefaultNumerical();
    this.rowIds[newName] = newRowId;
    this.propertyOrder.push(newName);

    this.resetRows();
    this.render();
    this.emitChange();
  }

  updateProperty(name: string, patch: Partial<PropertyDesirability>): void {
    const ctx = this.rowCtx.get(this.rowIds[name]);
    if (ctx) {
      ctx.mpoRow.applyPatch(patch);
      return;
    }
    const prop = this.profile?.properties[name];
    if (!prop)
      return;
    applyDesirabilityPatch(prop, patch);
    this.emitChange();
  }

  setPropertyColumn(name: string, col: DG.Column | null): void {
    const rowId = this.rowIds[name];
    this.columnMapping[name] = col?.name ?? null;
    if (col && this.switchPropertyType(name, rowId, col))
      return;
    this.rowCtx.get(rowId)?.editor.setColumn?.(col);
    this.emitChange();
  }

  deleteProperty(name: string): void {
    this.deleteRow(this.rowIds[name]);
  }

  private get showsPriority(): boolean {
    return this.weightingControl && this.weighting === 'correlation';
  }

  private buildHeader(): HTMLElement {
    return ui.divH([
      ui.divText('Property', 'statistics-mpo-header-property'),
      this.showsPriority ? ui.divText('Priority', 'statistics-mpo-header-priority') : null,
      ui.divText('Weight', this.showsPriority ? 'statistics-mpo-header-weight-campo' : 'statistics-mpo-header-weight'),
      ui.divText('Desirability', 'statistics-mpo-header-desirability'),
    ].filter(Boolean) as HTMLElement[], 'statistics-mpo-header');
  }

  private buildRow(
    name: string,
    rowId: string,
    prop: PropertyDesirability,
  ): HTMLElement {
    const col = this.resolveColumn(name);
    if (col) {
      const corrected = this.correctPropertyType(prop, col);
      if (corrected) {
        prop = corrected;
        this.profile!.properties[name] = prop;
      }
    }

    const currentName = (): string => this.rowCtx.get(rowId)?.name ?? name;
    const mpoRow = createMpoRow(prop, {
      name: currentName,
      column: () => this.mappedColumn(currentName()),
      design: this.design,
      onChanged: () => this.emitChange(),
      onReplaced: (replacement) => {
        const current = currentName();
        const props = this.profile!.properties;
        props[current] = MpoProfileEditor.keepWeighting(props[current], replacement);
        this.rebuildRow(current, rowId);
      },
    });
    const {root: row, editor} = mpoRow;
    row.dataset.rowId = rowId;

    const ctx: RowCtx = {mpoRow, row, editor, sub: mpoRow.sub, name};
    ctx.propertyCell = this.buildPropertyCell(ctx);
    const columnCell = this.buildColumnSelector(ctx);
    mpoRow.propertyCell.append(...[ctx.propertyCell, columnCell].filter(Boolean) as HTMLElement[]);
    this.buildWeightCell(ctx, prop);
    if (this.design)
      row.append(this.buildRowControls(rowId, prop));
    this.rowCtx.set(rowId, ctx);

    return row;
  }

  private removeRow(rowId: string): void {
    const ctx = this.rowCtx.get(rowId);
    if (ctx)
      this.disposeRowCtx(ctx);
    this.rowCtx.delete(rowId);
  }

  private buildPropertyCell(ctx: RowCtx): HTMLElement {
    ctx.nameInputSub?.unsubscribe();
    ctx.nameInputSub = undefined;

    if (!this.design) {
      ctx.nameInput = undefined;
      return this.buildPropertyNameText(ctx);
    }
    ctx.nameInput = this.buildPropertyNameInput(ctx);
    return ctx.nameInput.root;
  }

  private buildPropertyNameText(ctx: RowCtx): HTMLElement {
    const el = ui.divText(ctx.name);
    ui.tooltip.bind(el, () => ctx.name);
    return el;
  }

  private buildPropertyNameInput(ctx: RowCtx): DG.InputBase<string> {
    const propNameInp = ui.input.string('', {value: ctx.name, nullable: false,
      onValueChanged: (v) => {
        const oldName = ctx.name;
        const trimmed = v.trim();
        if (trimmed === oldName || this.validatePropertyName(ctx, trimmed))
          return;
        if (this.updatePropertyName(oldName, trimmed))
          this.emitChange();
      },
    });
    propNameInp.addValidator((v) => this.validatePropertyName(ctx, v));
    ctx.nameInputSub = propNameInp.onInput.subscribe(() => this.revalidateNameInputs());

    return propNameInp;
  }

  private isPropertyNameTaken(ctx: RowCtx | null, name: string): boolean {
    for (const other of this.rowCtx.values()) {
      if (other === ctx)
        continue;
      const {nameInput, name: otherName} = other;
      if ((nameInput ? nameInput.value.trim() : otherName) === name)
        return true;
    }
    return false;
  }

  private validatePropertyName(ctx: RowCtx, v: string | null | undefined): string | null {
    const newName = v?.trim();
    if (!newName)
      return 'Property name cannot be empty';
    return this.isPropertyNameTaken(ctx, newName) ? 'A property with this name already exists' : null;
  }

  private revalidateNameInputs(): void {
    for (const ctx of this.rowCtx.values())
      ctx.nameInput?.validate();
  }

  private updatePropertyCell(rowId: string): void {
    const ctx = this.rowCtx.get(rowId);
    if (!ctx)
      return;

    const newCell = this.buildPropertyCell(ctx);
    ctx.propertyCell?.replaceWith(newCell);
    ctx.propertyCell = newCell;
  }

  private appendWeightColumnControls(
    cell: HTMLElement,
    weightInput: DG.InputBase<number | null>,
    prop: PropertyDesirability,
  ): void {
    if (!this.dataFrame)
      return;

    const update = (fn: (p: PropertyDesirability) => void): void => {
      fn(prop);
      this.emitChange();
    };

    const numCols = this.getNumericalColumnNames();
    let isColumn = !!prop.weightColumn && numCols.includes(prop.weightColumn);

    const colInput = ui.input.choice('', {items: numCols, nullable: true, value: prop.weightColumn ?? '',
      onValueChanged: (v) => update((p) => p.weightColumn = v || undefined),
    });

    const syncToggle = () => {
      weightInput.root.classList.toggle('statistics-mpo-hidden', isColumn);
      colInput.root.classList.toggle('statistics-mpo-hidden', !isColumn);
      toggle.classList.toggle('statistics-mpo-weight-toggle-active', isColumn);
    };

    const toggle = ui.iconFA('exchange-alt', () => {
      isColumn = !isColumn;
      syncToggle();
      if (!isColumn)
        update((p) => delete p.weightColumn);
      else if (colInput.value)
        update((p) => p.weightColumn = colInput.value || undefined);
    }, 'Toggle');
    toggle.classList.add('statistics-mpo-weight-toggle');
    ui.tooltip.bind(toggle, () => isColumn ? 'Switch to manual weight' : 'Use weight from column');
    syncToggle();

    cell.append(colInput.root, toggle);
  }

  /** Copies the Correlation-aware inputs of a property onto its replacement (a type switch builds a new one). */
  private static keepWeighting(from: PropertyDesirability, to: PropertyDesirability): PropertyDesirability {
    if (from.priority != null)
      to.priority = from.priority;
    if (from.pinnedWeight != null)
      to.pinnedWeight = from.pinnedWeight;
    return to;
  }

  private propOf(ctx: RowCtx): PropertyDesirability {
    return this.profile!.properties[ctx.name];
  }

  private buildWeightCell(ctx: RowCtx, prop: PropertyDesirability): void {
    const {weightCell, weightInput} = ctx.mpoRow;
    const mode = this.weighting;
    if (mode === 'manual')
      return this.appendWeightColumnControls(weightCell, weightInput, prop);

    weightInput.enabled = false;
    ctx.refreshWeight = () => {
      weightInput.notify = false;
      weightInput.value = this.propOf(ctx).weight;
      weightInput.notify = true;
    };
    if (!this.weightingControl) {
      weightInput.setTooltip(`Set by ${MPO_WEIGHTING_CAPTIONS[mode]} weighting. Change it in MPO Profiles.`);
      return;
    }
    if (mode === 'equal')
      return;

    weightInput.root.classList.add('statistics-mpo-hidden');
    const priority = this.buildPriorityCell(ctx);
    weightCell.before(priority.root);
    ui.tooltip.bind(ctx.mpoRow.propertyCell, () => this.explainWeight?.(ctx.name) ?? null);

    if (prop.pinnedWeight != null)
      ctx.pinInput = this.createPinInput(prop.pinnedWeight);
    const host = ui.divH([], 'statistics-mpo-campo-weight');
    const valueEl = ui.divText('', 'statistics-mpo-campo-value');
    // the name input keeps the hover for its own validation tooltip, so the weight explains itself too
    ui.tooltip.bind(valueEl, () => this.explainWeight?.(ctx.name) ?? null);
    const render = () => {
      ui.empty(host);
      const icon = ui.iconFA('thumbtack', () => this.togglePin(ctx, render));
      icon.classList.add('statistics-mpo-pin');
      if (ctx.pinInput)
        icon.classList.replace('fal', 'fas');
      icon.classList.toggle('statistics-mpo-pin-active', !!ctx.pinInput);
      ui.tooltip.bind(icon, () => ctx.pinInput ? 'Unpin weight' : 'Pin weight');
      host.append(ctx.pinInput?.root ?? valueEl, icon);
      ctx.refreshWeight!();
    };
    ctx.refreshWeight = () => {
      valueEl.textContent = this.propOf(ctx).weight.toFixed(3);
      priority.refresh();
    };
    render();
    weightCell.append(host);
  }

  private buildPriorityCell(ctx: RowCtx): {root: HTMLElement, refresh: () => void} {
    const input = ui.input.float('', {value: this.propOf(ctx).priority ?? 1, min: 0, step: 0.1, format: '#0.00',
      nullable: true,
      onValueChanged: (v) => {
        if (!input.validate())
          return;
        const p = this.propOf(ctx);
        if (v === 1)
          delete p.priority;
        else
          p.priority = v!;
        refresh();
        this.weightingInputChanged();
      },
    });
    input.addValidator((s) => s == null || s.trim() === '' ? 'Enter a priority' :
      !(Number(s) >= 0) ? 'Must be 0 or greater' : null);
    input.captionLabel.classList.add('statistics-mpo-hidden');
    ui.tooltip.bind(input.input, 'Multiplies the computed weight');

    const hint = ui.divText('', 'statistics-mpo-priority-hint');
    const refresh = () => {
      const p = this.propOf(ctx);
      const changed = (p.priority ?? 1) !== 1;
      hint.textContent = !changed ? '' : p.pinnedWeight != null ? 'Applies when the weight is unpinned.' :
        this.baseWeights?.[ctx.name] === 0 ? 'No effect: computed weight is 0. Pin the weight instead.' : '';
      hint.classList.toggle('statistics-mpo-hidden', !hint.textContent);
    };
    refresh();
    return {root: ui.divV([input.root, hint], 'statistics-mpo-priority-cell'), refresh};
  }

  private createPinInput(value?: number): DG.InputBase<number | null> {
    const input = ui.input.float('', {value, min: 0, max: 1, step: 0.01, format: '#0.000', nullable: true,
      onValueChanged: () => {
        this.syncPins();
        this.weightingInputChanged();
      },
    });
    input.addValidator((s) => {
      if (s == null || s.trim() === '')
        return 'Enter a weight';
      const v = Number(s);
      if (!Number.isFinite(v) || v < 0 || v > 1)
        return 'Must be between 0 and 1';
      const sum = this.pinnedSum();
      return sum > 1 + PIN_SUM_TOLERANCE ? `Pinned weights sum to ${sum.toFixed(2)} (> 1)` : null;
    });
    input.root.classList.add('statistics-mpo-weight-input');
    input.captionLabel.classList.add('statistics-mpo-hidden');
    return input;
  }

  private pinnedSum(): number {
    let s = 0;
    for (const ctx of this.rowCtx.values()) {
      const v = ctx.pinInput?.value;
      if (v != null && Number.isFinite(v))
        s += v;
    }
    return s;
  }

  /** The sum condition spans rows, so every pin is revalidated. A property is pinned while its own value is a valid
   *  weight; a sum over 1 keeps the pins and is reported by the computation, which then applies nothing. */
  private syncPins(): void {
    for (const ctx of this.rowCtx.values()) {
      const p = this.propOf(ctx);
      ctx.pinInput?.validate();
      const v = ctx.pinInput?.value;
      if (v != null && v >= 0 && v <= 1)
        p.pinnedWeight = v;
      else
        delete p.pinnedWeight;
    }
  }

  private togglePin(ctx: RowCtx, render: () => void): void {
    const pinning = !ctx.pinInput;
    // pinning starts from the computed weight, so it changes nothing until the user edits it
    ctx.pinInput = pinning ? this.createPinInput(Number(this.propOf(ctx).weight.toFixed(3))) : undefined;
    render();
    if (pinning)
      setTimeout(() => (ctx.pinInput?.input as HTMLInputElement | undefined)?.select());
    this.syncPins();
    this.weightingInputChanged();
  }

  private weightingInputChanged(): void {
    this.onWeightingInputChanged.next();
    for (const ctx of this.rowCtx.values())
      ctx.refreshWeight?.();
    this.emitChange();
  }

  private buildColumnSelector(ctx: RowCtx): HTMLElement | null {
    if (!this.dataFrame)
      return null;

    const matchedCol = this.mappedColumn(ctx.name);

    const commit = (v: DG.Column | null): void => this.setPropertyColumn(ctx.name, v);

    // Local workaround to avoid blocking statistics releases on a js-api change.
    // Proper fix: expose `changeOnHover` on IColumnInputInitOptions.
    const isPopupOpen = (): boolean => !!document.querySelector(COLUMN_DROPDOWN_OPEN_SELECTOR);

    let pending: {value: DG.Column | null} | null = null;
    let observer: MutationObserver | null = null;

    const stopObserving = (): void => {
      observer?.disconnect();
      observer = null;
    };

    const flush = (): void => {
      if (!pending)
        return;
      const v = pending.value;
      pending = null;
      commit(v);
    };

    const watchForClose = (): void => {
      if (observer)
        return;
      observer = new MutationObserver(() => {
        if (isPopupOpen())
          return;
        stopObserving();
        flush();
      });
      observer.observe(document.body, {childList: true, subtree: true});
    };

    const input = ui.input.column('', {
      table: this.dataFrame,
      nullable: true,
      value: matchedCol ?? undefined,
      filter: (c: DG.Column) => isMpoEligibleColumn(c),
      onValueChanged: (v: DG.Column | null) => {
        if (isPopupOpen()) {
          pending = {value: v};
          watchForClose();
          return;
        }
        pending = null;
        stopObserving();
        commit(v);
      },
    });

    // Block d4's mouse-wheel column cycling so page scroll over the input doesn't
    // mutate the value. No preventDefault — the browser still scrolls.
    input.root.addEventListener('wheel', (e: Event) => e.stopImmediatePropagation(), {capture: true, passive: true});

    return input.root;
  }

  private getEligibleColumnNames(): string[] {
    if (!this.dataFrame)
      return [];
    return Array.from(this.dataFrame.columns)
      .filter((c) => isMpoEligibleColumn(c))
      .map((c) => c.name);
  }

  private getNumericalColumnNames(): string[] {
    if (!this.dataFrame)
      return [];
    return Array.from(this.dataFrame.columns).filter((c) => isMpoNumericColumn(c)).map((c) => c.name);
  }

  private mappedColumn(name: string): DG.Column | null {
    const colName = this.columnMapping[name];
    return colName ? this.dataFrame?.col(colName) ?? null : null;
  }

  private resolveColumn(name: string): DG.Column | null {
    if (!this.dataFrame)
      return null;
    const eligible = this.getEligibleColumnNames();
    const colName = this.columnMapping[name] ?? eligible.find((c) => c.toLowerCase() === name.toLowerCase()) ?? null;
    if (colName)
      this.columnMapping[name] = colName;
    return colName ? this.dataFrame.col(colName) ?? null : null;
  }

  private correctPropertyType(prop: PropertyDesirability, col: DG.Column): PropertyDesirability | null {
    if (isNumerical(prop) && col.isCategorical) {
      const categories = col.categories.map((c: string) => ({name: c, desirability: 1}));
      return MpoProfileEditor.keepWeighting(prop, createDefaultCategorical(prop.weight, categories));
    }

    if (!isNumerical(prop) && col.isNumerical)
      return MpoProfileEditor.keepWeighting(prop, createDefaultNumerical(prop.weight, col.min, col.max));

    if (isNumerical(prop) && col.isNumerical && !prop.rangeUserSet && (prop.min !== col.min || prop.max !== col.max)) {
      const updated = {...prop};
      rangeNumericalToColumn(updated, col);
      return updated;
    }

    return null;
  }

  private switchPropertyType(
    name: string,
    rowId: string,
    col: DG.Column,
  ): boolean {
    if (!this.profile)
      return false;

    const corrected = this.correctPropertyType(this.profile.properties[name], col);
    if (!corrected)
      return false;

    this.profile.properties[name] = corrected;
    this.rebuildRow(name, rowId);
    return true;
  }

  private rebuildRow(name: string, rowId: string): void {
    const oldRow = this.rowCtx.get(rowId)?.row;
    this.removeRow(rowId);
    const newRow = this.buildRow(name, rowId, this.profile!.properties[name]);

    if (oldRow?.parentNode)
      oldRow.replaceWith(newRow);

    this.revalidateNameInputs();
    this.emitChange();
  }

  private buildRowControls(rowId: string, prop: PropertyDesirability): HTMLElement {
    const compute = this.buildComputeIcon(prop);
    const add = ui.icons.add(() => this.insertRowAfterRow(rowId), 'Add');
    const del = ui.icons.delete(() => this.deleteRow(rowId), 'Delete');
    return ui.divH([compute, add, del], 'statistics-mpo-control-buttons');
  }

  renameProperty(oldName: string, newName: string): boolean {
    if (!this.updatePropertyName(oldName, newName))
      return false;

    this.updatePropertyCell(this.rowIds[newName]);
    this.emitChange();
    return true;
  }

  private updatePropertyName(oldName: string, newName: string): boolean {
    if (!this.profile || oldName === newName)
      return false;

    const rowId = this.rowIds[oldName];
    if (this.isPropertyNameTaken(this.rowCtx.get(rowId) ?? null, newName))
      return false;

    this.profile.properties[newName] = this.profile.properties[oldName];
    delete this.profile.properties[oldName];

    this.columnMapping[newName] = this.columnMapping[oldName] ?? null;
    delete this.columnMapping[oldName];

    const idx = this.propertyOrder.indexOf(oldName);
    if (idx >= 0)
      this.propertyOrder[idx] = newName;

    delete this.rowIds[oldName];
    this.rowIds[newName] = rowId;

    const ctx = this.rowCtx.get(rowId);
    if (ctx)
      ctx.name = newName;
    return true;
  }

  private deleteRow(rowId: string): void {
    const name = this.rowCtx.get(rowId)?.name;
    if (!name || !this.profile)
      return;

    delete this.profile.properties[name];
    delete this.columnMapping[name];
    delete this.rowIds[name];
    const rowEl = this.rowCtx.get(rowId)?.row;
    this.removeRow(rowId);

    const idx = this.propertyOrder.indexOf(name);
    if (idx >= 0)
      this.propertyOrder.splice(idx, 1);

    if (this.propertyOrder.length === 0)
      this.render();
    else
      rowEl?.remove();

    this.revalidateNameInputs();
    this.emitChange();
  }

  private insertRowAfterRow(rowId: string): void {
    if (!this.profile)
      return;

    const afterName = this.rowCtx.get(rowId)?.name;
    if (!afterName)
      return;

    const newName = this.uniquePropertyName();
    const newRowId = this.newRowId();

    this.profile.properties[newName] = createDefaultNumerical();

    this.rowIds[newName] = newRowId;

    const idx = this.propertyOrder.indexOf(afterName);
    this.propertyOrder.splice(idx + 1, 0, newName);

    const newRow = this.buildRow(newName, newRowId, this.profile.properties[newName]);
    this.rowCtx.get(rowId)?.row.after(newRow);
    this.emitChange();
  }

  private buildComputeIcon(prop: PropertyDesirability): HTMLElement {
    return ui.iconFA('calculator', async () => {
      const computeFunctions = discoverComputeFunctions('HitTriageFunction');
      const template = {
        compute: {
          descriptors: {enabled: false, args: []},
          functions: prop.computeFunction ? [prop.computeFunction] : [],
        },
      };
      await chemFunctionsDialog(
        computeFunctions,
        (result) => {
          const funcNames = Object.keys(result.externals);
          if (funcNames.length > 0) {
            const key = funcNames[0];
            const [pkg, name] = key.split(':');
            prop.computeFunction = {package: pkg, name, args: result.externals[key]};
          } else
            prop.computeFunction = undefined;
          this.runComputeFunction(prop);
          this.emitChange();
        },
        () => {},
        template,
        true,
        undefined,
        undefined,
        true,
      );
    }, 'Configure compute function');
  }

  private async runComputeFunction(prop: PropertyDesirability): Promise<void> {
    if (!this.dataFrame || !prop.computeFunction)
      return;
    const molCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    if (!molCol)
      return;
    const fs = DG.Func.find({package: prop.computeFunction.package, name: prop.computeFunction.name});
    if (!fs.length)
      return;
    const f = fs[0];
    const tablePropName = f.inputs[0].name;
    const colPropName = f.inputs[1].name;
    await f.apply({...prop.computeFunction.args, [tablePropName]: this.dataFrame, [colPropName]: molCol.name});
  }

  private emitChange(): void {
    this.onChanged.next();
    grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
  }
}
