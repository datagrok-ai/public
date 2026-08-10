import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subject, Subscription} from 'rxjs';
import {
  DEFAULT_AGGREGATION, WEIGHTED_AGGREGATIONS_LIST,
  DesirabilityMode, DesirabilityProfile, PropertyDesirability, WeightedAggregation,
  createDefaultCategorical, createDefaultNumerical, isMpoNumericColumn, isNumerical, migrateDesirability,
} from './mpo';
import {DesirabilityEditor, DesirabilityEditorFactory} from './editors/desirability-editor-factory';
import {DesirabilityModeDialog} from './dialogs/desirability-mode-dialog';
import {discoverComputeFunctions, chemFunctionsDialog} from '../compute-functions';

import '../../css/styles.css';

import {MPO_SCORE_CHANGED_EVENT} from './utils';

const MAX_CATEGORICAL_CATEGORIES = 20;
const COLUMN_DROPDOWN_OPEN_SELECTOR = '.d4-column-selector-backdrop';

interface RowCtx {
  row: HTMLElement;
  editor: DesirabilityEditor;
  editorHost: HTMLElement;
  prop: PropertyDesirability;
  sub: Subscription;
  propertyCell?: HTMLElement;
  nameInput?: DG.InputBase<string>;
  nameInputSub?: Subscription;
  name: string;
  modeGear?: HTMLElement;
  controls?: HTMLElement;
}

function isMpoEligibleColumn(c: DG.Column): boolean {
  return c.isCategorical ? c.categories.length <= MAX_CATEGORICAL_CATEGORIES : isMpoNumericColumn(c);
}

export class MpoProfileEditor {
  readonly root = ui.div([]);
  readonly onChanged = new Subject<void>();
  readonly aggregationInput: DG.ChoiceInput<WeightedAggregation | null>;

  profile?: DesirabilityProfile;
  dataFrame?: DG.DataFrame;
  design = false;
  preview = false;

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

  setDesignMode(on: boolean): void {
    if (this.design === on)
      return;
    this.design = on;

    if (!this.profile || this.propertyOrder.length === 0)
      return this.render();

    for (const name of this.propertyOrder) {
      const rowId = this.rowIds[name];
      this.rowCtx.get(rowId)?.editor.setDesignMode?.(on);
      this.updatePropertyCell(rowId);
      this.updateDesignControls(rowId);
    }
    this.revalidateNameInputs();
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

  addProperty(): void {
    if (!this.profile)
      return;

    const newName = this.uniquePropertyName();
    const newRowId = this.newRowId();

    this.profile.properties[newName] = createDefaultNumerical();
    this.rowIds[newName] = newRowId;
    this.propertyOrder.push(newName);

    this.resetRows();
    this.render();
    this.emitChange();
  }

  private buildHeader(): HTMLElement {
    return ui.divH([
      ui.divText('Property', 'statistics-mpo-header-property'),
      ui.divText('Weight', 'statistics-mpo-header-weight'),
      ui.divText('Desirability', 'statistics-mpo-header-desirability'),
    ], 'statistics-mpo-header');
  }

  private buildRow(
    name: string,
    rowId: string,
    prop: PropertyDesirability,
  ): HTMLElement {
    const row = ui.divH([], 'statistics-mpo-row');
    row.dataset.rowId = rowId;

    const col = this.resolveColumn(name);
    if (col) {
      const corrected = this.correctPropertyType(prop, col);
      if (corrected) {
        prop = corrected;
        this.profile!.properties[name] = prop;
      }
    }

    const editor = DesirabilityEditorFactory.create(prop, 300, 80, this.design);
    editor.root.classList.add('statistics-mpo-editor-fill');
    const sub = editor.onChanged.subscribe(() => this.emitChange());
    const editorHost = ui.divH([editor.root]);

    const ctx: RowCtx = {row, editor, editorHost, prop, sub, name};
    ctx.propertyCell = this.buildPropertyCell(ctx);
    const weightCell = this.buildWeightCell(prop);
    const columnCell = this.buildColumnSelector(rowId, ctx, editor);
    this.rowCtx.set(rowId, ctx);

    row.append(
      ui.divV([ctx.propertyCell, columnCell].filter(Boolean), 'statistics-mpo-property-cell'),
      weightCell,
      editorHost,
    );

    this.updateDesignControls(rowId);

    return row;
  }

  private updateDesignControls(rowId: string): void {
    const ctx = this.rowCtx.get(rowId);
    if (!ctx)
      return;

    if (this.design) {
      if (!ctx.modeGear) {
        ctx.modeGear = this.buildModeGear(rowId, ctx.prop, ctx.editor);
        ctx.editorHost.append(ctx.modeGear);
      }
      if (!ctx.controls) {
        ctx.controls = this.buildRowControls(rowId, ctx.prop);
        ctx.row.append(ctx.controls);
      }
      return;
    }

    ctx.modeGear?.remove();
    ctx.modeGear = undefined;
    ctx.controls?.remove();
    ctx.controls = undefined;
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
        if (this.renameProperty(oldName, trimmed))
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

  private buildWeightCell(prop: PropertyDesirability): HTMLElement {
    const children: HTMLElement[] = [];
    const update = (fn: (p: PropertyDesirability) => void): void => {
      fn(prop);
      this.emitChange();
    };

    const weightInput = ui.input.float('', {value: prop.weight, min: 0, max: 1, format: '#0.000',
      onValueChanged: (v) => update((p) => p.weight = Math.max(0, Math.min(1, v ?? 0))),
    });
    weightInput.root.classList.add('statistics-mpo-weight-input');
    children.push(weightInput.root);

    if (this.dataFrame) {
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

      children.push(colInput.root, toggle);
    }

    return ui.divH(children, 'statistics-mpo-weight-cell');
  }

  private buildColumnSelector(
    rowId: string,
    ctx: RowCtx,
    editor: DesirabilityEditor,
  ): HTMLElement | null {
    if (!this.dataFrame)
      return null;

    const matchedName = this.columnMapping[ctx.name] ?? null;
    const matchedCol = matchedName ? this.dataFrame.col(matchedName) : null;

    if (matchedCol)
      editor.setColumn?.(matchedCol);

    const commit = (v: DG.Column | null): void => {
      this.columnMapping[ctx.name] = v?.name ?? null;
      if (v && this.switchPropertyType(ctx.name, rowId, v))
        return;
      editor.setColumn?.(v);
      this.emitChange();
    };

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
      return createDefaultCategorical(prop.weight, categories);
    }

    if (!isNumerical(prop) && col.isNumerical)
      return createDefaultNumerical(prop.weight, col.min, col.max);

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

  private buildModeGear(
    rowId: string,
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
  ): HTMLElement {
    return ui.icons.settings(() => {
      const name = this.rowCtx.get(rowId)?.name;
      if (!name)
        return;
      const colName = this.columnMapping[name];
      const col = colName ? this.dataFrame?.col(colName) ?? null : null;
      this.openModeDialog(name, rowId, prop, editor, col);
    }, 'Settings');
  }

  private buildRowControls(rowId: string, prop: PropertyDesirability): HTMLElement {
    const compute = this.buildComputeIcon(prop);
    const add = ui.icons.add(() => this.insertRowAfterRow(rowId), 'Add');
    const del = ui.icons.delete(() => this.deleteRow(rowId), 'Delete');
    return ui.divH([compute, add, del], 'statistics-mpo-control-buttons');
  }

  private openModeDialog(
    name: string,
    rowId: string,
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
    mappedCol: DG.Column | null = null,
  ): void {
    new DesirabilityModeDialog(
      name,
      prop,
      (patch) => {
        const p = this.profile?.properties[name];
        if (p)
          Object.assign(p, patch);
        editor.redrawAll();
      },
      (newProp) => {
        if (!this.profile)
          return;
        this.profile.properties[name] = newProp;
        this.rebuildRow(name, rowId);
      },
      mappedCol,
    ).show();
  }

  private renameProperty(oldName: string, newName: string): boolean {
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

    this.emitChange();
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
