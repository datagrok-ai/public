import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subject, Subscription} from 'rxjs';
import {
  DEFAULT_AGGREGATION, WEIGHTED_AGGREGATIONS_LIST,
  DesirabilityProfile, PropertyDesirability, WeightedAggregation,
  createDefaultCategorical, createDefaultNumerical, isNumerical, migrateDesirability,
} from './mpo';
import {DesirabilityEditor, DesirabilityEditorFactory} from './editors/desirability-editor-factory';
import {DesirabilityModeDialog} from './dialogs/desirability-mode-dialog';
import {discoverComputeFunctions} from '../compute-functions/discovery';
import {chemFunctionsDialog} from '../compute-functions/dialog';

import '../../css/styles.css';

import {MPO_SCORE_CHANGED_EVENT} from './utils';

const MAX_CATEGORICAL_CATEGORIES = 20;

export class MpoProfileEditor {
  readonly root = ui.div([]);
  readonly onChanged = new Subject<void>();
  readonly aggregationInput: DG.ChoiceInput<WeightedAggregation | null>;

  profile?: DesirabilityProfile;
  dataFrame?: DG.DataFrame;
  design = false;
  preview = false;

  private rows: Record<string, HTMLElement> = {};
  private rowIds: Record<string, string> = {};
  private rowSubs = new Map<string, Subscription>();
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

  setProfile(profile?: DesirabilityProfile): void {
    if (profile) {
      for (const key of Object.keys(profile.properties)) {
        const prop = migrateDesirability(profile.properties[key]);
        if (isNumerical(prop))
          prop.mode ??= 'freeform';
        profile.properties[key] = prop;
      }
    }

    for (const sub of this.rowSubs.values())
      sub.unsubscribe();
    this.rowSubs.clear();

    this.profile = profile;
    this.aggregationInput.value = profile?.aggregation ?? DEFAULT_AGGREGATION;
    this.columnMapping = {};
    this.rows = {};
    this.rowIds = {};
    this.propertyOrder = profile ? Object.keys(profile.properties) : [];

    for (const name of this.propertyOrder)
      this.rowIds[name] = this.newRowId();

    this.render();
    this.runAllComputeFunctions();
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
    this.rows = {};
    this.render();
  }

  setPreviewMode(on: boolean): void {
    if (this.preview === on)
      return;
    this.preview = on;
    this.rows = {};
    this.render();
  }

  private render(): void {
    ui.empty(this.root);

    if (!this.profile)
      return this.renderEmpty('No profile specified.');

    const rows = this.propertyOrder.map((name) => {
      const rowId = this.rowIds[name];
      if (!this.rows[rowId])
        this.rows[rowId] = this.buildRow(name, rowId, this.profile!.properties[name]);
      return this.rows[rowId];
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

    const newName = `NewProperty${Object.keys(this.profile.properties).length + 1}`;
    const newRowId = this.newRowId();

    this.profile.properties[newName] = createDefaultNumerical();
    this.rowIds[newName] = newRowId;
    this.propertyOrder.push(newName);

    this.rows = {};
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
    this.rowSubs.get(rowId)?.unsubscribe();
    this.rowSubs.set(rowId, editor.onChanged.subscribe(() => this.emitChange()));

    const propertyCell = this.buildPropertyCell(rowId, name);
    const weightCell = this.buildWeightCell(rowId, prop);
    const columnCell = this.buildColumnSelector(rowId, name, editor);

    const modeGear = this.design ?
      this.buildModeGear(rowId, prop, editor) :
      null;

    const controls = this.design ? this.buildRowControls(rowId, prop) : null;

    row.append(
      ui.divV([propertyCell, columnCell].filter(Boolean), 'statistics-mpo-property-cell'),
      weightCell,
      ui.divH([editor.root, modeGear].filter(Boolean)),
    );

    if (controls)
      row.append(controls);

    return row;
  }

  private buildPropertyCell(rowId: string, name: string): HTMLElement | null {
    if (this.dataFrame) {
      const el = ui.divText(name);
      ui.tooltip.bind(el, () => name);
      return el;
    }

    let currentName = name;

    const propNameInp = ui.input.string('', {value: name, onValueChanged: (v) => {
      if (!v || v === currentName)
        return;
      this.renameProperty(currentName, v);
      currentName = v;
    }});

    ui.tooltip.bind(propNameInp.input, () => currentName);
    return propNameInp.root;
  }

  private buildWeightCell(
    rowId: string,
    prop: PropertyDesirability,
  ): HTMLElement {
    const name = this.getPropertyNameByRowId(rowId);
    const children: HTMLElement[] = [];

    const weightInput = ui.input.float('', {value: prop.weight, min: 0, max: 1, format: '#0.000',
      onValueChanged: (v) => {
        if (name)
          this.mutateProperty(name, (p) => p.weight = Math.max(0, Math.min(1, v ?? 0)));
      },
    });
    weightInput.root.classList.add('statistics-mpo-weight-input');
    children.push(weightInput.root);

    if (this.dataFrame) {
      const numCols = this.getNumericalColumnNames();
      let isColumn = !!prop.weightColumn && numCols.includes(prop.weightColumn);

      const colInput = ui.input.choice('', {items: numCols, nullable: true, value: prop.weightColumn ?? '',
        onValueChanged: (v) => {
          if (name)
            this.mutateProperty(name, (p) => p.weightColumn = v || undefined);
        },
      });

      const syncToggle = () => {
        weightInput.root.classList.toggle('statistics-mpo-hidden', isColumn);
        colInput.root.classList.toggle('statistics-mpo-hidden', !isColumn);
        toggle.classList.toggle('statistics-mpo-weight-toggle-active', isColumn);
      };

      const toggle = ui.iconFA('exchange-alt', () => {
        isColumn = !isColumn;
        syncToggle();
        if (!isColumn && name)
          this.mutateProperty(name, (p) => delete p.weightColumn);
        else if (isColumn && name && colInput.value)
          this.mutateProperty(name, (p) => p.weightColumn = colInput.value || undefined);
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
    name: string,
    editor: DesirabilityEditor,
  ): HTMLElement | null {
    if (!this.dataFrame)
      return null;

    const items = this.getEligibleColumnNames();
    const matched = this.columnMapping[name] ?? null;

    if (matched) {
      const col = this.dataFrame.col(matched);
      editor.setColumn?.(col);
    }

    const input = ui.input.choice('', {items, nullable: true, value: matched ?? '', onValueChanged: (v) => {
      this.columnMapping[name] = v ?? null;
      const col = v ? this.dataFrame!.col(v) : null;
      if (col && this.switchPropertyType(name, rowId, col))
        return;
      editor.setColumn?.(col);
      this.emitChange();
    }});
    return input.root;
  }

  private getEligibleColumnNames(): string[] {
    if (!this.dataFrame)
      return [];
    return Array.from(this.dataFrame.columns)
      .filter((c) => !c.isCategorical || c.categories.length <= MAX_CATEGORICAL_CATEGORIES)
      .map((c) => c.name);
  }

  private getNumericalColumnNames(): string[] {
    if (!this.dataFrame)
      return [];
    return Array.from(this.dataFrame.columns.numerical).map((c) => c.name);
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
    const oldRow = this.rows[rowId];
    const newRow = this.buildRow(name, rowId, this.profile!.properties[name]);
    this.rows[rowId] = newRow;

    if (oldRow?.parentNode)
      oldRow.replaceWith(newRow);

    this.emitChange();
  }

  private buildModeGear(
    rowId: string,
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
  ): HTMLElement {
    const gear = ui.icons.settings(() => {
      const name = this.getPropertyNameByRowId(rowId);
      if (!name)
        return;
      const colName = this.columnMapping[name];
      const col = colName ? this.dataFrame?.col(colName) ?? null : null;
      this.openModeDialog(name, rowId, prop, editor, col);
    }, 'Settings');

    return gear;
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

  private getPropertyNameByRowId(rowId: string): string | undefined {
    return Object.entries(this.rowIds)
      .find(([, id]) => id === rowId)?.[0];
  }

  private renameProperty(oldName: string, newName: string): void {
    if (!this.profile || this.profile.properties[newName])
      return;

    const rowId = this.rowIds[oldName];
    this.profile.properties[newName] = this.profile.properties[oldName];
    delete this.profile.properties[oldName];

    this.columnMapping[newName] = this.columnMapping[oldName] ?? null;
    delete this.columnMapping[oldName];

    const idx = this.propertyOrder.indexOf(oldName);
    if (idx >= 0)
      this.propertyOrder[idx] = newName;

    delete this.rowIds[oldName];
    this.rowIds[newName] = rowId;

    this.emitChange();
  }

  private deleteRow(rowId: string): void {
    const name = this.getPropertyNameByRowId(rowId);
    if (!name || !this.profile)
      return;

    delete this.profile.properties[name];
    delete this.columnMapping[name];
    delete this.rowIds[name];
    this.rowSubs.get(rowId)?.unsubscribe();
    this.rowSubs.delete(rowId);

    const idx = this.propertyOrder.indexOf(name);
    if (idx >= 0)
      this.propertyOrder.splice(idx, 1);

    const rowEl = this.rows[rowId];
    delete this.rows[rowId];

    if (this.propertyOrder.length === 0)
      this.render();
    else
      rowEl?.remove();

    this.emitChange();
  }

  private insertRowAfterRow(rowId: string): void {
    if (!this.profile)
      return;

    const afterName = this.getPropertyNameByRowId(rowId);
    if (!afterName)
      return;

    const newName = `NewProperty${Object.keys(this.profile.properties).length + 1}`;
    const newRowId = this.newRowId();

    this.profile.properties[newName] = createDefaultNumerical();

    this.rowIds[newName] = newRowId;

    const idx = this.propertyOrder.indexOf(afterName);
    this.propertyOrder.splice(idx + 1, 0, newName);

    const newRow = this.buildRow(newName, newRowId, this.profile.properties[newName]);
    this.rows[newRowId] = newRow;

    this.rows[rowId]?.after(newRow);
    this.emitChange();
  }

  private mutateProperty(
    name: string,
    updater: (p: PropertyDesirability) => void,
  ): void {
    const prop = this.profile?.properties[name];
    if (!prop)
      return;

    updater(prop);
    this.emitChange();
  }

  private buildComputeIcon(prop: PropertyDesirability): HTMLElement {
    const icon = ui.iconFA('calculator', async () => {
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
      );
    }, 'Configure compute function');
    return icon;
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
