import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subject, Subscription} from 'rxjs';
import {
  DesirabilityProfile, PropertyDesirability,
  createDefaultCategorical, createDefaultNumerical, isNumerical, migrateDesirability,
} from './mpo';
import {DesirabilityEditor, DesirabilityEditorFactory} from './editors/desirability-editor-factory';
import {DesirabilityModeDialog} from './dialogs/desirability-mode-dialog';

import '../../css/styles.css';

import {MPO_SCORE_CHANGED_EVENT} from './utils';

const MAX_CATEGORICAL_CATEGORIES = 20;

export class MpoProfileEditor {
  readonly root = ui.div([]);
  readonly onChanged = new Subject<void>();

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
    this.columnMapping = {};
    this.rows = {};
    this.rowIds = {};
    this.propertyOrder = profile ? Object.keys(profile.properties) : [];

    for (const name of this.propertyOrder)
      this.rowIds[name] = this.newRowId();

    this.render();
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

    if (!rows.length)
      return this.renderEmpty('No properties defined.');

    if (!this.preview)
      this.root.append(this.buildHeader());
    this.root.append(ui.divV(rows));
  }

  private renderEmpty(text: string): void {
    this.root.append(ui.divText(text));
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
    this.rowSubs.get(rowId)?.unsubscribe();
    this.rowSubs.set(rowId, editor.onChanged.subscribe(() => this.emitChange()));

    const propertyCell = this.buildPropertyCell(rowId, name);
    const weightCell = this.buildWeightCell(rowId, prop);
    const columnCell = this.buildColumnSelector(rowId, name, editor);

    const modeGear = this.design ?
      this.buildModeGear(rowId, prop, editor) :
      null;

    const controls = this.design ? this.buildRowControls(rowId) : null;

    row.append(
      ui.divV([propertyCell, columnCell].filter(Boolean)),
      weightCell,
      ui.divH([editor.root, modeGear].filter(Boolean)),
    );

    if (controls)
      row.append(controls);

    return row;
  }

  private buildPropertyCell(rowId: string, name: string): HTMLElement | null {
    if (!this.design)
      return ui.divText(name, 'statistics-mpo-property-name');

    if (!this.dataFrame) {
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

    return null;
  }

  private buildWeightCell(
    rowId: string,
    prop: PropertyDesirability,
  ): HTMLElement {
    const name = this.getPropertyNameByRowId(rowId);
    const weightInput = ui.input.float('', {value: prop.weight, min: 0, max: 1, format: '#0.000',
      onValueChanged: (v) => {
        if (name)
          this.mutateProperty(name, (p) => p.weight = Math.max(0, Math.min(1, v ?? 0)));
      },
    });

    weightInput.root.classList.add(
      'statistics-mpo-weight-input',
      this.design ?
        'statistics-mpo-weight-design' :
        'statistics-mpo-weight-view',
    );
    return weightInput.root;
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
      const col = this.dataFrame.col(matched) ?? null;
      editor.setColumn?.(col);
    }

    const input = ui.input.choice('', {items, nullable: false, value: matched ?? '', onValueChanged: (v) => {
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

  private resolveColumn(name: string): DG.Column | null {
    if (!this.dataFrame)
      return null;
    const eligible = this.getEligibleColumnNames();
    const colName = this.columnMapping[name] ??
      eligible.find((c) => c.toLowerCase() === name.toLowerCase()) ??
      eligible[0] ?? null;
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
    });

    gear.classList.add('statistics-mpo-gear');
    return gear;
  }

  private buildRowControls(rowId: string): HTMLElement {
    const add = ui.icons.add(() => this.insertRowAfterRow(rowId));
    const del = ui.icons.delete(() => this.deleteRow(rowId));
    return ui.divH([add, del], 'statistics-mpo-control-buttons');
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
        this.mutateProperty(name, (p) => Object.assign(p, patch));
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

    this.rows[rowId]?.remove();
    delete this.rows[rowId];

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

  private emitChange(): void {
    this.onChanged.next();
    grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
  }
}
