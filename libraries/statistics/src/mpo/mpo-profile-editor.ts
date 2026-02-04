import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subject} from 'rxjs';
import {DesirabilityProfile, PropertyDesirability, isNumerical, migrateDesirability} from './mpo';
import {DesirabilityEditor, DesirabilityEditorFactory} from './editors/desirability-editor-factory';
import {DesirabilityModeDialog} from './dialogs/desirability-mode-dialog';

import '../../css/styles.css';

export const MPO_SCORE_CHANGED_EVENT = 'grok-mpo-score-changed';

export class MpoProfileEditor {
  readonly root = ui.div([]);
  readonly onChanged = new Subject<void>();

  profile?: DesirabilityProfile;
  dataFrame?: DG.DataFrame;
  design = false;
  preview = false;

  private rows: Record<string, HTMLElement> = {};
  private rowIds: Record<string, string> = {};

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
      for (const key of Object.keys(profile.properties))
        profile.properties[key] = migrateDesirability(profile.properties[key]);
    }
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

    if (isNumerical(prop))
      prop.mode ??= 'freeform';
    const editor = DesirabilityEditorFactory.create(prop);
    editor.onChanged.subscribe(() => this.emitChange());

    const propertyCell = this.buildPropertyCell(rowId, name);
    const weightCell = this.buildWeightCell(rowId, prop);
    const columnCell = this.buildColumnSelector(rowId, editor);

    const modeGear = this.design && editor.supportsModeDialog ?
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

      const propNameInp = ui.input.string('', {
        value: name,
        onValueChanged: (v) => {
          if (!v || v === currentName)
            return;
          this.renameProperty(currentName, v);
          currentName = v;
        },
      });

      ui.tooltip.bind(propNameInp.input, () => currentName);
      return propNameInp.root;
    }

    return null;
  }

  private buildWeightCell(
    rowId: string,
    prop: PropertyDesirability,
  ): HTMLElement {
    const name = this.getPropertyNameByRowId(rowId)!;
    const weightInput = ui.input.float('', {
      value: prop.weight,
      min: 0,
      max: 1,
      format: '#0.000',
      onValueChanged: (v) =>
        this.mutateProperty(name, (p) =>
          p.weight = Math.max(0, Math.min(1, v ?? 0))),
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
    editor: DesirabilityEditor,
  ): HTMLElement | null {
    if (!this.dataFrame)
      return null;

    const name = this.getPropertyNameByRowId(rowId)!;
    const columns = this.dataFrame.columns.names();
    const matched =
      columns.find((c) => c.toLowerCase() === name.toLowerCase()) ?? null;

    const draw = (colName: string | null) => {
      const col = colName ? this.dataFrame!.col(colName) : null;
      editor.setColumn?.(col);
    };

    draw(matched);

    const input = ui.input.choice('', {
      items: columns,
      nullable: true,
      value: matched ?? '',
      onValueChanged: (v) => {
        this.columnMapping[name] = v ?? null;
        draw(v ?? null);
        this.emitChange();
      },
    });
    return input.root;
  }

  private buildModeGear(
    rowId: string,
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
  ): HTMLElement {
    const name = this.getPropertyNameByRowId(rowId)!;
    const gear = ui.icons.settings(() =>
      this.openModeDialog(name, prop, editor));

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
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
  ): void {
    if (!isNumerical(prop))
      return;

    new DesirabilityModeDialog(
      name,
      prop,
      (patch) => {
        this.mutateProperty(name, (p) => Object.assign(p, patch));
        editor.redrawAll();
      },
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

    this.profile.properties[newName] = {
      functionType: 'numerical',
      weight: 1,
      mode: 'freeform',
      min: 0,
      max: 1,
      line: [],
    };

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
