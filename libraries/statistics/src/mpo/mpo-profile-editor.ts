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
  private propertyOrder: string[] = [];
  columnMapping: Record<string, string | null> = {};

  constructor(dataFrame?: DG.DataFrame, design = false, preview = false) {
    this.dataFrame = dataFrame;
    this.design = design;
    this.preview = preview;
  }

  setProfile(profile?: DesirabilityProfile): void {
    if (profile) {
      for (const key of Object.keys(profile.properties))
        profile.properties[key] = migrateDesirability(profile.properties[key]);
    }
    this.profile = profile;
    this.columnMapping = {};
    this.rows = {};
    this.propertyOrder = profile ? Object.keys(profile.properties) : [];
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

    const rows = this.propertyOrder
      .map((name) => this.profile!.properties[name])
      .filter(Boolean)
      .map((prop, i) => {
        const name = this.propertyOrder[i];
        if (!this.rows[name])
          this.rows[name] = this.buildRow(name, prop);
        return this.rows[name];
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
    ].filter(Boolean), 'statistics-mpo-header');
  }

  private buildRow(
    name: string,
    prop: PropertyDesirability,
  ): HTMLElement {
    const row = ui.divH([], 'statistics-mpo-row');

    if (isNumerical(prop))
      prop.mode ??= 'freeform';
    const editor = DesirabilityEditorFactory.create(prop);
    editor.onChanged.subscribe(() => this.emitChange());

    const propertyCell = this.buildPropertyCell(name);
    const weightCell = this.buildWeightCell(name, prop);
    const columnCell = this.buildColumnSelector(name, editor);

    const modeGear = this.design && editor.supportsModeDialog ?
      this.buildModeGear(name, prop, editor) :
      null;

    const controls = this.design ? this.buildRowControls(name) : null;

    row.append(
      ui.divV([propertyCell, columnCell].filter(Boolean)),
      weightCell,
      ui.divH([editor.root, modeGear].filter(Boolean)),
    );

    if (controls)
      row.append(controls);

    return row;
  }

  private buildPropertyCell(name: string): HTMLElement | null {
    if (!this.design)
      return ui.divText(name, 'statistics-mpo-property-name');

    if (!this.dataFrame) {
      return ui.input.string('', {
        value: name,
        onValueChanged: (v) => {
          if (v && v !== name)
            this.renameProperty(name, v);
        },
      }).root;
    }

    return null;
  }

  private buildWeightCell(
    name: string,
    prop: PropertyDesirability,
  ): HTMLElement {
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
    name: string,
    editor: DesirabilityEditor,
  ): HTMLElement | null {
    if (!this.dataFrame)
      return null;

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
      },
    });
    return input.root;
  }

  private buildModeGear(
    name: string,
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
  ): HTMLElement {
    const gear = ui.icons.settings(() =>
      this.openModeDialog(name, prop, editor));

    gear.classList.add('statistics-mpo-gear');
    return gear;
  }

  private buildRowControls(name: string): HTMLElement {
    const add = ui.icons.add(() => this.insertRowAfter(name));
    const del = ui.icons.delete(() => this.deleteProperty(name));
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

  private renameProperty(oldName: string, newName: string): void {
    if (!this.profile || this.profile.properties[newName])
      return;

    this.profile.properties[newName] = this.profile.properties[oldName];
    delete this.profile.properties[oldName];

    this.columnMapping[newName] = this.columnMapping[oldName] ?? null;
    delete this.columnMapping[oldName];

    const idx = this.propertyOrder.indexOf(oldName);
    this.propertyOrder[idx] = newName;

    this.rows[newName] = this.rows[oldName];
    delete this.rows[oldName];

    this.emitChange();
  }

  private deleteProperty(name: string): void {
    if (!this.profile)
      return;

    delete this.profile.properties[name];
    delete this.columnMapping[name];

    const idx = this.propertyOrder.indexOf(name);
    if (idx >= 0)
      this.propertyOrder.splice(idx, 1);

    this.rows[name]?.remove();
    delete this.rows[name];

    this.emitChange();
  }

  private insertRowAfter(propertyName: string): void {
    if (!this.profile)
      return;

    const newName = `NewProperty${Object.keys(this.profile.properties).length + 1}`;
    this.profile.properties[newName] = {functionType: 'numerical', weight: 1, mode: 'freeform', min: 0, max: 1, line: []};

    const idx = this.propertyOrder.indexOf(propertyName);
    this.propertyOrder.splice(idx + 1, 0, newName);

    const newRow = this.buildRow(newName, this.profile.properties[newName]);
    this.rows[newName] = newRow;

    this.rows[propertyName]?.after(newRow);

    this.emitChange();
  }

  private emitChange(): void {
    this.onChanged.next();
    grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
  }
}
