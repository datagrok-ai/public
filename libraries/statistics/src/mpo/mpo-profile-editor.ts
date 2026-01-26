import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subject} from 'rxjs';
import {DesirabilityProfile, PropertyDesirability} from './mpo';
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
  columnMapping: Record<string, string | null> = {};

  constructor(dataFrame?: DG.DataFrame, design = false) {
    this.dataFrame = dataFrame;
    this.design = design;
  }

  setProfile(profile?: DesirabilityProfile): void {
    this.profile = profile;
    this.columnMapping = {};
    this.render();
  }

  getProfile(): DesirabilityProfile | undefined {
    return this.profile;
  }

  setDesignMode(on: boolean): void {
    if (this.design === on)
      return;
    this.design = on;
    this.render();
  }

  private render(): void {
    ui.empty(this.root);

    if (!this.profile)
      return this.renderEmpty('No profile specified.');

    const rows = Object.entries(this.profile.properties)
      .map(([name, prop]) => this.buildRow(name, prop));

    if (!rows.length)
      return this.renderEmpty('No properties defined.');

    this.root.append(
      this.buildHeader(),
      ui.divV(rows),
    );
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

    prop.mode ??= 'freeform';
    const editor = DesirabilityEditorFactory.create(prop);
    editor.onChanged.subscribe(() => this.emitChange());

    const propertyCell = this.buildPropertyCell(name);
    const weightCell = this.buildWeightAndRangeCell(name, prop);
    const columnCell = this.buildColumnSelector(name, editor);

    const modeGear = this.design ?
      this.buildModeGear(name, prop, editor) :
      null;

    row.append(
      ui.divV([propertyCell, columnCell].filter(Boolean)),
      weightCell,
      ui.divH([editor.root, modeGear].filter(Boolean)),
    );

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

  private buildWeightAndRangeCell(
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

    const draw = (col: string | null) => {
      if (!col) return;
      const values = this.dataFrame!.col(col)?.toList();
      (editor as any).drawBars?.(values);
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

    input.root.style.width = '100px';
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

  private openModeDialog(
    name: string,
    prop: PropertyDesirability,
    editor: DesirabilityEditor,
  ): void {
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
    if (!prop) return;

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

    this.render();
    this.emitChange();
  }

  private emitChange(): void {
    this.onChanged.next();
    grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
  }
}
