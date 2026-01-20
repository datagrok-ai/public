import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {DesirabilityProfile, PropertyDesirability} from './mpo';
import {MpoDesirabilityLineEditor} from './mpo-line-editor';
import {Subject} from 'rxjs';

import '../../css/styles.css';

export const MPO_SCORE_CHANGED_EVENT = 'grok-mpo-score-changed';

export class MpoProfileEditor {
  readonly root = ui.div([]);
  readonly onChanged = new Subject<void>();

  dataFrame?: DG.DataFrame;
  design = false;
  profile?: DesirabilityProfile;
  columnMapping: Record<string, string | null> = {};

  constructor(dataFrame?: DG.DataFrame, design = false) {
    this.dataFrame = dataFrame;
    this.design = design;
  }

  getProfile(): DesirabilityProfile | undefined {
    return this.profile;
  }

  setDesignMode(on: boolean): void {
    if (this.design === on)
      return;

    this.design = on;
    if (this.profile)
      this.setProfile(this.profile);
  }

  setProfile(profile?: DesirabilityProfile): void {
    this.profile = profile;
    this.columnMapping = {};
    ui.empty(this.root);

    if (!profile) {
      this.root.append(ui.divText('No profile specified.'));
      return;
    }

    const header = this.buildHeader();
    const rows = Object.entries(profile.properties).map(([name, prop]) => this.buildRow(name, prop));

    if (rows.length === 0) {
      this.root.append(ui.divText('No properties defined.'));
      return;
    }

    this.root.append(header, ui.divV(rows));
  }

  private buildRow(propertyName: string, prop: PropertyDesirability): HTMLElement {
    const row = ui.divH([], 'statistics-mpo-row');

    const lineEditor = new MpoDesirabilityLineEditor(prop, 300, 80);
    lineEditor.onChanged.subscribe(() => this.notifyChanged());

    const propertyCell = this.buildPropertyCell(propertyName);
    const weightCell = this.buildWeightAndRangeCell(propertyName, prop, lineEditor);
    const columnCell = this.buildColumnSelector(propertyName, lineEditor);
    const controls = this.design ? this.buildRowControls(row, propertyName) : null;

    row.append(
      ui.divV([propertyCell, columnCell].filter(Boolean)),
      weightCell,
      lineEditor.root,
    );

    if (controls)
      row.append(controls);

    return row;
  }

  private buildPropertyCell(propertyName: string): HTMLElement | null {
    if (!this.design)
      return ui.divText(propertyName, 'statistics-mpo-property-name');

    if (this.dataFrame)
      return null;

    return ui.input.string('', {
      value: propertyName,
      onValueChanged: (v) => {
        if (v && v !== propertyName)
          this.renameProperty(propertyName, v);
      },
    }).root;
  }

  private buildWeightAndRangeCell(
    propertyName: string,
    prop: PropertyDesirability,
    lineEditor: MpoDesirabilityLineEditor,
  ): HTMLElement {
    const weightInput = ui.input.float('', {
      value: prop.weight,
      min: 0,
      max: 1,
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.weight = Math.max(0, Math.min(1, v ?? 0));
        });
      },
    });

    weightInput.root.classList.add(
      'statistics-mpo-weight-input',
      this.design ? 'statistics-mpo-weight-design' : 'statistics-mpo-weight-view',
    );


    if (!this.design)
      return weightInput.root;

    const minInput = ui.input.float('', {
      value: prop.min ?? lineEditor.getMinX(),
      showSlider: true,
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.min = v ?? 0;
          lineEditor.redrawAll();
        });
      },
    }).root;

    const maxInput = ui.input.float('', {
      value: prop.max ?? lineEditor.getMaxX(),
      showSlider: true,
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.max = v ?? 1;
          lineEditor.redrawAll();
        });
      },
    }).root;

    minInput.classList.add('statistics-mpo-range-input');
    maxInput.classList.add('statistics-mpo-range-input');

    return ui.divH(
      [weightInput.root, ui.divH([minInput, maxInput])],
    );
  }

  private buildColumnSelector(
    propertyName: string,
    lineEditor: MpoDesirabilityLineEditor,
  ): HTMLElement | null {
    if (!this.dataFrame)
      return null;

    const columnNames = this.dataFrame.columns.names();
    const matched = columnNames.find(
      (c) => c.toLowerCase() === propertyName.toLowerCase(),
    ) ?? null;

    const drawFromColumn = (name: string | null) => {
      if (!name)
        return;

      const values = this.dataFrame!.col(name)?.toList();
      lineEditor.drawBars(values);
    };

    drawFromColumn(matched);

    const input = ui.input.choice('', {
      items: columnNames,
      nullable: true,
      value: matched ?? '',
      onValueChanged: (v) => {
        this.columnMapping[propertyName] = v ?? null;
        drawFromColumn(v ?? null);
      },
    });

    input.root.style.width = '100px';
    return input.root;
  }

  private buildRowControls(row: HTMLElement, propertyName: string): HTMLElement {
    const add = ui.icons.add(() => this.insertRowAfter(row));
    const del = ui.icons.delete(() => this.deleteProperty(propertyName));
    return ui.divH([add, del], 'statistics-mpo-control-buttons');
  }

  private buildHeader(): HTMLElement {
    return ui.divH([
      ui.divText('Property', 'statistics-mpo-header-property'),
      ui.divText('Weight', 'statistics-mpo-header-weight'),
      this.design ? ui.divText('Min', 'statistics-mpo-header-weight') : null,
      this.design ? ui.divText('Max', 'statistics-mpo-header-weight') : null,
      ui.divText('Desirability', 'statistics-mpo-header-desirability'),
    ].filter(Boolean), 'statistics-mpo-header');
  }

  private updateProperty(
    name: string,
    updater: (p: PropertyDesirability) => void,
  ): void {
    const prop = this.profile?.properties[name];
    if (!prop)
      return;

    updater(prop);
    this.notifyChanged();
  }

  private renameProperty(oldName: string, newName: string): void {
    if (!this.profile || this.profile.properties[newName])
      return;

    this.profile.properties[newName] = this.profile.properties[oldName];
    delete this.profile.properties[oldName];

    this.columnMapping[newName] = this.columnMapping[oldName] ?? null;
    delete this.columnMapping[oldName];

    this.setProfile(this.profile);
    this.notifyChanged();
  }

  private deleteProperty(name: string): void {
    if (!this.profile)
      return;

    delete this.profile.properties[name];
    delete this.columnMapping[name];

    this.setProfile(this.profile);
    this.notifyChanged();
  }

  private insertRowAfter(_: HTMLElement): void {
    if (!this.profile)
      return;

    const name = `NewProperty${Object.keys(this.profile.properties).length + 1}`;
    this.profile.properties[name] = {weight: 1, min: 0, max: 1, line: []};

    this.setProfile(this.profile);
    this.notifyChanged();
  }

  private notifyChanged(): void {
    this.onChanged.next();
    grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
  }
}
