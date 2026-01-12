import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {DesirabilityProfile, PropertyDesirability} from './mpo';
import {MpoDesirabilityLineEditor} from './mpo-line-editor';
import {Subject} from 'rxjs';

import '../../css/styles.css';

export class MpoProfileEditor {
  root = ui.div([]);
  dataFrame?: DG.DataFrame;
  design: boolean = false;
  onChanged = new Subject();
  profile?: DesirabilityProfile;

  columnMapping: Record<string, string | null> = {};

  constructor(dataFrame?: DG.DataFrame, design: boolean = false) {
    this.dataFrame = dataFrame;
    this.design = design;
    this.setProfile();
  }

  private buildRow(propertyName: string, prop: PropertyDesirability): HTMLElement {
    const row = ui.divH([], 'statistics-mpo-row');
    const lineEditor = new MpoDesirabilityLineEditor(prop, 300, 80);
    lineEditor.onChanged.subscribe((_) => this.onChanged.next());

    // Input for weight - updates the *copy* of the template data
    const weightInput = ui.input.float('', { // No label needed here
      value: prop.weight, min: 0, max: 1,
      onValueChanged: (newValue) => { // Changed parameter name for clarity
        if (this.profile && this.profile.properties[propertyName]) {
          // Update the weight in the temporary template object
          let clampedWeight = newValue ?? 0;
          clampedWeight = Math.max(0, Math.min(1, clampedWeight)); // Clamp 0-1
          this.profile.properties[propertyName].weight = clampedWeight;
        }
      },
    });
    weightInput.root.classList.add('statistics-mpo-weight-input');

    const matchedColumnName = this.dataFrame ?
      this.dataFrame.columns.names().find((name) => name.toLowerCase() == propertyName.toLowerCase()) :
      null;

    const drawBarsFromColumnValues = (columnName: string | null | undefined) => {
      if (!this.dataFrame || !columnName)
        return null;
      const values = this.dataFrame.col(columnName)?.toList();
      lineEditor.drawBars(values);
    };

    drawBarsFromColumnValues(matchedColumnName);
    const columnInput = ui.input.choice('', {
      nullable: true,
      items: this.dataFrame?.columns?.names() ?? [''],
      value: matchedColumnName ?? '',
      onValueChanged: (value) => {
        drawBarsFromColumnValues(value);
        this.columnMapping[propertyName] = value ?? null;
      },
    });

    const controls = this.design ? (() => {
      const add = ui.icons.add(() => this.insertRowAfter(row));
      const del = ui.icons.delete(() => this.deleteRow(row, propertyName));
      return ui.divH([add, del], 'statistics-mpo-control-buttons');
    })() : null;

    const propertyLabel = this.design ?
      ui.input.string('', {
        value: propertyName,
        onValueChanged: (v) => {
          if (!v || !this.profile) return;
          if (v === propertyName) return;

          this.profile.properties[v] = this.profile.properties[propertyName];
          delete this.profile.properties[propertyName];
        },
      }).root :
      ui.divText(propertyName, 'statistics-mpo-property-name');


    row.append(
      ui.divV([
        propertyLabel,
        this.dataFrame ? columnInput.root : null,
      ]),
      weightInput.root,
      lineEditor.root, // Add the Konva container div
    );

    if (controls)
      row.append(controls);

    return row;
  }

  insertRowAfter(afterRow: HTMLElement) {
    if (!this.profile) return;
    const name = ``;
    const prop: PropertyDesirability = {weight: 0.5, line: []};
    this.profile.properties[name] = prop;
    afterRow.after(this.buildRow(name, prop));
  }

  deleteRow(row: HTMLElement, name: string) {
    delete this.profile?.properties[name];
    row.remove();
  }

  getProfile(): DesirabilityProfile | undefined {
    return this.profile;
  }

  setDesignMode(on: boolean) {
    this.design = on;
    if (this.profile)
      this.setProfile(this.profile);
  }

  setProfile(profile?: DesirabilityProfile) {
    this.profile = profile;
    this.columnMapping = {};
    ui.empty(this.root);
    if (!profile) {
      this.root.append(ui.divText('No profile specified.'));
      return;
    }

    // Create header row for the table
    const header = ui.divH([
      ui.divText('Property', 'statistics-mpo-header-property'),
      ui.divText('Weight', 'statistics-mpo-header-weight'),
      ui.divText('Desirability', 'statistics-mpo-header-desirability'), // Let editor take space
    ], 'statistics-mpo-header');

    const propertyRows = Object.entries(profile.properties).map(([propertyName, prop]) => {
      return this.buildRow(propertyName, prop);
    }).filter((el) => el !== null); // Filter out skipped properties

    if (propertyRows.length > 0) {
      this.root.append(header);
      this.root.append(ui.divV(propertyRows as HTMLElement[])); // Cast needed after filter
    } else
      this.root.append(ui.divText('No matching properties found in the table for this template.'));
  }
}
