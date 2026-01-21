import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {DesirabilityMode, DesirabilityProfile, PropertyDesirability} from './mpo';
import {MpoDesirabilityLineEditor} from './mpo-line-editor';
import {Subject} from 'rxjs';

import '../../css/styles.css';
import {MpoCategoricalEditor} from './mpo-categorical-editor';

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

    // Ensure default mode exists
    Object.values(profile.properties).forEach((p) => {
      p.mode ??= 'freeform';
    });

    const header = this.buildHeader();
    const rows = Object.entries(profile.properties).map(([name, prop]) => this.buildRow(name, prop));

    if (rows.length === 0) {
      this.root.append(ui.divText('No properties defined.'));
      return;
    }

    this.root.append(header, ui.divV(rows));
  }

  private buildEditor(prop: PropertyDesirability): any {
    if (prop.mode === 'categorical')
      return new MpoCategoricalEditor(prop);
    return new MpoDesirabilityLineEditor(prop, 300, 80);
  }


  private buildRow(propertyName: string, prop: PropertyDesirability): HTMLElement {
    const row = ui.divH([], 'statistics-mpo-row');

    // const lineEditor = new MpoDesirabilityLineEditor(prop, 300, 80);
    const editor = this.buildEditor(prop);
    editor.onChanged.subscribe(() => this.notifyChanged());

    const propertyCell = this.buildPropertyCell(propertyName);
    const weightCell = this.buildWeightAndRangeCell(propertyName, prop, editor);
    const columnCell = this.buildColumnSelector(propertyName, editor);
    const modeControl = this.design ? this.buildModeGear(propertyName, prop, editor) : null;
    const controls = this.design ? this.buildRowControls(row, propertyName) : null;

    row.append(
      ui.divV([propertyCell, columnCell].filter(Boolean)),
      weightCell,
      ui.divV([modeControl, editor.root].filter(Boolean)),
    );

    if (controls)
      row.append(controls);

    return row;
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

  private buildModeGear(
    propertyName: string,
    prop: PropertyDesirability,
    lineEditor: MpoDesirabilityLineEditor,
  ): HTMLElement {
    const gear = ui.icons.settings(() => {
      this.openModeDialog(propertyName, prop, lineEditor);
    });

    gear.classList.add('statistics-mpo-gear');
    return gear;
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

  private ensureCategories(prop: PropertyDesirability, values: any[]) {
    const unique = Array.from(new Set(values)).map((v) => String(v));
    const existing = prop.categories ?? [];

    prop.categories = unique.map((name) => ({
      name,
      weight: existing.find((c) => c.name === name)?.weight ?? 0.5,
    }));
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
    this.profile.properties[name] = {
      weight: 1,
      min: 0,
      max: 1,
      line: [],
      mode: 'freeform',
    };

    this.setProfile(this.profile);
    this.notifyChanged();
  }

  private notifyChanged(): void {
    this.onChanged.next();
    grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
  }

  private openModeDialog(
    propertyName: string,
    prop: PropertyDesirability,
    mainLineEditor: MpoDesirabilityLineEditor,
  ) {
    const dialog = ui.dialog({title: `Desirability Settings: ${propertyName}`});

    const modeChoice = ui.input.choice('Mode', {
      items: ['freeform', 'gaussian', 'sigmoid', 'categorical'],
      value: prop.mode ?? 'freeform',
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.mode = (v ?? 'freeform') as DesirabilityMode;
        });
        updatePanel();
        dialogLineEditor.redrawAll();
        mainLineEditor.redrawAll();
      },
    });

    const meanInput = ui.input.float('Mean', {
      value: prop.mean ?? 0,
      //@ts-ignore
      format: '#0.000',
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.mean = v ?? 0;
        });
        dialogLineEditor.redrawAll();
        mainLineEditor.redrawAll();
      },
    });

    const sigmaInput = ui.input.float('Sigma', {
      value: prop.sigma ?? 1,
      //@ts-ignore
      format: '#0.000',
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.sigma = Math.max(0.01, v ?? 1);
        });
        dialogLineEditor.redrawAll();
        mainLineEditor.redrawAll();
      },
    });

    const x0Input = ui.input.float('x0', {
      value: prop.x0 ?? 0,
      //@ts-ignore
      format: '#0.000',
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.x0 = v ?? 0;
        });
        dialogLineEditor.redrawAll();
        mainLineEditor.redrawAll();
      },
    });

    const kInput = ui.input.float('k', {
      value: prop.k ?? 10,
      //@ts-ignore
      format: '#0.000',
      onValueChanged: (v) => {
        this.updateProperty(propertyName, (p) => {
          p.k = Math.max(0.1, v ?? 10);
        });
        dialogLineEditor.redrawAll();
        mainLineEditor.redrawAll();
      },
    });

    const paramPanel = ui.divV([]);

    const updatePanel = () => {
      ui.empty(paramPanel);

      if (prop.mode === 'freeform')
        return;

      paramPanel.append(ui.h3('Parameters'));

      if (prop.mode === 'gaussian')
        paramPanel.append(ui.form([meanInput, sigmaInput]));
      else
        paramPanel.append(ui.form([x0Input, kInput]));
    };

    updatePanel();

    const dialogLineEditor = new MpoDesirabilityLineEditor(prop, 300, 80);

    dialogLineEditor.onParamsChanged = (p) => {
      meanInput.value = p.mean ?? 0;
      sigmaInput.value = p.sigma ?? 1;
      x0Input.value = p.x0 ?? 0;
      kInput.value = p.k ?? 10;

      this.updateProperty(propertyName, (prop) => {
        Object.assign(prop, p);
      });

      mainLineEditor.redrawAll();
      this.notifyChanged();
    };


    dialog.add(
      ui.divV([
        modeChoice.root,
        paramPanel,
        dialogLineEditor.root,
      ]),
    );

    dialog.addButton('Close', () => dialog.close());
    dialog.show();
  }
}
