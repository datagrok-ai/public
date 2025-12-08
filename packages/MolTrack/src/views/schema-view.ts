import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {SchemaEditor} from '@datagrok-libraries/utils/src/schema-editor';
import {merge} from 'rxjs';

import {fetchSchema, registerMolTrackProperties} from '../package';
import {GroupedProperties, MolTrackProp} from '../utils/constants';
import {buildPropertyOptions} from '../utils/utils';

let openedView: DG.ViewBase | null = null;

export class PropertySchemaView {
  view: DG.View;
  extraPropertiesDiv: HTMLDivElement;
  grouped: GroupedProperties = {};
  editors: Record<string, SchemaEditor> = {};
  saveButton: HTMLButtonElement;

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Register a schema';
    this.extraPropertiesDiv = ui.div([], 'ui-form');
    this.saveButton = ui.bigButton('SAVE', () => this.onSave());
    this.saveButton.disabled = true;
  }

  async init(): Promise<DG.View> {
    const parsed: any = JSON.parse(await fetchSchema());
    const propArray = parsed.properties ?? parsed;

    this.grouped = this.groupByEntityType(propArray);
    this.render();
    this.view.setRibbonPanels([[this.saveButton]]);

    this.initSubscriptions();
    return this.view;
  }

  private initSubscriptions(): void {
    const editorsArray = Object.values(this.editors);

    merge(...editorsArray.map((ed) => ed.table.onChanged))
      .subscribe(() => this.saveButton.disabled = false);

    merge(...editorsArray.map((ed) => DG.debounce(ed.table.onSelected, 5)))
      .subscribe((item) => {
        if (!item.userEditable) this.disableInputsIn(this.extraPropertiesDiv);
      });

    merge(...editorsArray.map((ed) => ed.table.onItemAdded))
      .subscribe((item) => item.userEditable = true);
  }

  private disableInputsIn(element: HTMLElement): void {
    element.querySelectorAll('.ui-input-editor').forEach((el) => {
      if (el instanceof HTMLSelectElement)
        el.disabled = true;
      else if (el instanceof HTMLInputElement)
        el.type === 'checkbox' ? el.disabled = true : el.readOnly = true;
    });
  }

  private disableAllInputs(): void {
    for (const editor of Object.values(this.editors))
      this.disableInputsIn(editor.root);

    this.disableInputsIn(this.extraPropertiesDiv);
  }

  private async onSave(): Promise<void> {
    this.saveButton.disabled = true;
    try {
      const allProps = [];

      for (const [entityType, editor] of Object.entries(this.editors)) {
        for (const p of editor.properties) {
          if (!p || Object.keys(p).length === 0)
            continue;

          p.userEditable = false;
          const newProp = {
            ...p,
            value_type: p.type,
            entity_type: entityType,
            property_class: 'measured',
            unit: p.units,
            friendly_name: p.friendlyName,
          };
          allProps.push(newProp);
        }
      }

      await registerMolTrackProperties(JSON.stringify({properties: allProps}));
      grok.shell.info('Schema saved');

      this.disableAllInputs();
    } catch (e: any) {
      grok.shell.error('Failed to save');
      console.error(e);
    } finally {
      this.saveButton.disabled = true;
    }
  }

  private groupByEntityType(propArray: any[]): GroupedProperties {
    return propArray.reduce((acc: GroupedProperties, p: MolTrackProp) => {
      const entityType = p.entity_type ?? 'Unknown';
      (acc[entityType] ??= []).push(buildPropertyOptions(p));
      return acc;
    }, {});
  }

  private render(): void {
    ui.empty(this.view.root);

    const accordion = ui.accordion();

    for (const [entityType, props] of Object.entries(this.grouped)) {
      const editor = new SchemaEditor({
        properties: props,
        extraPropertiesDiv: this.extraPropertiesDiv,
        allowRemove: false,
      });
      editor.extraProperties = ['nullable', 'description', 'units', 'min', 'max', 'friendlyName'];
      this.editors[entityType] = editor;

      accordion.addPane(entityType, () => editor.root, true);
    }

    this.view.root.appendChild(
      ui.divH([
        ui.divV([
          ui.h2('MolTrack properties'),
          accordion.root,
        ]),
        ui.div([], {style: {width: '20px'}}),
        this.extraPropertiesDiv,
      ]),
    );

    this.disableAllInputs();
  }

  show(): void {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
