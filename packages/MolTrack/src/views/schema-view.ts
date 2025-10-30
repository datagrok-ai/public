import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { SchemaEditor } from '@datagrok-libraries/utils/src/schema-editor';
import { merge } from 'rxjs';

import { fetchSchema, registerMolTrackProperties } from '../package';
import { GroupedProperties, MolTrackProp } from '../utils/constants';
import { buildPropertyOptions } from '../utils/utils';

let openedView: DG.ViewBase | null = null;

export class PropertySchemaView {
  view: DG.View;
  extraPropertiesDiv: HTMLDivElement;
  grouped: GroupedProperties = {};
  editors: Record<string, SchemaEditor> = {};
  saveButton: HTMLButtonElement;

  constructor() {
    this.view = DG.View.create();
    this.extraPropertiesDiv = ui.div([]);
    this.saveButton = ui.bigButton('SAVE', () => this.onSave());
    this.saveButton.disabled = true;
  }

  async init(): Promise<DG.View> {
    const parsed: any = JSON.parse(await fetchSchema());
    const propArray = parsed.properties ?? parsed;

    this.grouped = this.groupByEntityType(propArray);
    this.render();
    this.view.setRibbonPanels([[this.saveButton]]);

    this.subscribeToEditorChanges();
    return this.view;
  }

  private subscribeToEditorChanges(): void {
    const editorChanges = Object.values(this.editors).map((ed) => ed.table.onChanged);
    merge(...editorChanges).subscribe(() => this.saveButton.disabled = false);
  }

  private async onSave(): Promise<void> {
    this.saveButton.disabled = true;
    try {
      const allProps = [];

      for (const [entityType, editor] of Object.entries(this.editors)) {
        for (const p of editor.properties) {
          if (!p || Object.keys(p).length === 0)
            continue;

          const newProp = { ...p, value_type: p.type, entity_type: entityType, property_class: 'measured' };
          delete newProp.type;
          allProps.push(newProp);
        }
      }

      const jsonObject = {'properties': allProps};
      const jsonPayload = JSON.stringify(jsonObject);
      await registerMolTrackProperties(jsonPayload);
      grok.shell.info('Schema saved');
    } catch (e: any) {
      grok.shell.error('Failed to save');
      console.log(e.message ?? e);
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
      });
      this.editors[entityType] = editor;

      accordion.addPane(entityType, () => editor.root, true);
    }

    this.view.root.appendChild(
      ui.divH([
        ui.divV([
          ui.h2('MolTrack properties'),
          accordion.root,
        ]),
        ui.div([], { style: { width: '20px' } }),
        this.extraPropertiesDiv,
      ]),
    );
  }

  show(): void {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
