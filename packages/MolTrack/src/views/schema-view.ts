import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { SchemaEditor } from '@datagrok-libraries/utils/src/schema-editor';

import { fetchSchema } from '../package';
import { MolTrackProp } from '../utils/constants';
import { buildPropertyOptions } from '../utils/utils';


type GroupedProperties = Record<string, any[]>;
let openedView: DG.ViewBase | null = null;

export class PropertySchemaView {
  view: DG.View;
  extraPropertiesDiv: HTMLDivElement;
  grouped: GroupedProperties = {};

  constructor() {
    this.view = DG.View.create();
    this.extraPropertiesDiv = ui.div([]);
  }

  async init(): Promise<DG.View> {
    const parsed: any = JSON.parse(await fetchSchema());
    const propArray = parsed.properties ?? parsed;

    this.grouped = this.groupByEntityType(propArray);
    this.render();

    return this.view;
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
      accordion.addPane(entityType, () => {
        return new SchemaEditor({
          properties: props,
          extraPropertiesDiv: this.extraPropertiesDiv,
        }).root;
      }, true);
    }

    this.view.root.appendChild(ui.divH([
      ui.divV([
        ui.h2('MolTrack properties'),
        accordion.root,
      ]),
      ui.div([], { style: { width: '20px' } }),
      this.extraPropertiesDiv,
    ]));
  }

  show(): void {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
