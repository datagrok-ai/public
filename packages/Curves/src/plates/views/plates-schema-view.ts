import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateProperty, PlateTemplate, plateTemplates, plateTypes, savePlate} from '../plates-crud';
import {SchemaEditor} from '@datagrok-libraries/utils/src/schema-editor';
import {PlateWidget} from '../../plate/plate-widget';


function platePropertyToOptions(p: Partial<PlateProperty>): DG.IProperty {
  return { name: p.name, type: p.value_type}
}

export function propertySchemaView(template: PlateTemplate): DG.View {
  const view = DG.View.create();
  view.name = 'Templates / ' + template.name;

  const platePropEditor = new SchemaEditor({properties: template.plateProperties.map(platePropertyToOptions)});
  const wellPropEditor = new SchemaEditor({properties: template.wellProperties.map(platePropertyToOptions)});

  platePropEditor.table.onChanged.subscribe((_) => grok.shell.info('x'));
  platePropEditor.table.onItemChanged.subscribe((_) => grok.shell.info('changed'));
  platePropEditor.table.onItemAdded.subscribe((_) => grok.shell.info('added'));
  platePropEditor.table.onItemRemoved.subscribe((_) => grok.shell.info('removed'));

  view.root.appendChild(ui.divV([
    ui.h2('Plate properties'),
    //ui.list(template.plateProperties.map(p => p.name)),
    platePropEditor.root,

    ui.h2('Well properties'),
    //ui.list(template.wellProperties.map(p => p.name)),
    wellPropEditor.root,
  ]));

  return view;
}
