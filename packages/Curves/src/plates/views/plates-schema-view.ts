import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateProperty, PlateTemplate, plateTemplates, plateTypes, savePlate} from '../plates-crud';
import {SchemaEditor} from '@datagrok-libraries/utils/src/schema-editor';
import {PlateWidget} from '../../plate/plate-widget';


function platePropertyToOptions(p: PlateProperty): DG.PropertyOptions {
  return { name: p.name, type: p.value_type}
}

export function propertySchemaView(template: PlateTemplate): DG.View {
  const view = DG.View.create();
  view.name = 'Templates / ' + template.name;

  //@ts-ignore
  const platePropEditor = new SchemaEditor({properties: template.plateProperties.map(platePropertyToOptions)});
  //@ts-ignore
  const wellPropEditor = new SchemaEditor({properties: template.wellProperties.map(platePropertyToOptions)});

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
