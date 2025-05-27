import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate, plateTemplates, plateTypes, savePlate} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';


export function propertySchemaView(template: PlateTemplate): DG.View {
  const view = DG.View.create();
  view.name = 'Templates / ' + template.name;

  view.root.appendChild(ui.divV([
    ui.h2('Plate properties'),
    ui.list(template.plateProperties.map(p => p.name)),
    ui.h2('Well properties'),
    ui.list(template.wellProperties.map(p => p.name)),
  ]));

  return view;
}
