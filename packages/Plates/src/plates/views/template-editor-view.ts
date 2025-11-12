// src/plates/views/plates-schema-view.ts

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {PlateTemplate} from '../plates-crud';
// import {PlateWidget} from '../../plate/plate-widget';
import {Plate} from '../../plate/plate';
import {PlateWidget} from '../../plate/plate-widget/plate-widget';

/**
 * Creates a view for editing a template's properties and
 * displaying a prototype PlateWidget in the adjacent panel.
 */
export function propertySchemaView(template: PlateTemplate): DG.View {
  const view = DG.View.create();
  view.name = template.id === -1 ? 'New Template' : `Edit: ${template.name}`;
  view.box = true; // Adds a bit of padding around the view

  // --- Left Panel: Template Property Editor ---
  const nameInput = ui.input.string('Name', {value: template.name});
  const descInput = ui.input.string('Description', {value: template.description});

  // Simple display for properties (not editable, as per request)
  const platePropsList = ui.divV(
    template.plateProperties.map((p) => ui.divText(p.name ?? ''))
  );
  const wellPropsList = ui.divV(
    template.wellProperties.map((p) => ui.divText(p.name ?? ''))
  );

  const leftPanel = ui.divV([
    nameInput.root,
    descInput.root,
    ui.h3('Plate properties'),
    platePropsList,
    ui.h3('Well properties'),
    wellPropsList,
  ], 'ui-form');


  // --- Right Panel: PlateWidget Sandbox ---
  const plateWidget = new PlateWidget();
  // Provide a default 96-well plate for the widget to display
  plateWidget.plate = new Plate(8, 12);
  // **Enable editing to allow interaction and prototyping**
  plateWidget.editable = true;

  const rightPanel = ui.divV([
    plateWidget.root
  ], {style: {width: '100%', height: '100%'}});

  // Set the widget to take up all available space in its panel
  plateWidget.root.style.width = '100%';
  plateWidget.root.style.height = '100%';

  // --- Final Layout ---
  view.root.appendChild(ui.splitH([leftPanel, rightPanel]));

  view.setRibbonPanels([[
    ui.bigButton('SAVE', () => {
      // This button is just a placeholder for now
    //   grok.shell.info('Save functionality is not connected in this prototype view.');
    }),
  ]]);

  return view;
}
