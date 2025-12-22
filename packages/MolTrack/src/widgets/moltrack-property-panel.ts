/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {EntityBaseView} from '../views/registration-entity-base';
import {compoundView} from '../utils/view-utils';
import {getBatchInfoBySynonym} from '../utils/utils';

function extractPropValue(prop: Record<string, any>): any {
  return prop?.value_uuid ?? prop?.value_num ?? prop?.value_string ?? prop?.value_qualifier ?? prop?.value_datetime ?? null;
}

function addPlusIconToValue(column: DG.Column, key: string, value: any): HTMLDivElement {
  const div = ui.divH([ui.divText(value ?? '')]);

  const addColumnIcon = ui.iconFA(
    'plus',
    async () => {
      try {
        const batchData = await getBatchInfoBySynonym('compounds.details.corporate_compound_id', column.toList(), [key]);
        const {dataFrame} = column;
        const unusedColumnName = dataFrame.columns.getUnusedName(key);
        const newColumn = DG.Column.fromType(DG.TYPE.OBJECT, unusedColumnName, column.length);
        dataFrame.columns.add(newColumn);
        newColumn.init((i) => batchData[i][key.toLowerCase()] ?? null);
      } catch (err) {
        console.error(err);
      }
    },
    `Calculate ${key} for the whole table`,
  );

  div.prepend(addColumnIcon);
  ui.tools.setHoverVisibility(div, [addColumnIcon]);
  return div;
}

function createTableFromObject(column: DG.Column, obj: Record<string, unknown>, excludeKeys: string[] = []) {
  const map: Record<string, unknown> = {};
  for (const [key, value] of Object.entries(obj)) {
    if (!excludeKeys.includes(key))
      map[key] = addPlusIconToValue(column, `compounds.${key}`, value);
  }
  return ui.tableFromMap(map);
}

function createPropertiesTable(column: DG.Column, properties: any[]) {
  if (!properties?.length) return ui.divText('No properties');
  const map: Record<string, unknown> = {};
  for (const prop of properties) {
    const value = extractPropValue(prop);
    map[prop.name] = addPlusIconToValue(column, `compounds.details.${prop.name}`, value);
  }
  return ui.tableFromMap(map);
}

export async function molTrackPropPanel(
  retrievedCompound: any,
  column: DG.Column,
  initialStructure?: string,
): Promise<DG.Widget> {
  const {detail, properties = []} = retrievedCompound;

  if (detail) {
    const registerButton = ui.bigButton('Register', async () => {
      const registrationView = new EntityBaseView(false);
      if (initialStructure) registrationView.initialSmiles = initialStructure;
      await registrationView.buildUIMethod();
      registrationView.show();
    });
    registerButton.classList.add('moltrack-panel-register-button');

    const div = ui.divV([ui.divText('No matches'), registerButton], 'moltrack-panel-register-container');
    return DG.Widget.fromRoot(div);
  }

  const acc = ui.accordion();
  const corpProp = properties.find((p: any) => p.name === 'corporate_compound_id');
  const corpValue = extractPropValue(corpProp);
  const corpLink = ui.link(corpValue?.toString() ?? '', async () => {
    await compoundView(corpValue?.toString());
  });

  const compoundContent = ui.div(createTableFromObject(column, retrievedCompound, ['properties']));
  acc.addPane('Properties', () => createPropertiesTable(column, properties));

  return DG.Widget.fromRoot(ui.divV([corpLink, compoundContent, acc.root]));
}
