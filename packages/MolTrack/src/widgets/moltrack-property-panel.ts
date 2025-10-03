import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import { EntityBaseView } from '../views/registration-entity-base';
import { compoundView } from '../utils/view-utils';

export async function molTrackPropPanel(retrievedCompound: any, initialStructure?: string): Promise<DG.Widget> {
  const {canonical_smiles: smiles, detail, batches = [], properties = []} = retrievedCompound;
  if (detail) {
    const registerButton = ui.bigButton('Register', async () => {
      const registrationView = new EntityBaseView(false);
      if (initialStructure)
        registrationView.initialSmiles = initialStructure;
      await registrationView.buildUIMethod();
      registrationView.show();
    });
    registerButton.classList.add('moltrack-panel-register-button');
    const div = ui.divV([ui.divText('No matches'), registerButton], 'moltrack-panel-register-container');
    return DG.Widget.fromRoot(div);
  }

  const acc = ui.accordion();

  const createTableFromObject = (obj: Record<string, unknown>, excludeKeys: string[] = []) => {
    const map: Record<string, unknown> = {};
    for (const [key, value] of Object.entries(obj))
      if (!excludeKeys.includes(key)) map[key] = value;

    return ui.tableFromMap(map);
  };

  const extractPropValue = (prop: Record<string, any>) =>
    prop.value_uuid ?? prop.value_num ?? prop.value_string ?? prop.value_qualifier ?? prop.value_datetime ?? null;

  const createPropertiesTable = (properties: any[]) => {
    if (!properties || !properties.length) return ui.divText('No properties');
    const map: Record<string, unknown> = {};
    for (const prop of properties) map[prop.name] = extractPropValue(prop);
    return ui.tableFromMap(map);
  };

  const createNestedAccordion = (items: any[], idKey: string) => {
    if (!items || !items.length) return ui.divText('No items');
    const nested = ui.accordion();
    for (const item of items)
      nested.addPane(item[idKey], () => createTableFromObject(item));

    return nested.root;
  };

  const panes = [
    {
      title: 'Batches',
      content: () => createNestedAccordion(batches, 'id'),
    },
    {
      title: 'Properties',
      content: () => createPropertiesTable(properties),
    },
  ];

  // --- Create the corporate_compound_id link first ---
  const corpProp = properties.find((p: any) => p.name === 'corporate_compound_id');
  const value = extractPropValue(corpProp);
  const corpLink = ui.link(value.toString() ?? '', async () => {
    await compoundView(value.toString());
  });

  const compoundContent = ui.div(createTableFromObject(retrievedCompound, ['properties', 'batches']));

  for (const pane of panes) acc.addPane(pane.title, pane.content);

  return DG.Widget.fromRoot(ui.divV([corpLink, compoundContent, acc.root]));
}
