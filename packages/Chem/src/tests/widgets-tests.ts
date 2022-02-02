import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {category, test, testExpectFinish} from "@datagrok-libraries/utils/src/test";
import {drugLikenessWidget} from '../widgets/drug-likeness';
import {identifiersWidget} from '../widgets/identifiers';
import {molfileWidget} from '../widgets/molfile';
import {propertiesWidget} from '../widgets/properties';
import {structuralAlertsWidget} from '../widgets/structural-alerts';
import {structure2dWidget} from '../widgets/structure2d';
import {structure3dWidget} from '../widgets/structure3d';
import {toxicityWidget} from '../widgets/toxicity';

category('Chem: Widgets', () => {
  const molStr = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';

  test('drug-likeness', async () => {
    drugLikenessWidget(molStr);
  });

  test('identifiers', async () => {
    identifiersWidget(molStr);
  });

  test('molfile', async () => {
    molfileWidget(molStr);
  });

  test('properties', async () => {
    propertiesWidget(molStr);
  });

  test('structural-alerts', async () => {
    structuralAlertsWidget(molStr);
  });

  test('structure-2d', async () => {
    await grok.functions.call('structure2d', {smiles: molStr});
  });

  test('structure-3d', async () => {
    await grok.functions.call('structure3d', {smiles: molStr});
  });

  test('toxicity', async () => {
    toxicityWidget(molStr);
  });

  testExpectFinish('manual-substructure-filter', async () => {
    let df = grok.data.demo.molecules(1000);
    await grok.data.detectSemanticTypes(df);
    // previously: let filter = await grok.functions.call("Chem:substructureFilter");
    //@ts-ignore
    let filter = chem.substructureFilter();
    filter.attach(df);
    grok.shell.addTableView(df);
    let colChoice = ui.columnInput('Column', filter.dataFrame, filter.column, (col: DG.Column) => {
      filter.column = col;
      filter.dataFrame.filter.setAll(true, false);
      filter.dataFrame.rows.requestFilter();
    });
    ui.dialog({title: 'Chem Filter'})
      .add(colChoice)
      .add(filter.root)
      .show();
  });
});
