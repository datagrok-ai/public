import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {category, test, expect, delay} from "@datagrok-libraries/utils/src/test";
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

  test('substructure-filter-manual', async () => {
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

  test('substructure-filter-panel', async () => {
    let df = DG.DataFrame.fromColumns([grok.data.demo.molecules(1000).columns[0]]);
    let view = grok.shell.addTableView(df);
    view.filters();
    await delay(100);
    (document.getElementsByClassName(
      'panel-titlebar disable-selection panel-titlebar-tabhost')[0].childNodes[0] as HTMLElement).click();
    await delay(100);
    let menuItem = document.getElementsByClassName(
      'd4-menu-item-container d4-vert-menu d4-menu-popup')[0].childNodes[3];
    menuItem.dispatchEvent(new MouseEvent('mouseenter')); await delay(100);
    menuItem.dispatchEvent(new MouseEvent('mousemove')); await delay(100);
    let submenuItem = menuItem.childNodes[0].childNodes[0].childNodes[1];
    submenuItem.dispatchEvent(new MouseEvent('mouseenter')); await delay(100);
    submenuItem.dispatchEvent(new MouseEvent('mousedown')); await delay(100);
    // https://developer.chrome.com/docs/devtools/console/utilities/#getEventListeners-function
    let dialogContents = document.getElementsByClassName('d4-dialog-contents')[0];
    let allButton = dialogContents.childNodes[2].childNodes[1].childNodes[0];
    (allButton as HTMLElement).click();
    let dialogFooter = document.getElementsByClassName('d4-dialog-footer')[0];
    let okButton = dialogFooter.childNodes[0].childNodes[0];
    (okButton as HTMLElement).click(); await delay(100);
    let smilesInput = document.getElementsByClassName('grok-sketcher-input ui-div')[0].childNodes[0];
    (smilesInput as HTMLInputElement).value = 'c1ccccc1'; await delay(100);
    smilesInput.dispatchEvent(new Event('focus')); await delay(100);
    // Only this combination of parameters worked:
    // https://tutorial.eyehunts.com/js/how-to-press-enter-key-programmatically-in-javascript-example-code/
    smilesInput.dispatchEvent(new KeyboardEvent('keydown', {
      altKey: false, bubbles: true, cancelable: true,
      charCode: 0, code: "Enter", composed: true, ctrlKey: false,
      detail: 0, isComposing: false, key: "Enter", keyCode: 13, location: 0,
      metaKey: false, repeat: false, shiftKey: false}));
    await delay(100);
    expect(df.filter.trueCount, 700);
  });
});
