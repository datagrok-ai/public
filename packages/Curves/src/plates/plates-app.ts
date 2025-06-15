import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {createPlatesView} from "./views/plates-create-view";

import {
  initPlates,
  plateTemplates,
  plateUniquePropertyValues,
  queryPlates
} from "./plates-crud";
import { searchWellsView } from './views/plates-search-view';
import { searchPlatesView } from './views/plates-search-view';
import { propertySchemaView as templateView } from './views/plates-schema-view';
import { createTemplatesView } from './views/plates-templates-view';

export function platesAppView(): DG.View {
  const dummy = DG.DataFrame.create(5);
  const view = DG.TableView.create(dummy);

  queryPlates({plateMatchers: [], wellMatchers: []}).then((df: DG.DataFrame) => {
    df.col('barcode')!.semType = 'Plate Barcode';
    view.dataFrame = df;
    view.grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  })

  return view;
}


export async function initPlatesAppTree(treeNode: DG.TreeViewGroup): Promise<void> {
  await initPlates();

  const createPlatesNode = treeNode.item('Create');
  const searchPlatesNode = treeNode.item('Search plates');
  const searchWellsNode = treeNode.item('Search wells');
  createPlatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(createPlatesView()));
  searchPlatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(searchPlatesView()));
  searchWellsNode.onSelected.subscribe(async (_) => grok.shell.addPreview(searchWellsView()));

  const templatesNode = treeNode.group('Templates');
  templatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(createTemplatesView()));
  
  const plateTemplatesNode = templatesNode.group('Plates');
  const wellTemplatesNode = templatesNode.group('Wells');
  
  for (const template of plateTemplates) {
    const plateTemplateNode = plateTemplatesNode.item(template.name);
    plateTemplateNode.onSelected.subscribe(async (_) => grok.shell.addPreview(templateView(template)));
  }
  
  // const queriesNode = treeNode.group('Queries');
  // const nameCol = plateUniquePropertyValues.col('name')!;
  // const valueCol = plateUniquePropertyValues.col('value_string')!;
  // for (const propertyName of nameCol.categories) {
  //   const propertyNode = queriesNode.group(propertyName);
  //   for (const i of plateUniquePropertyValues.rows.where(i => nameCol.get(i) == propertyName))
  //     propertyNode.item(valueCol.get(i), () => grok.shell.info('foo!'));
  // }
}



