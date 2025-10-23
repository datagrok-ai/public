/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
// import {createPlatesView} from './views/viewcreate';

import * as crud from './plates-crud';
import {searchWellsView} from './views/plates-search-view';
import {searchPlatesView} from './views/plates-search-view';
import {propertySchemaView as templateView} from './views/plates-schema-view';
import {createTemplatesView} from './views/plates-templates-view';
import {filter} from 'rxjs/operators';
import {createPlatesView} from './views/plates-create-view';
import {AnalysisManager} from '../plate/analyses/analysis-manager';
import {searchAnalysesView} from './views/analyses-search-view';

export function platesAppView(): DG.View {
  const dummy = DG.DataFrame.create(5);
  const view = DG.TableView.create(dummy);
  view.name = 'Plates';


  crud.queryPlates({plateMatchers: [], wellMatchers: [], analysisMatchers: []}).then((df: DG.DataFrame) => {
    df.col('barcode')!.semType = 'Plate Barcode';
    view.dataFrame = df;
    view.grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await plates.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  });

  return view;
}


export async function initPlatesAppTree(treeNode: DG.TreeViewGroup): Promise<void> {
  await crud.initPlates();

  try {
    console.log('Plates: Initializing AnalysisManager...');
    await AnalysisManager.instance.init();
    console.log('Plates: AnalysisManager initialized successfully.');
  } catch (e) {
    console.error('Plates: Failed to initialize AnalysisManager:', e);
  }

  const createPlatesNode = treeNode.item('Create');
  const searchPlatesNode = treeNode.item('Search plates');
  const searchWellsNode = treeNode.item('Search wells');
  const searchAnalysesNode = treeNode.item('Search analyses');

  createPlatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(createPlatesView()));
  searchPlatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(await searchPlatesView()));
  searchWellsNode.onSelected.subscribe(async (_) => grok.shell.addPreview(await searchWellsView()));
  searchAnalysesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(await searchAnalysesView())); // <-- AND THIS ONE


  const templatesNode = treeNode.group('Templates');
  templatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(createTemplatesView()));

  const plateTemplatesNode = templatesNode.group('Plates');
  const wellTemplatesNode = templatesNode.group('Wells');

  for (const template of crud.plateTemplates) {
    const plateTemplateNode = plateTemplatesNode.item(template.name, DG.SemanticValue.fromValueType(template, crud.TYPE.TEMPLATE));
    plateTemplateNode.onSelected.subscribe(async (_) => grok.shell.addPreview(templateView(template)));
  }

  // auto-delete templates from the tree when they are deleted
  crud.events.pipe(filter((event) => event.eventType === 'deleted')).subscribe((event) => {
    treeNode.removeChildrenWhere((node) => node.value instanceof DG.SemanticValue &&
       crud.entityTypes.includes(node.value.semType) &&
       node.value.value.id && node.value.value.id === event.object.id);
  });

  // const queriesNode = treeNode.group('Queries');
  // const nameCol = plateUniquePropertyValues.col('name')!;
  // const valueCol = plateUniquePropertyValues.col('value_string')!;
  // for (const propertyName of nameCol.categories) {
  //   const propertyNode = queriesNode.group(propertyName);
  //   for (const i of plateUniquePropertyValues.rows.where(i => nameCol.get(i) == propertyName))
  //     propertyNode.item(valueCol.get(i), () => grok.shell.info('foo!'));
  // }
}
