import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {initPlates, plateProperties} from "./plates_crud";

export function platesAppView(): DG.View {
  const dummy = DG.DataFrame.create(5);
  const view = DG.TableView.create(dummy);

  grok.functions.call('Curves:getPlates').then((df: DG.DataFrame) => {
    df.col('barcode')!.semType = 'Plate Barcode';
    view.dataFrame = df;
    view.grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  })

  return view;
}


export async function initPlatesAppTree(treeNode: DG.TreeViewGroup): Promise<void> {
  await initPlates();

  const df: DG.DataFrame = await grok.functions.call('Curves:getUniquePlatePropertyValues');
  const queriesNode = treeNode.group('Queries');
  const nameCol = df.col('name')!;
  const valueCol = df.col('value_string')!;
  for (const propertyName of nameCol.categories) {
    const propertyNode = queriesNode.group(propertyName);
    for (const i of df.rows.where(i => nameCol.get(i) == propertyName))
      propertyNode.item(valueCol.get(i), () => grok.shell.info('foo!'));
  }
}