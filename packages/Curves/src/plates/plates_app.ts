import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {debounceTime} from "rxjs/operators";
import {NumericMatcher} from "./numeric_matcher";

import {
  getPlateProperty,
  getPlateUniquePropertyValues,
  initPlates,
  plateProperties, PlateQuery,
  plateUniquePropertyValues,
  queryPlates
} from "./plates_crud";

export function platesAppView(): DG.View {
  const dummy = DG.DataFrame.create(5);
  const view = DG.TableView.create(dummy);

  queryPlates({platePropertyValues: {}, platePropertyConditions: []}).then((df: DG.DataFrame) => {
    df.col('barcode')!.semType = 'Plate Barcode';
    view.dataFrame = df;
    view.grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  })

  return view;
}


export async function initPlatesAppTree(treeNode: DG.TreeViewGroup): Promise<void> {
  await initPlates();

  const searchPlatesNode = treeNode.item('Search plates');
  const searchWellsNode = treeNode.item('Search wells');
  searchPlatesNode.onSelected.subscribe(async (_) => grok.shell.addPreview(await searchView()));
  const queriesNode = treeNode.group('Queries');
  const nameCol = plateUniquePropertyValues.col('name')!;
  const valueCol = plateUniquePropertyValues.col('value_string')!;
  for (const propertyName of nameCol.categories) {
    const propertyNode = queriesNode.group(propertyName);
    for (const i of plateUniquePropertyValues.rows.where(i => nameCol.get(i) == propertyName))
      propertyNode.item(valueCol.get(i), () => grok.shell.info('foo!'));
  }
}


async function searchView(): Promise<DG.View> {
  const dummy = DG.DataFrame.create(5, 'Search plates');
  const view = DG.TableView.create(dummy);

  const inputs: DG.InputBase[] = [];
  for (const prop of plateProperties) {
    if (prop.value_type == DG.TYPE.STRING) {
      const items = getPlateUniquePropertyValues(prop);
      const input: DG.InputBase = ui.input.multiChoice(prop.name, { items: items})
      inputs.push(input);
    }
    else if (prop.value_type == DG.TYPE.FLOAT) {
      const input: DG.InputBase = ui.input.string(prop.name);
      inputs.push(input);
    }
  }

  const refresh = () => {
    const query: PlateQuery = { platePropertyValues: {}, platePropertyConditions: []};

    for (const input of inputs.filter(input => input.inputType == DG.InputType.MultiChoice && input.value.length > 0))
      query.platePropertyValues[getPlateProperty(input.caption).id] = input.value;

    for (const input of inputs.filter(input => input.inputType == DG.InputType.Text && input.value !== ''))
      query.platePropertyConditions.push({
        propertyId: getPlateProperty(input.caption).id,
        matcher: NumericMatcher.parse(input.value)!
      });

    queryPlates(query).then((df) => {
      view.dataFrame = df;
      view.grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
        .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
    });
  }

  const form = DG.InputForm.forInputs(inputs);
  form.onInputChanged.pipe(debounceTime(500)).subscribe((_) => refresh());

  refresh();

  //@ts-ignore
  view.dockManager.dock(form.root, DG.DOCK_TYPE.TOP);

  return view;
}