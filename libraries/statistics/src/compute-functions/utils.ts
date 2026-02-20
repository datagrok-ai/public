import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ComputeQueryMolColName} from './consts';

// #region Function ordering

export type FunctionOrdering = {
  order: string[],
  hidden: string[]
}

export function getSavedFunctionOrdering(storageKey: string): FunctionOrdering {
  const orderingString = localStorage.getItem(storageKey) ?? '{}';
  try {
    const orderingP = JSON.parse(orderingString);
    const ordering: FunctionOrdering = {
      order: orderingP?.order ?? [],
      hidden: orderingP?.hidden ?? [],
    };
    return ordering;
  } catch (e) {
    console.error('error parsing function ordering', e);
  }
  return {order: [], hidden: []};
}

export function setSavedFunctionOrdering(ordering: FunctionOrdering, storageKey: string) {
  localStorage.setItem(storageKey, JSON.stringify(ordering));
}

export function getReorderingInput(
  functions: string[], onOk: (ordering: FunctionOrdering) => void, storageKey: string,
) {
  const order = getSavedFunctionOrdering(storageKey);
  const dataFrame = DG.DataFrame.fromColumns(functions.map((f) => DG.Column.fromStrings(f, [f])));
  dataFrame.columns.setOrder(order.order ?? []);
  const columnsEditor = ui.input.columns('reorder', {table: dataFrame, value: dataFrame.columns.toList().filter((c) => !order.hidden.includes(c.name))});
  columnsEditor.onChanged.subscribe(() => {
    try {
      const chosenColumns = columnsEditor.value.map((c: DG.Column) => c.name);
      const hiddenColumns = functions.filter((f) => !chosenColumns.includes(f));
      const newOrdering = {
        order: columnsEditor.value.map((c: DG.Column) => c.name),
        hidden: hiddenColumns,
      };
      dataFrame.columns.setOrder(newOrdering.order);
      setSavedFunctionOrdering(newOrdering, storageKey);
      onOk(newOrdering);
    } catch (e) {
      console.error(e);
    }
  });

  const children = Array.from(columnsEditor.root.children) as HTMLElement[];
  setTimeout(() => {
    children.forEach((child) => {
      child.style.maxWidth = '0px';
      child.style.overflow = 'hidden';
      child.style.padding = '0px';
      child.style.paddingRight = '0px';
      child.style.visibility = 'hidden';
      if (child instanceof HTMLLabelElement)
        child.style.display = 'none';
    });
  }, 200);
  columnsEditor.root.style.justifyContent = 'end';
  columnsEditor.root.style.width = '40px';
  columnsEditor.root.style.height = '0px';
  columnsEditor.root.style.overflow = 'visible';
  columnsEditor.root.style.padding = '0px';


  const editIcon = ui.icons.edit(() => {
    children.forEach((child) => {
      child.click();
    });
  }, 'Order or hide functions');

  columnsEditor.addOptions(editIcon);
  (Array.from(columnsEditor.root.children) as HTMLElement[]).forEach((child) => {
    child.style.borderBottom = 'unset';
  });
  editIcon.style.fontSize = '16px';
  return columnsEditor.root;
}

// #endregion

export function getFuncPackageNameSafe(func: DG.Func): string | undefined {
  try {
    return func.package?.name;
  } catch (e) {
    return undefined;
  }
}

export async function joinQueryResults(df: DG.DataFrame, molColName: string, qRes: DG.DataFrame) {
  if (qRes.rowCount === 0)
    return;
  const molCol = df.col(molColName);
  if (!molCol)
    throw new Error('There is no molecule column in dataframe');
  let resOriginalCol = qRes.col(ComputeQueryMolColName);
  if (!resOriginalCol) {
    await qRes.meta.detectSemanticTypes();
    resOriginalCol = qRes.columns.bySemType(DG.SEMTYPE.MOLECULE);
  }
  if (!resOriginalCol)
    throw new Error('There is no original molecule column in query result dataframe');
  const newColNames = qRes.columns.names().filter((name) => name !== ComputeQueryMolColName);
  df.join(qRes, [molColName], [resOriginalCol.name], undefined, newColNames, undefined, true);
}
