/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {CampaignGroupingType, CampaignJsonName, CampaignTableColumns, ComputeQueryMolColName, HDCampaignsGroupingLSKey, HDCampaignTableColumnsLSKey, HTFunctionOrderingLSKey, i18n} from './consts';
import {AppName, CampaignsType, HitDesignCampaign, HitTriageCampaign, TriagePermissions} from './types';
import {_package} from '../package';

export const toFormatedDateString = (d: Date): string => {
  return `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;
};

/**
 * Modifies the current URL by updating or adding a query parameter with the specified key and value.
 * Updates the browser's history state without reloading the page.
 *
 * @param key - The query parameter key to be added or updated in the URL.
 * @param value - The value to be assigned to the specified query parameter key.
 *
 * @remarks
 * This function uses the `history.replaceState` method to update the URL and browser history state.
 * It ensures that the base URL remains unchanged and appends the query parameter.
 * If `history.replaceState` is not supported, the function will not modify the URL.
 *
 * @example
 * ```typescript
 * modifyUrl('page', '2');
 * // Updates the URL to something like: http://example.com/?page=2
 * ```
 */
export function modifyUrl(key: string, value: string) {
  const title = document.title;
  const url = window.location.href.split('?')[0] + '?' + key + '=' + value;
  if (history.replaceState) {
    const obj = {
      Title: title,
      Url: url,
    };
    history.replaceState(obj, obj.Title, obj.Url);
  }
}

export function hideComponents(styles: CSSStyleDeclaration[]) {
  styles.forEach((style) => {
    style.display = 'none';
  });
}

export function showComponents(styles: CSSStyleDeclaration[]) {
  styles.forEach((style) => {
    style.display = 'flex';
  });
}

export function addBreadCrumbsToRibbons(view: DG.ViewBase, firstPath: string, secondPath: string, handler: Function):
  {breadcrumbs: DG.Breadcrumbs, sub: Subscription} {
  const breadcrumbs = ui.breadcrumbs([firstPath, secondPath]);
  //breadcrumbs.root.id = 'hit-triage-breadcrumbs';
  const sub = breadcrumbs.onPathClick.subscribe((path) => {
    if (path.length === 1) {
      sub.unsubscribe();
      popRibbonPannels(view);
      handler();
    }
  });
  const viewNameRoot = grok.shell.v.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
  if (viewNameRoot) {
    viewNameRoot.textContent = '';
    viewNameRoot.appendChild(breadcrumbs.root);
  }
  /*
  let ribbons = view.getRibbonPanels();
  if (!ribbons?.length)
    ribbons = [];
  ribbons.push([breadcrumbs.root]);

  if (document.getElementById('hit-triage-breadcrumbs') === null) {
    view.setRibbonPanels(ribbons);
    breadcrumbs.root.parentElement?.classList.add('d4-ribbon-name');
    breadcrumbs.root.parentElement?.classList.remove('d4-ribbon-item');
  }*/
  return {breadcrumbs, sub};
}

export function popRibbonPannels(view: DG.ViewBase) {
  const ribbons = view.getRibbonPanels();
  if (!ribbons?.length)
    return;
  ribbons.pop();
  view.setRibbonPanels(ribbons);
};

export function checkRibbonsHaveSubmit(ribbons: HTMLElement[][]): boolean {
  return ribbons.reduce((prev, cur) => {
    return prev || cur.reduce((p, c) =>
      p || c.classList.contains('hit-design-submit-button') ||
        Array.from(c.children).some((child) =>child.classList.contains('hit-design-submit-button')), false);
  }, false);
}

/** Executes {@link func} while showing the "running" indicator on {@link root}.
     * Handles and logs exceptions.
     * @param {HTMLElement} root
     * @param {Function} func
     *  */
export async function runAsync<T>(root: HTMLElement, func: () => Promise<T>) {
  ui.setUpdateIndicator(root, true);
  try {
    return await func();
  } catch (e: any) {
    grok.log.error(e);
  } finally {
    ui.setUpdateIndicator(root, false);
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
};

export async function loadCampaigns<T extends AppName>(
  appName: T, deletedCampaigns: string[],
): Promise<{[name: string]: CampaignsType[T]}> {
  const campaignFolders = (await _package.files.list(`${appName}/campaigns`))
    .filter((f) => deletedCampaigns.indexOf(f.name) === -1);
  const campaignNamesMap: {[name: string]: CampaignsType[T]} = {};
  const campaignPromises = campaignFolders.map(async (folder) => {
    try {
      const campaignJson: CampaignsType[T] = JSON.parse(await _package.files
        .readAsText(`${appName}/campaigns/${folder.name}/${CampaignJsonName}`));
      if (campaignJson.authorUserId && campaignJson.permissions &&
        !await checkViewPermissions(campaignJson.authorUserId, campaignJson.permissions)
      )
        return null;
      if (campaignJson.authorUserId && !campaignJson.authorUserFriendlyName) {
        const user = await grok.dapi.users.find(campaignJson.authorUserId);
        if (user)
          campaignJson.authorUserFriendlyName = user.friendlyName;
      }
      campaignNamesMap[campaignJson.name] = campaignJson;
    } catch (e) {
      console.error(e);
    }
    return null;
  });
  await Promise.all(campaignPromises);
  return campaignNamesMap;
}

async function checkPermissions(authorId: string, groupIdList: string[]): Promise<boolean> {
  const userId = DG.User.current().id;
  const userGroupId = DG.User.current().group?.id;
  if (authorId === userId)
    return true;
  for (const groupId of groupIdList) {
    const group = await grok.dapi.groups.find(groupId);
    if (!group)
      continue;
    if (group.personal) {
      const user = await grok.dapi.groups.getUser(group);
      if (user?.id === userId)
        return true;
    } else if (userGroupId) {
      if ((group.members?.length ?? 0) > 0 && group.members.some((memberGroup) => memberGroup.id === userGroupId))
        return true;
      if ((group.adminMembers?.length ?? 0) > 0 && group.adminMembers.some((memberGroup) => memberGroup.id === userGroupId))
        return true;
    }
  }
  return false;
}

export async function checkEditPermissions(authorId: string, permissions: TriagePermissions): Promise<boolean> {
  return checkPermissions(authorId, permissions.edit);
}

export async function checkViewPermissions(authorId: string, permissions: TriagePermissions): Promise<boolean> {
  return checkPermissions(authorId, Array.from(new Set([...permissions.view, ...permissions.edit])));
}

export function getLocalStorageValue<T = string>(key: string): T | null {
  return localStorage.getItem(key) as unknown as T;
}

export function setLocalStorageValue(key: string, value: string) {
  localStorage.setItem(key, value as unknown as string);
}

export function getSavedCampaignsGrouping(): CampaignGroupingType {
  return getLocalStorageValue<CampaignGroupingType>(HDCampaignsGroupingLSKey) ?? CampaignGroupingType.None;
}

export function setSavedCampaignsGrouping(value: CampaignGroupingType) {
  setLocalStorageValue(HDCampaignsGroupingLSKey, value);
}

export function getSavedCampaignTableColumns(): CampaignTableColumns[] {
  const columnsString = getLocalStorageValue<string>(HDCampaignTableColumnsLSKey);
  const defaultCols: CampaignTableColumns[] = ['Code', 'Created', 'Author', 'Last Modified by', 'Molecules', 'Status'];
  if (!columnsString)
    return defaultCols;
  try {
    const columnsP = JSON.parse(columnsString);
    if (Array.isArray(columnsP) && columnsP.every((c) => typeof c === 'string'))
      return columnsP as CampaignTableColumns[];
    return defaultCols;
  } catch (e) {
    console.error('error parsing campaign table columns', e);
    return defaultCols;
  }
}

export function setSavedCampaignTableColumns(columns: CampaignTableColumns[]) {
  setLocalStorageValue(HDCampaignTableColumnsLSKey, JSON.stringify(columns));
}


export const getGroupingKey = <T extends HitDesignCampaign | HitTriageCampaign = HitDesignCampaign>(grouping: CampaignGroupingType, campaign: T): string => {
  switch (grouping) {
  case CampaignGroupingType.Template:
    return campaign.template?.key ?? campaign.templateName;
  case CampaignGroupingType.Status:
    return campaign.status;
  case CampaignGroupingType.Author:
    return campaign.authorUserFriendlyName ?? i18n.noInformation;
  case CampaignGroupingType.LastModifiedUser:
    return campaign.lastModifiedUserName ?? i18n.noInformation;
  default:
    return '';
  }
};

export function getGroupedCampaigns<T extends HitDesignCampaign | HitTriageCampaign = HitDesignCampaign>(campaigns: T[], grouping: CampaignGroupingType):
  {[key: string]: T[]} {
  if (grouping === CampaignGroupingType.None)
    return {'': campaigns};
  const groupedCampaigns: {[key: string]: T[]} = {};
  for (const campaign of campaigns) {
    const key = getGroupingKey(grouping, campaign);
    if (!groupedCampaigns[key])
      groupedCampaigns[key] = [];
    groupedCampaigns[key].push(campaign);
  }
  return groupedCampaigns;
}

export function processGroupingTable<T extends HitDesignCampaign | HitTriageCampaign = HitDesignCampaign>(table: HTMLTableElement, groupedCampaigns: {[key: string]: T[]}, numCols = 9) {
  table.classList.add('hit-design-groupped-campaigns-table');
  const keys = Object.keys(groupedCampaigns);
  if (keys.length < 2)
    return table;
  const body = table.getElementsByTagName('tbody')[0];
  const rows = Array.from(table.getElementsByTagName('tr')).filter((row) => !row.classList.contains('header'));
  let curRow = 0;

  const setState = (expanded: boolean, start: number, end: number) => {
    for (let i = start; i < end; i++)
      rows[i]?.style && (rows[i].style.display = expanded ? 'table-row' : 'none');
  };


  for (const key of keys) {
    const row = rows[curRow];
    if (!row)
      break;
    const l = groupedCampaigns[key].length;
    const startRow = curRow;
    const endRow = curRow + l;
    const acc = ui.accordion(`Hit-Design-campaigns-group-${key}`);
    const pane = acc.addPane(key, () => {return ui.div();}, undefined, undefined, false);
    pane.root.style.marginLeft = '-24px';
    pane.root.onclick = () => {
      setState(pane.expanded, startRow, endRow);
    };
    setState(pane.expanded, startRow, endRow);
    const newRow = body.insertRow(0);
    const newCell = newRow.insertCell(0);
    newCell.appendChild(pane.root);
    newCell.colSpan = numCols;
    body.insertBefore(newRow, row);
    curRow += l;
  }
}

export type EditableFieldOptions = {
  onChange?: (val?: string | null) => void,
  onOk?: (val?: string | null) => Promise<void>,
  validator?: (val?: string | null) => Promise<string | null>,
  onCancel?: () => void,
  afterEditTextContent?: (val: string) => string,
  beforeEditTextContent?: (val: string) => string,
  nullable?: boolean,
  tooltip?: string,
}

export function editableTableField(field: HTMLElement, options?: EditableFieldOptions) {
  const editIcon = ui.icons.edit(() => {
    let tooltipMsg = options?.tooltip ?? field.textContent ?? '';
    const beforeEditTextContent = options?.beforeEditTextContent ?? ((val) => val);
    const afterEditTextContent = options?.afterEditTextContent ?? ((val) => val);
    const input = ui.input.string('smth', {value: beforeEditTextContent(field.textContent ?? ''), onValueChanged: async () => {
      const vr = await internalValidator();
      setTimeout(() => {
        if (vr) {
          input.input.classList.add('d4-invalid');
          tooltipMsg = vr;
          return;
        } else {
          input.input.classList.remove('d4-invalid');
          tooltipMsg = options?.tooltip ?? field.textContent ?? '';
        }
      }, 100);
      options?.onChange?.(input.value);
    }});
    ui.tooltip.bind(input.input, () => tooltipMsg);
    const labelElement = input.root.getElementsByTagName('label').item(0);
    if (labelElement)
      labelElement.remove();
    input.root.style.width = '100%';
    const internalValidator = async () => {
      const initialRes = !!input.value?.trim() ? null : !!options?.nullable ? null :'Field cannot be empty';
      return initialRes ?? (options?.validator ? await options.validator(input.value) : null);
    };

    const saveButton = ui.button('Save', async () => {
      const vr = await internalValidator();
      if (vr) {
        grok.shell.error(vr);
        return;
      }
      await options?.onOk?.(input.value);
      field.textContent = afterEditTextContent(input.value) ?? '';
      ui.empty(container);
      container.appendChild(field);
      container.appendChild(editIcon);
    });
    const cancelButton = ui.button('Cancel', () => {
      ui.empty(container);
      container.appendChild(field);
      container.appendChild(editIcon);
      options?.onCancel?.();
    });

    ui.empty(container);
    input.addOptions(cancelButton);
    input.addOptions(saveButton);
    container.appendChild(input.root);
    input.input.focus();
  }, options?.tooltip ?? 'Edit');
  const container = ui.divH([field, editIcon], {style: {display: 'flex', alignItems: 'center'}});
  editIcon.style.marginLeft = '5px';
  return container;
}

export async function checkFileExists(path: string) {
  if (!path || path.trim() === '') {
    grok.shell.error('Path can not be empty');
    return false;
  }
  const exists = await grok.dapi.files.exists(path);
  if (!exists) {
    grok.shell.error('Given folder does not exist');
    return false;
  }
  return true;
}

export function timeoutOneTimeEventListener(element: HTMLElement, eventName: string, callback: () => void, timeout = 4000) {
  function listener() {
    callback();
    element.removeEventListener(eventName, listener);
  };
  element.addEventListener(eventName, listener);
  setTimeout(() => {
    element.removeEventListener(eventName, listener);
  }, timeout);
}

// //name: Demo Design with reinvent
// //input: int numberOfMolecules
// //tags: HitDesignerFunction
// //output: dataframe result
// export async function demoDesignWithReinvent(numberOfMolecules: number): Promise<DG.DataFrame> {
//   const df = grok.data.demo.molecules(numberOfMolecules);
//   df.name = 'Reinvent Design';
//   return df;
// }

// #region Function ordering

export type FunctionOrdering = {
  order: string[],
  hidden: string[]
}

export function getSavedFunctionOrdering(): FunctionOrdering {
  const orderingString = getLocalStorageValue<string>(HTFunctionOrderingLSKey) ?? '{}';
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

export function setSavedFunctionOrdering(ordering: FunctionOrdering) {
  setLocalStorageValue(HTFunctionOrderingLSKey, JSON.stringify(ordering));
}

export function getReorderedFunctionTabArgs(args: {[key: string]: HTMLElement}) {
  const ordering = getSavedFunctionOrdering();
  const orderedArgs: {[key: string]: HTMLElement} = {};
  for (const key of ordering.order) {
    if (args[key])
      orderedArgs[key] = args[key];
  }
  for (const key in args) {
    if (!orderedArgs[key] && !ordering.hidden.includes(key))
      orderedArgs[key] = args[key];
  }
  return orderedArgs;
}

export function getReorderingInput(functions: string[], onOk: (ordering: FunctionOrdering) => void) {
  const order = getSavedFunctionOrdering();
  const dataFrame = DG.DataFrame.fromColumns(functions.map((f) => DG.Column.fromStrings(f, [f])));
  dataFrame.columns.setOrder(order.order ?? []);
  const columnsEditor = ui.input.columns('reorder', {table: dataFrame, value: dataFrame.columns.toList().filter((c) => !order.hidden.includes(c.name))});
  columnsEditor.onChanged.subscribe(() => {
    try {
      const chosenColumns = columnsEditor.value.map((c) => c.name);
      const hiddenColumns = functions.filter((f) => !chosenColumns.includes(f));
      const newOrdering = {
        order: columnsEditor.value.map((c) => c.name),
        hidden: hiddenColumns,
      };
      dataFrame.columns.setOrder(newOrdering.order);
      setSavedFunctionOrdering(newOrdering);
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

export function getFuncPackageNameSafe(func: DG.Func): string | undefined {
  try {
    return func.package?.name;
  } catch (e) {
    return undefined;
  }
}
