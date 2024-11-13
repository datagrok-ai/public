/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {CampaignGroupingType, CampaignJsonName, ComputeQueryMolColName, HDCampaignsGroupingLSKey, i18n} from './consts';
import {AppName, CampaignsType, HitDesignCampaign, HitTriageCampaign, TriagePermissions} from './types';
import {_package} from '../package';

export const toFormatedDateString = (d: Date): string => {
  return `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;
};

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
  for (const folder of campaignFolders) {
    try {
      const campaignJson: CampaignsType[T] = JSON.parse(await _package.files
        .readAsText(`${appName}/campaigns/${folder.name}/${CampaignJsonName}`));
      if (campaignJson.authorUserId && campaignJson.permissions &&
        !await checkViewPermissions(campaignJson.authorUserId, campaignJson.permissions)
      )
        continue;
      if (campaignJson.authorUserId && !campaignJson.authorUserFriendlyName) {
        const user = await grok.dapi.users.find(campaignJson.authorUserId);
        if (user)
          campaignJson.authorUserFriendlyName = user.friendlyName;
      }
      campaignNamesMap[campaignJson.name] = campaignJson;
    } catch (e) {
      continue;
    }
  }
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
    } else if (userGroupId && group.members.length > 0) {
      if (group.members.some((memberGroup) => memberGroup.id === userGroupId))
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

export function processGroupingTable<T extends HitDesignCampaign | HitTriageCampaign = HitDesignCampaign>(table: HTMLTableElement, groupedCampaigns: {[key: string]: T[]}, numCols = 8) {
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
