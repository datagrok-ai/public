/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {CampaignJsonName, HTFunctionOrderingLSKey} from './consts';
import {AppName, CampaignsType, TriagePermissions} from './types';
import {_package} from '../package';

// Re-export shared compute-function utilities from the statistics library
export {joinQueryResults, getFuncPackageNameSafe} from '@datagrok-libraries/statistics/src/compute-functions/utils';
import {
  getSavedFunctionOrdering as _getSavedFunctionOrdering,
  setSavedFunctionOrdering as _setSavedFunctionOrdering,
  getReorderingInput as _getReorderingInput,
} from '@datagrok-libraries/statistics/src/compute-functions/utils';
import type {FunctionOrdering} from '@datagrok-libraries/statistics/src/compute-functions/utils';
export type {FunctionOrdering};
import type {TemplateCompute} from '@datagrok-libraries/statistics/src/compute-functions/types';

// Wrappers that bind the HitTriage-specific localStorage key
export function getSavedFunctionOrdering(): FunctionOrdering {
  return _getSavedFunctionOrdering(HTFunctionOrderingLSKey);
}

export function setSavedFunctionOrdering(ordering: FunctionOrdering) {
  _setSavedFunctionOrdering(ordering, HTFunctionOrderingLSKey);
}

export function getReorderingInput(functions: string[], onOk: (ordering: FunctionOrdering) => void) {
  return _getReorderingInput(functions, onOk, HTFunctionOrderingLSKey);
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

export async function loadCampaigns<T extends AppName>(
  appName: T, deletedCampaigns: string[],
): Promise<{[name: string]: CampaignsType[T]}> {
  const campaignFolders = (await _package.files.list(`${appName}/campaigns`))
    .filter((f) => deletedCampaigns.indexOf(f.name) === -1);
  const campaignNamesMap: {[name: string]: CampaignsType[T]} = {};
  const membership: MembershipCache = new Map();
  const campaignPromises = campaignFolders.map(async (folder) => {
    try {
      const campaignJson: CampaignsType[T] = JSON.parse(await _package.files
        .readAsText(`${appName}/campaigns/${folder.name}/${CampaignJsonName}`));
      if (campaignJson.authorUserId && campaignJson.permissions &&
        !await checkViewPermissions(campaignJson.authorUserId, campaignJson.permissions, membership)
      )
        return null;
      if (campaignJson.authorUserId && !campaignJson.authorUserFriendlyName) {
        const user = await grok.dapi.users.find(campaignJson.authorUserId);
        if (user)
          campaignJson.authorUserFriendlyName = user.friendlyName;
      }
      campaignNamesMap[campaignJson.name] = campaignJson;
    } catch (e) {
      console.warn(`Skipped campaign ${appName}/${folder.name}:`, e);
    }
    return null;
  });
  await Promise.all(campaignPromises);
  return campaignNamesMap;
}

/** The groups the current user was looked up in, each once: the campaigns of one list share it. */
export type MembershipCache = Map<string, Promise<boolean>>;

async function isMember(groupId: string): Promise<boolean> {
  const user = DG.User.current();
  const group = await grok.dapi.groups.find(groupId);
  if (!group)
    return false;
  if (group.personal)
    return (await grok.dapi.groups.getUser(group))?.id === user.id;
  const own = user.group?.id;
  return !!own && [...group.members ?? [], ...group.adminMembers ?? []].some((member) => member?.id === own);
}

async function checkPermissions(authorId: string, groupIdList: string[], cache?: MembershipCache): Promise<boolean> {
  try {
    // the server does not list the members of All users, so the lookup below would never find the user there
    if (authorId === DG.User.current().id || groupIdList.includes(DG.Group.defaultGroupsIds['All users']))
      return true;
    for (const groupId of groupIdList) {
      let member = cache?.get(groupId);
      if (member === undefined) {
        member = isMember(groupId);
        cache?.set(groupId, member);
      }
      if (await member)
        return true;
    }
  } catch (e) {
    _package.logger.error(e);
  }
  return false;
}

export async function checkEditPermissions(authorId: string, permissions: TriagePermissions,
  cache?: MembershipCache): Promise<boolean> {
  return checkPermissions(authorId, permissions.edit, cache);
}

export async function checkViewPermissions(authorId: string, permissions: TriagePermissions,
  cache?: MembershipCache): Promise<boolean> {
  return checkPermissions(authorId, Array.from(new Set([...permissions.view, ...permissions.edit])), cache);
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

export async function loadTemplate<T extends {compute: TemplateCompute}>(path: string): Promise<T> {
  return JSON.parse(await _package.files.readAsText(path)) as T;
}
