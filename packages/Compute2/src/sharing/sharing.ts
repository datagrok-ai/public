import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-instance';
import {findPublishDialogFunc, publishLiveRun, publishSavedRun} from './artifact-alignment';
import {isWorkspaceSharingAvailable, shareRunToWorkspace} from './workspace-share';

export type SharingMethod = 'none' | 'artifact-alignment' | 'workspaces';

/** What a host offers the share action; the sharing module owns dialogs and server calls. */
export interface ShareTarget {
  /** Live call for methods that can share unsaved state (RFV standalone, TreeWizard step). */
  liveCall?: () => DG.FuncCall;
  /** Persisted run/meta-call id, if the current state is saved. */
  savedCallId?: () => string | null;
  /** Save-then-share; resolves the saved id, null when cancelled or failed. */
  saveRun?: () => Promise<string | null>;
  defaultName?: () => string | undefined;
}

export interface ShareAction {
  tooltip: string;
  run: (target: ShareTarget) => Promise<void>;
}

export function getSharingMethod(): SharingMethod {
  const value = _package.settings?.['sharingMethod'];
  return value === 'artifact-alignment' || value === 'workspaces' ? value : 'none';
}

// The share action for the configured method, or null when sharing is off or unavailable
export function getShareAction(): ShareAction | null {
  const method = getSharingMethod();
  if (method === 'artifact-alignment') {
    const func = findPublishDialogFunc();
    if (func == null)
      return null;
    return {
      tooltip: 'Publish to program',
      run: async (target) => {
        const defaultName = target.defaultName?.();
        if (target.liveCall != null) {
          await publishLiveRun(func, target.liveCall(), defaultName);
          return;
        }
        const savedId = target.savedCallId?.();
        if (savedId != null)
          await publishSavedRun(func, savedId, defaultName);
        else
          grok.shell.warning('Publishing works on a saved run — save the workflow first');
      },
    };
  }
  if (method === 'workspaces') {
    if (!isWorkspaceSharingAvailable())
      return null;
    return {tooltip: 'Share to workspace', run: (target) => shareRunToWorkspace(target)};
  }
  return null;
}
