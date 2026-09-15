import * as DG from 'datagrok-api/dg';
import {_package} from '../package-instance';
import {isWorkspaceSharingAvailable, shareRunToWorkspace} from './workspace-share';

export type SharingMethod = 'none' | 'workspaces';

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
  return value === 'workspaces' ? value : 'none';
}

// The share action for the configured method, or null when sharing is off or unavailable
export function getShareAction(): ShareAction | null {
  const method = getSharingMethod();
  if (method === 'workspaces') {
    if (!isWorkspaceSharingAvailable())
      return null;
    return {tooltip: 'Share to workspace', run: (target) => shareRunToWorkspace(target)};
  }
  return null;
}
