/** Who may run the portfolio-keeping actions: sync, building or parking the workspace, export and
 * advance by data. The functions check the role themselves, not only the buttons, but package code
 * runs in the user's browser: a check here is the app's, and the server's own grants are listed
 * beside each action, so it is clear which actions a console call could still run.
 *
 * Pure: the caller says whether the user is one of the portfolio's approvers (admins are). */

export type KeeperAction = 'sync' | 'workspace' | 'park' | 'export' | 'auto_advance';

export interface KeeperRule {
  /** The verb phrase in a refusal: "Only EX-SMD's approvers sync it". */
  what: string;
  /** What the server enforces on its own, whatever the app checks. */
  server: string;
}

export const KEEPER_ACTIONS: Record<KeeperAction, KeeperRule> = {
  sync: {what: 'sync it',
    server: 'program rows and their work units take writes from their author and admins only; sync_log from the ' +
      'approvers groups. Compound links and gate-log rows are secured through the program'},
  workspace: {what: 'build its workspace',
    server: 'none beyond the space grants: the contributors group holds Edit on the root, so a console call by a ' +
      'contributor could still create folders'},
  park: {what: 'park or unpark its folders',
    server: 'parked_folder takes writes from the approvers groups; moving the folder needs Edit on the root, which ' +
      'contributors hold'},
  export: {what: 'export it',
    server: 'none: the portfolio table is readable by all users, so the file could be assembled from the console'},
  auto_advance: {what: 'advance its programs by data',
    server: 'the gate log is secured through the program row: author and admins only'},
};

/** Board buttons that keep a portfolio, and the action each runs. */
export const KEEPER_BUTTONS: Record<string, KeeperAction> = {sync: 'sync', workspace: 'workspace', export: 'export'};

/** Why [action] is refused in portfolio [code], or null when it is allowed. */
export function keeperRefusal(action: KeeperAction, code: string, isApprover: boolean): string | null {
  return isApprover ? null : `Only ${code}'s approvers ${KEEPER_ACTIONS[action].what}`;
}

/** Whether a board button shows: a keeper's button only for the portfolio's approvers, and the
 * advanced ones only in Advanced mode. */
export function showsBoardButton(key: string, opts: {simple: boolean, advanced: string[], isApprover: boolean}): boolean {
  if (KEEPER_BUTTONS[key] && !opts.isApprover)
    return false;
  return !(opts.simple && opts.advanced.includes(key));
}
