import * as grok from 'datagrok-api/grok';
import {PortfolioRow} from '../generated/db';
import {portfolioByCode} from '../portal';
import {withRestrictedUser} from './restricted-user';

/** A plain user to check what someone sees: the stand's named test user when its password is in the
 * test runner's user settings, else a throwaway restricted user put into the same groups.
 *
 * - `viewer`: in EX-SMD's and EX-BIOC's viewers groups (the stand's `ppviewer`);
 * - `outsider`: in EX-BIOD's viewers group only (the stand's `ppoutsider`).
 *
 * The passwords stay out of the repository: an admin stores them once on the stand with
 * `grok.userSettings.add('ProjectPortal tests', 'ppviewer', '<password>')`. */
export type StandUserKind = 'viewer' | 'outsider';

const SETTINGS = 'ProjectPortal tests';
const LOGINS: Record<StandUserKind, string> = {viewer: 'ppviewer', outsider: 'ppoutsider'};
const GROUPS: Record<StandUserKind, string[]> = {viewer: ['EX-SMD', 'EX-BIOC'], outsider: ['EX-BIOD']};

export interface StandUser {
  /** The login the session belongs to. */
  login: string;
  /** The user's id. */
  id: string;
  /** Runs [action] as that user, restoring the admin session after. */
  asUser<T>(action: () => Promise<T>): Promise<T>;
}

/** A session token for [login], or null when the password is not stored or the login fails. The
 * session cookie is the admin's and stays untouched (credentials: 'omit'). */
async function tokenOf(login: string): Promise<string | null> {
  let password: string | null = null;
  try {
    password = grok.userSettings.getValue(SETTINGS, login) ?? null;
  } catch (_) {
    password = null;
  }
  if (!password)
    return null;
  const r = await fetch(`${grok.dapi.root}/users/login`, {method: 'POST', credentials: 'omit',
    headers: {'Content-Type': 'application/json'}, body: JSON.stringify({login, password: btoa(password)})});
  const body = await r.json().catch(() => null);
  return body?.isSuccess && body.token ? body.token : null;
}

/** Runs [body] as a [kind] user. Resolves to null when no such user can be had on this stand. */
export async function withStandUser<T>(kind: StandUserKind, body: (u: StandUser) => Promise<T>): Promise<T | null> {
  const login = LOGINS[kind];
  const token = await tokenOf(login);
  if (token) {
    const user = await grok.dapi.users.filter(`login = "${login}"`).first();
    return body({login, id: user.id, asUser: async <R>(action: () => Promise<R>): Promise<R> => {
      grok.dapi.impersonationToken = token;
      try {
        return await action();
      } finally {
        grok.dapi.impersonationToken = null;
      }
    }});
  }
  // no named user: a throwaway one in the same viewers groups
  return withRestrictedUser(`pp${kind}`, async (u) => {
    const personal = await grok.dapi.groups.find(u.group);
    const groups = [];
    for (const code of GROUPS[kind]) {
      const p: PortfolioRow = await portfolioByCode(code);
      if (p.viewers)
        groups.push(await grok.dapi.groups.find(p.viewers));
    }
    for (const g of groups)
      await grok.dapi.groups.includeTo(personal, g);
    try {
      return await body({login: u.login, id: u.id, asUser: u.asUser});
    } finally {
      for (const g of groups)
        await grok.dapi.groups.excludeFrom(personal, g).catch(() => undefined);
    }
  });
}
