/**
 * A throwaway restricted (non-admin) user for the live permission tests — a copy of ApiTests'
 * `withRestrictedUser` (public/packages/ApiTests/src/dapi/domain-lifecycle.ts), which another
 * package cannot import.
 *
 * The session cookie is HttpOnly, so impersonation goes through `grok.dapi.impersonationToken`,
 * which puts the probe's token in the Authorization header the server reads ahead of the cookie.
 * The finally clears it and BLOCKS the user (users are not API-deletable). Resolves to null, with
 * the reason logged, where self-signup is unavailable (SSO-only or email-confirm setups).
 */

import * as grok from 'datagrok-api/grok';

import {isLoopbackServer, SHARED_STAND_SKIP_REASON} from './shared-stand-guard';

export interface RestrictedUser {
  readonly login: string;
  readonly id: string;
  /** The user's personal group — what a grant is addressed to. */
  readonly group: string;
  /** Runs `action` under THAT user's session, restoring the admin one after. */
  asUser<T>(action: () => Promise<T>): Promise<T>;
}

export async function withRestrictedUser<T>(prefix: string,
  body: (user: RestrictedUser) => Promise<T>): Promise<T | null> {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const dapi = grok.dapi as any;
  // A probe user can only be blocked, never deleted: refuse to make one on a shared server.
  if (!isLoopbackServer(grok.dapi.root))
    throw new Error(`Probe user refused on a shared server (${grok.dapi.root}): ${SHARED_STAND_SKIP_REASON}`);
  const restoreAdmin = (): void => {dapi.impersonationToken = null;};
  const stamp = `${Date.now()}${Math.floor(Math.random() * 1e4)}`;
  const login = `${prefix}${stamp}`;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  let user: any = null;
  try {
    // credentials: 'omit' is load-bearing: a stored Set-Cookie of the NEW user would overwrite the
    // admin's HttpOnly session for the rest of the run.
    const signup = await (await fetch(`${grok.dapi.root}/users/signup`, {
      method: 'POST', credentials: 'omit', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({login, email: `${login}@test.datagrok.ai`, password: btoa(`Pw-${stamp}`),
        firstName: 'NcaTests', lastName: 'Probe'}),
    })).json();
    if (!signup?.token) {
      console.log(`skipped: self-signup unavailable (${signup?.reason ?? 'no token'})`);
      return null;
    }
    user = await grok.dapi.users.filter(`login = "${login}"`).include('group').first();
    if (user == null) throw new Error('signed-up user not found');
    return await body({
      login, id: user.id, group: user.group.id,
      asUser: async <R>(action: () => Promise<R>): Promise<R> => {
        dapi.impersonationToken = signup.token;
        try {
          return await action();
        } finally {
          restoreAdmin();
        }
      },
    });
  } finally {
    restoreAdmin();
    if (user?.id == null) {
      try {
        user = await grok.dapi.users.filter(`login = "${login}"`).first();
      } catch (x) {
        console.error(`test user ${login} could not be looked up — block it manually: ${x}`);
      }
    }
    if (user?.id != null) {
      // The public endpoint (what `grok s users block` calls): `/users/block` answers 200 and
      // blocks nothing — every run used to leave its probe user active (found 2026-11-02).
      await fetch(`${grok.dapi.root}/public/v1/users/block`, {
        method: 'POST', credentials: 'include', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({'#type': 'User', 'id': user.id, 'login': login}),
      });
    }
  }
}
