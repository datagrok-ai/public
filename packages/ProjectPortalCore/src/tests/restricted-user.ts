import * as grok from 'datagrok-api/grok';
import {expect} from '@datagrok-libraries/test/src/test';

// Copied from packages/ApiTests/src/dapi/domain-lifecycle.ts (withRestrictedUser): a throwaway
// self-signed-up user, impersonated through grok.dapi.impersonationToken and blocked afterwards.

/** A freshly signed-up user with no privileges, and the way to act as them. */
export interface RestrictedUser {
  login: string;
  id: string;
  /** The user's personal group — what a grant is addressed to. */
  group: string;
  /** Runs [action] under THAT user's session, restoring the admin one after. */
  asUser<T>(action: () => Promise<T>): Promise<T>;
}

/** Runs [body] with a throwaway restricted user (login prefix [prefix]).
 *
 * The session cookie is HttpOnly, so script can neither read nor overwrite it:
 * impersonation goes through `grok.dapi.impersonationToken`, which puts the
 * probe's token in the Authorization header the server reads ahead of the
 * cookie. Clearing it restores the admin session — the cookie never moved. The
 * finally clears it again and BLOCKS the user: users are not API-deletable, and
 * blocking revokes their sessions including the signup one. The signup itself
 * runs inside that try and the finally resolves the user by login when the body
 * never did, so neither a failed setup nor a tokenless signup can leave an
 * unblocked user behind.
 *
 * Resolves to null, reason logged, where the harness cannot support it:
 * unavailable self-signup (SSO-only or email-confirm setups). */
export async function withRestrictedUser<T>(prefix: string,
  body: (user: RestrictedUser) => Promise<T>): Promise<T | null> {
  const restoreAdmin = () => grok.dapi.impersonationToken = null;
  const stamp = `${Date.now()}${Math.floor(Math.random() * 1e4)}`;
  const login = `${prefix}${stamp}`;
  let user: any = null;
  try {
    // credentials: 'omit' is load-bearing. Signup mints a session and the server
    // answers with a Set-Cookie for the NEW user; same-origin fetch defaults to
    // sending and STORING cookies, so the admin's HttpOnly `auth` cookie would be
    // overwritten and — being HttpOnly — unrestorable. The whole run would
    // silently continue as the probe user, which the finally then blocks.
    const signup = await (await fetch(`${grok.dapi.root}/users/signup`, {
      method: 'POST', credentials: 'omit', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({login: login, email: `${login}@test.datagrok.ai`,
        password: btoa(`Pw-${stamp}`), firstName: 'ApiTests', lastName: 'Probe'}),
    })).json();
    // No token, but the signup may still have created the user — the finally
    // below looks it up by login and blocks it.
    if (!signup?.token) {
      console.log(`skipped: self-signup unavailable (${signup?.reason ?? 'no token'})`);
      return null;
    }
    user = await grok.dapi.users.filter(`login = "${login}"`).include('group').first();
    expect(user != null, true, 'signed-up user not found');
    return await body({
      login: login,
      id: user.id,
      group: user.group.id,
      asUser: async <R>(action: () => Promise<R>): Promise<R> => {
        grok.dapi.impersonationToken = signup.token;
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
    if (user?.id == null)
      console.log(`no test user ${login} to block`);
    else {
      const blocked = await fetch(`${grok.dapi.root}/users/block`, {
        method: 'POST', credentials: 'include', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({'#type': 'User', 'id': user.id, 'login': login}),
      });
      expect(blocked.ok, true, `test user ${login} not blocked: ${blocked.status}`);
    }
  }
}
