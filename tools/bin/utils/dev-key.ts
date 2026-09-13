import * as keypair from './keypair';

const fetch = require('node-fetch');

/** `Authorization` header carrying the developer key. */
export function devKeyHeaders(key: string, headers: any = {}): any {
  return {...headers, 'Authorization': `Dev ${key}`};
}

/** One login per API root per process: every publish step reuses the token. */
const tokenCache = new Map<string, Promise<string>>();

/**
 * Bearer token obtained by signing a nonce with the server's registered keypair,
 * or `null` when no keypair is configured for [apiRoot] and the caller should
 * fall back to the developer key.
 */
export async function keypairToken(apiRoot: string, devKey?: string): Promise<string | null> {
  const privateKey = keypair.keypairFor(apiRoot, devKey);
  if (!privateKey)
    return null;
  if (!tokenCache.has(apiRoot))
    tokenCache.set(apiRoot, keypair.keyLogin(apiRoot, privateKey));
  return await tokenCache.get(apiRoot)!;
}

/**
 * Calls [url] as the configured user. With a keypair (`grok login`) that is a
 * session token; otherwise the developer key rides in the `Authorization`
 * header, falling back to [legacyUrl] - which carries the key as a path segment -
 * for servers that predate the header form. Those answer 404 (no such route) or
 * 401 (the route-less path is not on their anonymous allow-list); a server that
 * knows the header form answers neither.
 *
 * [apiRoot] enables the keypair path; without it this stays dev-key only.
 */
export async function devKeyFetch(url: string, legacyUrl: string, key: string,
  init: any = {}, apiRoot?: string): Promise<any> {
  if (apiRoot != null) {
    const token = await keypairToken(apiRoot, key);
    if (token)
      return await fetch(url, {...init, headers: {...init.headers, 'Authorization': token}});
  }
  const response = await fetch(url, {...init, headers: devKeyHeaders(key, init.headers)});
  if (response.status !== 404 && response.status !== 401)
    return response;
  return await fetch(legacyUrl, init);
}
