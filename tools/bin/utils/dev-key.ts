const fetch = require('node-fetch');

/** `Authorization` header carrying the developer key. */
export function devKeyHeaders(key: string, headers: any = {}): any {
  return {...headers, 'Authorization': `Dev ${key}`};
}

/**
 * Calls [url] with the developer key in the `Authorization` header, falling back to
 * [legacyUrl] - which carries the key as a path segment - for servers that predate the
 * header form. Those answer 404 (no such route) or 401 (the route-less path is not on
 * their anonymous allow-list); a server that knows the header form answers neither.
 */
export async function devKeyFetch(url: string, legacyUrl: string, key: string, init: any = {}): Promise<any> {
  const response = await fetch(url, {...init, headers: devKeyHeaders(key, init.headers)});
  if (response.status !== 404 && response.status !== 401)
    return response;
  return await fetch(legacyUrl, init);
}
