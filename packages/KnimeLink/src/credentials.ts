import {_package} from './package';

let cachedId = '';
let cachedSecret = '';

export async function getCredentials(): Promise<{id: string; secret: string}> {
  if (cachedId !== '')
    return {id: cachedId, secret: cachedSecret};

  const credentials = await _package.getCredentials();
  if (!credentials)
    throw new Error('KNIME credentials are not set in package credentials');

  const id = credentials.parameters['appPasswordId'];
  const secret = credentials.parameters['appPasswordValue'];
  if (!id || !secret)
    throw new Error('KNIME credentials incomplete: expected "appPasswordId" and "appPasswordValue"');

  cachedId = id;
  cachedSecret = secret;
  return {id: cachedId, secret: cachedSecret};
}

export function getBasicAuthHeader(id: string, secret: string): string {
  return 'Basic ' + btoa(`${id}:${secret}`);
}
