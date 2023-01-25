import * as grok from 'datagrok-api/grok';

import * as constants from './const';
import * as alationApi from './alation-api';
import * as types from './types';
import {getPackage} from './package';


export function filterDuplicates(objects: types.baseEntity[]): types.baseEntity[] {
  const ids = new Set<number>();
  return objects.filter(({id}) => {
    const isNotKnownId = !ids.has(id);
    if (isNotKnownId)
      ids.add(id);
    return isNotKnownId;
  });
}

export async function retrieveKeys() {
  let tokenMap = await getAllTokensFromStorage();

  const apiResponse = await alationApi.testToken(constants.API_TOKEN_KEY, tokenMap.apiToken, tokenMap.userId);
  if (apiResponse.token_status !== constants.TOKEN_STATUS.ACTIVE) {
    grok.shell.warning('Creating API Acess Token...');

    const createApiTokenResponse = await alationApi.createAPIAccessToken(tokenMap.refreshToken, tokenMap.userId);

    updateApiToken(createApiTokenResponse.api_access_token);
    tokenMap = await getAllTokensFromStorage();
  }

  return tokenMap;
}

async function getAllTokensFromStorage(): Promise<{userId: number, refreshToken: string, apiToken: string}> {
  const creds = await getPackage().getCredentials();
  const credParams = creds.parameters;
  return {
    userId: parseInt(credParams[constants.USER_ID]),
    refreshToken: credParams[constants.REFRESH_TOKEN_KEY],
    apiToken: getApiToken(),
  };
}

function updateApiToken(apiToken: string): void {
  localStorage.setItem(constants.API_TOKEN_KEY, apiToken);
}

export function getApiToken(): string {
  const apiToken = localStorage.getItem(constants.API_TOKEN_KEY);
  return apiToken ?? '';
}
