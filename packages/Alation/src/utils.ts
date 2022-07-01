import * as grok from 'datagrok-api/grok';

import * as constants from './const';
import * as alationApi from './alation-api';
import * as types from './types';

const UDS = grok.dapi.userDataStorage;

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
  let userId = parseInt(tokenMap.userId === '' ? '0' : tokenMap.userId);

  if (tokenMap.refreshToken == '') {
    grok.shell.warning('Creating refresh token...');

    const credentials = await getCredentialsFromUDS();
    if (credentials.username == '' || credentials.password == '')
      throw new Error('Service account credentails are not set');

    const createRefreshTokenResponse =
      await alationApi.createRefreshToken(credentials.username, credentials.password, constants.REFRESH_TOKEN_KEY);
    userId = createRefreshTokenResponse.user_id;

    await updateUserStorage(createRefreshTokenResponse.refresh_token, userId, tokenMap.apiToken);
    tokenMap = await getAllTokensFromStorage();
  }

  const refreshResponse = await alationApi.testToken(constants.REFRESH_TOKEN_KEY, tokenMap.refreshToken, userId);
  if (refreshResponse.token_status !== constants.TOKEN_STATUS.ACTIVE) {
    grok.shell.warning('Regenerating refresh token...');

    const regenerateRefreshTokenResponse = await alationApi.regenerateRefreshToken(tokenMap.refreshToken, userId);

    await updateUserStorage(regenerateRefreshTokenResponse.refresh_token, userId, tokenMap.apiToken);
    tokenMap = await getAllTokensFromStorage();
  }

  const apiResponse = await alationApi.testToken(constants.API_TOKEN_KEY, tokenMap.apiToken, userId);
  if (apiResponse.token_status !== constants.TOKEN_STATUS.ACTIVE) {
    grok.shell.warning('Creating API Acess Token...');

    const createApiTokenResponse = await alationApi.createAPIAccessToken(tokenMap.refreshToken, userId);

    await updateUserStorage(tokenMap.refreshToken, userId, createApiTokenResponse.api_access_token);
    tokenMap = await getAllTokensFromStorage();
  }

  return tokenMap;
}

async function updateUserStorage(
  refreshToken: string, userId: number, apiToken: string, currentUser: boolean = false): Promise<void> {
  await UDS.postValue(constants.STORAGE_NAME, constants.REFRESH_TOKEN_KEY, refreshToken, currentUser);
  await UDS.postValue(constants.STORAGE_NAME, constants.USER_ID, `${userId}`, currentUser);
  // const apiToken = (await alationApi.createAPIAccessToken(refreshToken, userId)).api_access_token;
  await UDS.postValue(constants.STORAGE_NAME, constants.API_TOKEN_KEY, apiToken, currentUser);
}

async function getAllTokensFromStorage(currentUser: boolean = false) {
  const refreshToken = await UDS.getValue(constants.STORAGE_NAME, constants.REFRESH_TOKEN_KEY, currentUser);
  const apiToken = await UDS.getValue(constants.STORAGE_NAME, constants.API_TOKEN_KEY, currentUser);
  const userIdStr = await UDS.getValue(constants.STORAGE_NAME, constants.USER_ID, currentUser);

  return {userId: userIdStr, refreshToken: refreshToken, apiToken: apiToken};
}

async function getCredentialsFromUDS(currentUser: boolean = false) {
  const username = await UDS.getValue(constants.STORAGE_NAME, constants.SERVICE_USERNAME, currentUser);
  const password = await UDS.getValue(constants.STORAGE_NAME, constants.SERVICE_PASSWORD, currentUser);

  return {username: username, password: password};
}

export async function getApiToken(currentUser: boolean = false) {
  return UDS.getValue(constants.STORAGE_NAME, constants.API_TOKEN_KEY, currentUser);
}
