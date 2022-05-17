import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {take} from 'rxjs/operators';
import $ from 'cash-dom';

import * as constants from './const';
import * as alationApi from './alation-api';
import * as types from './types';
import { getBaseURL } from './package';

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

  if (tokenMap.userId === '' || tokenMap.refreshToken === '' || tokenMap.apiToken === '') {
    await updateTokensDialog(tokenMap.refreshToken, userId);
    tokenMap = await getAllTokensFromStorage();
    userId = parseInt(tokenMap.userId === '' ? '0' : tokenMap.userId);
  }

  const refreshResponse = await alationApi.testToken(
    constants.REFRESH_TOKEN_KEY, tokenMap.refreshToken, userId) as types.refreshTokenResponse;
  if (refreshResponse.token_status !== constants.TOKEN_STATUS.ACTIVE) {
    grok.shell.error(`Refresh token status is ${refreshResponse.token_status ?? 'expired'}`);
    await updateTokensDialog(tokenMap.refreshToken, userId);
    tokenMap = await getAllTokensFromStorage();
    userId = parseInt(tokenMap.userId === '' ? '0' : tokenMap.userId);
  }

  const apiResponse = await alationApi.testToken(constants.API_TOKEN_KEY, tokenMap.apiToken, userId);
  if (apiResponse.token_status !== constants.TOKEN_STATUS.ACTIVE) {
    grok.shell.warning(`API access token status is ${apiResponse.token_status ?? 'expired, regenerating...'}`);
    await updateUserStorage(tokenMap.refreshToken, userId);
    tokenMap = await getAllTokensFromStorage();
  }

  return tokenMap;
}

async function updateTokensDialog(refreshToken: string, userId: number) {
  const authLink = ui.link(
    'My Account -> Account Settings -> Authentication', `${await getBaseURL()}${constants.URI_MAP.account_auth})`);
  const refreshTokenHelpText = ui.inlineText(['Refresh token is a long living token the used to manage and create ',
    'API Access Tokens which can be used to interact with the other Alation APIs. It can be found in ',
    authLink, '.']);
  const userIdHelpText = ui.inlineText(['User ID is an ID associated with your account on Alation instance. It can ',
    'be found in My Account -> User Profile. The User ID then will be shown in address bar of your browser.']);
  const dataStorageHelpText = ui.inlineText(['Datagrok stores this data in a secure environment called ',
    ui.link('User Data Storage', 'https://datagrok.ai/help/develop/how-to/user-data-storage'), '.']);
  const helpHost = ui.divV([refreshTokenHelpText, userIdHelpText, dataStorageHelpText], 'alation-help-host');
  $(helpHost).width(500);

  const refreshTokenInput = ui.stringInput('Refresh token', refreshToken);
  const userIdInput = ui.intInput('User ID', userId);
  const dialog = ui.dialog('Update keys')
    .add(ui.divV([helpHost, userIdInput.root, refreshTokenInput.root]))
    .onOK(async () => {
      const refreshTokenInputValue = refreshTokenInput.stringValue;
      const userIdInputValue = userIdInput.value as number;
      const refreshResponse = await alationApi.testToken(
        constants.REFRESH_TOKEN_KEY, refreshTokenInputValue, userIdInputValue) as types.refreshTokenResponse;
      if (refreshResponse.token_status !== constants.TOKEN_STATUS.ACTIVE) {
        grok.shell.error(`Refresh token status is ${refreshResponse.token_status ?? 'expired'}`);
        await updateTokensDialog(refreshTokenInputValue, userIdInputValue);
        return;
      }
      await updateUserStorage(refreshTokenInputValue, userIdInputValue);
    })
    .showModal(false);
  return dialog.onClose.pipe(take(1)).toPromise();
}

async function updateUserStorage(refreshToken: string, userId: number) {
  await UDS.postValue(constants.STORAGE_NAME, constants.REFRESH_TOKEN_KEY, refreshToken, true);
  await UDS.postValue(constants.STORAGE_NAME, constants.USER_ID, `${userId}`, true);
  const apiToken = (await alationApi.createAPIAccessToken(refreshToken, userId)).api_access_token;
  await UDS.postValue(constants.STORAGE_NAME, constants.API_TOKEN_KEY, apiToken, true);
}

async function getAllTokensFromStorage() {
  const refreshToken = await UDS.getValue(constants.STORAGE_NAME, constants.REFRESH_TOKEN_KEY, true);
  const apiToken = await UDS.getValue(constants.STORAGE_NAME, constants.API_TOKEN_KEY, true);
  const userIdStr = await UDS.getValue(constants.STORAGE_NAME, constants.USER_ID, true);
  return {userId: userIdStr, refreshToken: refreshToken, apiToken: apiToken};
}

export async function getApiToken() {
  return await UDS.getValue(constants.STORAGE_NAME, constants.API_TOKEN_KEY, true);
}
