import * as grok from 'datagrok-api/grok';

import * as types from './types';
import * as constants from './const';
import * as utils from './utils';
import {getBaseURL} from './package';

export async function testToken(
  tokenType: types.alationTokenType, token: string, userId: number): Promise<types.tokenResponse> {
  const url = `${await getBaseURL()}integration/v1/${constants.URI_MAP[tokenType]}/`;
  const data: types.alationDataType = {user_id: userId};
  data[tokenType] = token;
  const params = {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    },
    body: JSON.stringify(data),
  };

  const response = await grok.dapi.fetchProxy(url, params);
  return response.json();
}

export async function createAPIAccessToken(
  refreshToken: string, userId: number): Promise<types.createApiTokenResponse> {
  const url = `${await getBaseURL()}integration/v1/${constants.URI_MAP.create_api_token}/`;
  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify({user_id: userId, refresh_token: refreshToken}),
  };

  const response = await grok.dapi.fetchProxy(url, params);
  return await response.json();
}

export async function getDataSources(): Promise<types.dataSource[]> {
  const url = `${await getBaseURL()}integration/v1/${constants.URI_MAP['datasource']}/`;
  return getResponseFor(url) as Promise<types.dataSource[]>;
}

export async function getSchemas(dataSourceId: number): Promise<types.schema[]> {
  const url = `${await getBaseURL()}integration/v2/${constants.URI_MAP.schema}/?ds_id=${dataSourceId}&limit=100&skip=0`;
  return getResponseFor(url) as Promise<types.schema[]>;
}

export async function getTables(schemaId: number): Promise<types.table[]> {
  const url = `${await getBaseURL()}integration/v2/${constants.URI_MAP.table}/?schema_id=${schemaId}&limit=100&skip=0`;
  return getResponseFor(url) as Promise<types.table[]>;
}

export async function getColumns(tableId: number): Promise<types.column[]> {
  const url = `${await getBaseURL()}integration/v2/${constants.URI_MAP.column}/?table_id=${tableId}&limit=100&skip=0`;
  return getResponseFor(url) as Promise<types.column[]>;
}

async function getResponseFor(url: string): Promise<types.baseEntity[]>;
async function getResponseFor(url: string, params: RequestInit): Promise<types.baseEntity[]>;
async function getResponseFor(url: string, returnJson: boolean): Promise<types.baseEntity[] | string>;
async function getResponseFor(
  url: string, params: RequestInit, returnJson: boolean): Promise<types.baseEntity[] | string>;
async function getResponseFor(
  url: string, paramsOrReturnJson?: RequestInit | boolean, returnJson?: boolean,
): Promise<types.baseEntity[] | string> {
  let params: RequestInit;
  if (typeof paramsOrReturnJson === 'boolean' || typeof paramsOrReturnJson === 'undefined') {
    returnJson = paramsOrReturnJson;
    params = {
      method: 'GET',
      headers: {
        'Accept': 'application/json',
        'TOKEN': utils.getApiToken(),
      },
    };
  }

  params ??= paramsOrReturnJson as RequestInit;
  returnJson ??= true;

  const response = await grok.dapi.fetchProxy(url, params);
  return returnJson ? response.json() : response.text();
}

export async function getQueries(
  dataSourceId: number, saved: boolean = true, published: boolean = true): Promise<types.query[]> {
  const urlParams = `datasource_id=${dataSourceId}&saved=${saved}&published=${published}`;
  const url = `${await getBaseURL()}integration/v1/${constants.URI_MAP.query}/?${urlParams}`;
  const params = {
    method: 'GET',
    headers: {
      'Accept': 'application/json',
      'TOKEN': utils.getApiToken(),
    },
  };

  const response = await grok.dapi.fetchProxy(url, params);
  return response.json();
}

export async function testQuery(execResultId: number): Promise<string> {
  const url =
    `${await getBaseURL()}integration/v1/${constants.URI_MAP.result}/${execResultId}/${constants.FORMAT.csv}/`;
  return getResponseFor(url, false) as Promise<string>;
}

export async function getDataSourceById(dataSourceId: number): Promise<types.dataSource> {
  const url = `${await getBaseURL()}integration/v1/${constants.URI_MAP.datasource}/${dataSourceId}`;
  const params = {
    method: 'GET',
    headers: {
      'Accept': 'application/json',
      'TOKEN': utils.getApiToken(),
    },
  };

  const response = await grok.dapi.fetchProxy(url, params);
  return response.json();
}

export async function getTableObject(table_id: string): Promise<types.table> {
  const url = `${await getBaseURL()}integration/v2/table/?id=${table_id}&limit=1`;
  const params = {
    method: 'GET',
    headers: {
      'Accept': 'application/json',
      'TOKEN': utils.getApiToken(),
    },
  };
  const response = await grok.dapi.fetchProxy(url, params);
  const json = await response.json();
  return json[0];
}

export async function getQueryObject(query_id: string): Promise<types.query> {
  const url = `${await getBaseURL()}integration/v1/query/${query_id}/`;
  const params = {
    method: 'GET',
    headers: {
      'Accept': 'application/json',
      'TOKEN': utils.getApiToken(),
    },
  };
  const response = await grok.dapi.fetchProxy(url, params);
  return response.json();
}
