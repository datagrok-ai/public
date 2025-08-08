// Revvity Signals REST API - Response Interfaces

import { getApiKey, getApiUrl } from "./credentials-utils";
import * as grok from 'datagrok-api/grok';
import { SignalsSearchParams, SignalsSearchQuery } from "./signals-search-query";

// Top-level response interface
export interface RevvityApiResponse<T = any, I = any> {
  meta?: RevvityMeta;
  links?: RevvityLinks;
  data?: RevvityData<T> | RevvityData<T>[];
  included?: RevvityIncluded<I>[];
  errors?: RevvityApiError[];
}

//Meta data included into response
export interface RevvityMeta {
  "took-ms": number,
  "query-timed-out": boolean,
  "query-reached-limit": boolean,
  count: number,
  total: number
}

// Links object (can be extended as needed)
export interface RevvityLinks {
  self: string;
  [key: string]: any;
}

// Data object
export interface RevvityData<T = any> {
  type: string;
  id: string;
  links?: RevvityLinks;
  attributes?: T;
  relationships?: RevvityRelationships;
}

// Included object (same structure as data, but can be more generic)
export type RevvityIncluded<I = any> = RevvityData<I>;

// Relationships object (can be extended as needed)
export interface RevvityRelationships {
  [key: string]: {
    links?: RevvityLinks;
    data?: { type: string; id: string } | { type: string; id: string }[];
    meta?: any;
  };
}

// Error object
export interface RevvityApiError {
  status: string;
  code: string;
  title: string;
  detail: string;
}

export interface RevvityUserLicense {
  id: string;
  name: string;
  active: boolean;
}

export interface RevvityUser {
  lastLoginAt?: string;
  createdAt?: string;
  licenses?: RevvityUserLicense[];
  userId?: string;
  userName?: string;
  email?: string;
  firstName?: string;
  lastName?: string;
  country?: string;
  organization?: string;
  isEnabled?: boolean;
}


async function request<T>(
  method: string,
  path: string,
  body?: any,
  text?: boolean,
): Promise<RevvityApiResponse<T>> {
  const headers: any = {
    'x-api-key': await getApiKey(),
    'Accept': 'application/vnd.api+json',
  };
  if (method === 'POST')
    headers['Content-Type'] = 'application/json';

  const url = `${await getApiUrl()}${path}`;
  const response = await grok.dapi.fetchProxy(url, {
    method,
    headers,
    body: body ? JSON.stringify(body) : undefined,
  });

  const res: RevvityApiResponse<T> = text ? await response.bytes() : await response.json();

  if (!response.ok || res.errors) {
    if (res.errors)
      throw res.errors;
    throw new Error(`HTTP error!: ${response.status}`, { cause: response.status });
  }
  return res;
}

export async function queryUsers(): Promise<RevvityApiResponse> {
  return request('GET', `/users`);
}

export async function queryEntities(body: SignalsSearchQuery, queryParams?: SignalsSearchParams): Promise<RevvityApiResponse> {
  const paramsStr = queryParams ? `?${encodeURI(paramsStringFromObj(queryParams))}` : '';
  return request('POST', `/entities/search${paramsStr}`, body);
}

export async function queryEntityById(id: string, queryParams?: any): Promise<RevvityApiResponse> {
  const paramsStr = queryParams ? `?${encodeURI(paramsStringFromObj(queryParams))}` : '';
  return request('GET', `/entities/${id}${paramsStr}`);
}

export async function queryMaterialById(id: string, queryParams?: any): Promise<RevvityApiResponse> {
  const paramsStr = queryParams ? `?${encodeURI(paramsStringFromObj(queryParams))}` : '';
  return request('GET', `/materials/${id}${paramsStr}`);
}


export async function queryStructureById(id: string): Promise<string> {

  const url = `${await getApiUrl()}/materials/${id}/drawing?format=smiles`;
  const response = await grok.dapi.fetchProxy(url, {
    method: 'GET',
    headers: {
      'x-api-key': await getApiKey(),
    }
  });

  const res: string = await response.text();

  if (!response.ok) {
    throw new Error(`HTTP error!: ${response.status}`, { cause: response.status });
  }
  return res;
}

export function paramsStringFromObj(params: any): string {
  let paramsStr = '';
  const paramNames = Object.keys(params);
  for (let i = 0; i < paramNames.length; i++) {
    let paramVal = params[paramNames[i]];
    if (paramVal) {
      if (Array.isArray(paramVal))
        paramVal = paramVal.join(',')
      paramsStr += paramsStr === '' ? `${paramNames[i]}=${paramVal}` : `&${paramNames[i]}=${paramVal}`;
    }
  }
  return paramsStr;
}
