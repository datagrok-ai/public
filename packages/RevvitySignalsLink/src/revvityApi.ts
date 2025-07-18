// Revvity Signals REST API - Response Interfaces

import { getApiKey, getApiUrl } from "./credentials-utils";
import * as grok from 'datagrok-api/grok';
import { SignalsSearchParams, SignalsSearchQuery } from "./signalsSearchQuery";

// Top-level response interface
export interface RevvityApiResponse<T = any, I = any> {
  links?: RevvityLinks;
  data?: RevvityData<T>;
  included?: RevvityIncluded<I>[];
  errors?: RevvityApiError[];
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


async function request<T>(
  method: string,
  path: string,
  body?: any,
  text?: boolean,
): Promise<RevvityApiResponse<T>> {
  const headers: any = {
    'x-api-key': await getApiKey(),
    'Accept': 'application/json',
  };
  if (method === 'POST')
    headers['Content-Type'] = 'application/json';

  const url = `${await getApiUrl()}${path}`;
  const response = await grok.dapi.fetchProxy(url, {
    method,
    headers,
    body: body ? JSON.stringify(body) : undefined,
  });

  const data = text ? await response.bytes() : await response.json();

  if (!response.ok) {
    throw new Error(data.error ?? `HTTP error!: ${response.status}`, { cause: response.status });
  }
  return { data };
}


/** Search entities */
export async function searchEntities(body: SignalsSearchQuery, queryParams: SignalsSearchParams): Promise<RevvityApiResponse> {
  const paramsStr = encodeURI(paramsStringFromObj(queryParams));
  return request('POST', `/entities/search?${paramsStr}`, body);
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
