// Revvity Signals REST API - Response Interfaces

import { getApiKey, getApiUrl } from "./credentials-utils";
import * as grok from 'datagrok-api/grok';
import { SignalsSearchParams, SignalsSearchQuery } from "./signals-search-query";
import { MAX_RETURN_ROWS } from "./constants";
import { delay } from '@datagrok-libraries/test/src/test';

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
  customHeaders?: any,
  text?: boolean,
): Promise<RevvityApiResponse<T>> {
  const maxRetries = 9;
  let lastError: Error | string | undefined;
  let lastResponseStatus: number | undefined;
  
  for (let attempt = 0; attempt <= maxRetries; attempt++) {
    try {
      const headers: any = customHeaders ?? {
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

      lastResponseStatus = response.status;

      // Handle 429 (Too Many Requests) - retry after delay
      if (response.status === 429 && attempt < maxRetries) {
       // console.log(`${attempt} attempt for ${path} failed`);
        await delay(1000);
        continue; // Retry the request
      }

      const res: RevvityApiResponse<T> = text ? await response.text() : await response.json();

      if (!response.ok || res.errors) {
        if (res.errors) {
          //check for 403 error to further check in getUsers method
          if (res.errors.length === 1 && res.errors[0].status === '403')
            throw '403';
          throw res.errors.map((error) => error.detail).join(';');
        }
        throw new Error(`HTTP error!: ${response.status}`, { cause: response.status });
      }
      
    //  console.log(`${attempt} attempt for ${path} succeeded`);
      return res;
    } catch (error) {
      lastError = error as Error | string;
      
      // Only retry if we got a 429 error and have retries left
      if (lastResponseStatus === 429 && attempt < maxRetries) {
        await delay(1000);
        continue;
      }
      
      // For any other error or if we've exhausted retries, throw immediately
      throw lastError;
    }
  }
  
  throw lastError || new Error('Request failed after retries');
}

export async function queryUsers(): Promise<RevvityApiResponse> {
  return request('GET', `/users`);
}

export async function search(body: SignalsSearchQuery, queryParams?: SignalsSearchParams): Promise<RevvityApiResponse> {
  const params = queryParams ?? {};
  if (!params['page[limit]'])
    params['page[limit]'] = MAX_RETURN_ROWS;
  const paramsStr = `?${encodeURI(paramsStringFromObj(params))}`;
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

export async function queryLibraries(): Promise<RevvityApiResponse> {
  return request('GET', `/materials/libraries`);
}

export async function queryTerms(body: SignalsSearchQuery): Promise<RevvityApiResponse> {
  return request('POST', `/entities/search/terms`, body);
}

export async function queryTags(body: SignalsSearchQuery): Promise<RevvityApiResponse> {
  return request('POST', `/entities/search/tags`, body);
}


export async function queryStructureById(id: string): Promise<RevvityApiResponse> {
  return request('GET', `/materials/${id}/drawing?format=smiles`, undefined,
    {
      'x-api-key': await getApiKey(),
    },
    true);
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
