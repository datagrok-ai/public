import {randomUUID} from 'node:crypto';

export {runWithContext, request, toGrokPath, listFiles, downloadFile, uploadFile} from './shared-api-client.js';
import {request} from './shared-api-client.js';

export function listFunctions(filter?: string): Promise<unknown[]> {
  const qs = filter ? `?text=${encodeURIComponent(filter)}` : '';
  return request('GET', `/public/v1/functions${qs}`);
}

export function getFunction(id: string): Promise<unknown> {
  const name = id.replace(/:/g, '.');
  return request('GET', `/public/v1/functions/${name}`);
}

export function callFunction(name: string, params?: Record<string, unknown>): Promise<unknown> {
  return request('POST', `/public/v1/functions/${name.replace(/:/g, '.')}/call`, params ?? {});
}

export function saveFunction(data: Record<string, unknown>): Promise<unknown> {
  return request('POST', '/public/v1/functions', data);
}

export function getCurrentUser(): Promise<unknown> {
  return request('GET', '/public/v1/users/current');
}

// --- Projects ---

export function listProjects(filter?: string): Promise<unknown[]> {
  const qs = filter ? `?text=${encodeURIComponent(filter)}` : '';
  return request('GET', `/projects${qs}`);
}

export function getProject(id: string): Promise<unknown> {
  return request('GET', `/projects/${id}`);
}

export function createProject(name: string, description?: string): Promise<unknown> {
  return request('POST', '/projects', {id: randomUUID(), name, description: description ?? ''});
}

export function deleteProject(id: string): Promise<unknown> {
  return request('DELETE', `/projects/${id}`);
}

export function searchProject(name: string): Promise<unknown[]> {
  return request('GET', `/projects/search?name=${encodeURIComponent(name)}`);
}

export function listRecentProjects(): Promise<unknown[]> {
  return request('GET', '/projects/recent');
}

// --- Spaces ---

export function listSpaces(filter?: string): Promise<unknown[]> {
  const qs = filter ? `?text=${encodeURIComponent(filter)}` : '';
  return request('GET', `/spaces${qs}`);
}

export function getSpace(id: string): Promise<unknown> {
  return request('GET', `/spaces/${id}`);
}

export function createRootSpace(name: string): Promise<unknown> {
  return request('POST', '/spaces', {name});
}

export function deleteSpace(id: string): Promise<unknown> {
  return request('DELETE', `/spaces/${id}`);
}

export function createSubspace(spaceId: string, name: string, link?: boolean): Promise<unknown> {
  const qs = link ? '?link=true' : '';
  return request('POST', `/spaces/${spaceId}/subspaces${qs}`, {name});
}

export function listSpaceChildren(
  spaceId: string, types?: string, includeLinked?: boolean,
): Promise<unknown[]> {
  const params = new URLSearchParams();
  if (types)
    params.set('types', types);
  if (includeLinked !== undefined)
    params.set('includeLinked', String(includeLinked));
  const qs = params.toString() ? `?${params}` : '';
  return request('GET', `/spaces/${spaceId}/children${qs}`);
}

export function addEntityToSpace(spaceId: string, entityId: string, link?: boolean): Promise<unknown> {
  const qs = link ? '?link=true' : '';
  return request('POST', `/spaces/${spaceId}/entities/${entityId}${qs}`);
}

export function removeEntityFromSpace(spaceId: string, entityId: string): Promise<unknown> {
  return request('DELETE', `/spaces/${spaceId}/entities/${entityId}`);
}

export function readSpaceFile(spaceId: string, path: string): Promise<unknown> {
  return request('GET', `/spaces/${spaceId}/files/${path.replace(/^\//, '')}`);
}

export function writeSpaceFile(spaceId: string, path: string, content: string): Promise<unknown> {
  return request('POST', `/spaces/${spaceId}/files/${path.replace(/^\//, '')}`, content,
    {'Content-Type': 'application/octet-stream'});
}

export function deleteSpaceFile(spaceId: string, path: string): Promise<unknown> {
  return request('DELETE', `/spaces/${spaceId}/files/${path.replace(/^\//, '')}`);
}
