import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';
import { UserSummary, Organization, ArchiveRecord } from './types';

export type ProjectOwner = Organization | UserSummary;

// Project interface
export interface Project {
  archiveRecord?: ArchiveRecord | null;
  id: string;
  name: string;
  owner?: ProjectOwner;
}

export interface ProjectsPaginatedList {
  projects: Project[];
  nextToken?: string;
}

export interface ProjectsQueryParams {
  nextToken?: string;
  pageSize?: number;
  sort?: string;
  archiveReason?: string;
  ids?: string;
  name?: string;
}

// Mock data
export const mockProjects: ProjectsPaginatedList = {
  projects: Array.from({length: 10}, (_, i) => ({
    id: `prj_${String(i+1).padStart(3, '0')}`,
    name: `Project ${i+1}`,
    archiveRecord: i % 3 === 0 ? { reason: 'Retired' } : null,
    owner: i % 2 === 0
      ? { handle: `org${i+1}`, id: `org_${i+1}`, name: `Organization ${i+1}` }
      : { handle: `user${i+1}`, id: `ent_${i+1}`, name: `User ${i+1}` },
  })),
  nextToken: undefined,
};

/**
 * Query projects (GET /projects)
 * @param params Query parameters
 * @returns DataFrame of projects
 */
export async function queryProjects(params: ProjectsQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/projects${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.projects ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockProjects.projects);
  return df;
} 