import * as grok from 'datagrok-api/grok';
import {GITHUB_BASE_URL} from './constants';

export async function fetchGithubFile(fileName: string): Promise<string> {
  const response = await grok.dapi.fetchProxy(`${GITHUB_BASE_URL}${fileName}`);
  if (!response.ok)
    throw new Error(`Failed to fetch ${fileName}: ${response.statusText}`);
  return await response.text();
}
