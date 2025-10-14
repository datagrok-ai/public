import { GITHUB_BASE_URL, Scope } from './constants';

export async function fetchSchema(fileName: string): Promise<any> {
  const url = `${GITHUB_BASE_URL}${fileName}`;

  const response = await fetch(url);
  if (!response.ok)
    throw new Error(`Failed to fetch schema for ${fileName}: ${response.statusText}`);
  return await response.json();
}

export async function fetchCsv(scope: Scope): Promise<string> {
  const url = `${GITHUB_BASE_URL}${scope}.csv`;
  const res = await fetch(url);
  if (!res.ok)
    throw new Error(`Failed to fetch CSV for ${scope}: ${res.statusText}`);
  return await res.text();
}
