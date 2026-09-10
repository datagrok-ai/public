/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import * as fs from 'fs';
import * as path from 'path';
import {NodeDapi} from '../node-dapi';

export interface Part {name: string; entities?: number; failed?: number; seconds?: number; error?: string}

/**
 * Every space worth migrating on its own: the personal space of each user, and each team space.
 * A whole instance moved in one bundle has every project competing for the same entities, and
 * placement is exclusive on the server — scoping a run to one space keeps that contention inside
 * it, which is the difference between a run that finishes and one that does not.
 */
export async function namespacesOf(dapi: NodeDapi): Promise<string[]> {
  const names = new Set<string>();
  for (const p of await dapi.internal('/projects').listAll({includeRoot: 'true'})) {
    if (p?.isEntity || p?.isPackage) continue;
    const own = String(p.namespace ?? '').split(':')[0];
    if (own) names.add(own);
    else if (p.isRoot && p.name) names.add(String(p.name));
  }
  return [...names].sort();
}

/** Users the source has and the target does not: their content lands under the pusher instead. */
export async function missingUsers(from: NodeDapi, to: NodeDapi): Promise<string[]> {
  const here = new Set((await to.internal('/users').listAll({})).map((u: any) => u.login));
  return (await from.internal('/users').listAll({}))
    .map((u: any) => u.login).filter((login: string) => login && !here.has(login)).sort();
}

/** Packages the source has and the target does not: their functions and connections cannot resolve. */
export async function missingPackages(from: NodeDapi, to: NodeDapi): Promise<string[]> {
  const here = new Set((await to.internal('/packages').listAll({})).map((p: any) => p.name));
  return (await from.internal('/packages').listAll({}))
    .map((p: any) => p.name).filter((name: string) => name && !here.has(name)).sort();
}

/** A run is resumable: what finished is remembered, so a repeat does not redo it. */
export function readState(file: string): Record<string, Part> {
  try { return JSON.parse(fs.readFileSync(file, 'utf8')); }
  catch { return {}; }
}

export function writeState(file: string, state: Record<string, Part>): void {
  fs.mkdirSync(path.dirname(file), {recursive: true});
  fs.writeFileSync(file, JSON.stringify(state, null, 2));
}
