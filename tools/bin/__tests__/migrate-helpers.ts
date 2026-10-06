/**
 * Test-only bundle surgery: a copy of a bundle under fresh ids and renamed entities, so it
 * can be pushed back to the server it was pulled from (Ruling R6). Not part of the CLI.
 */
import {randomUUID} from 'crypto';
import {BundleEntity} from '../utils/migrate/bundle';
import {nestedIds, rewrite} from '../utils/migrate/rewriter';

/**
 * A fresh identity for every entity and every nested row, with all cross-references
 * rewritten. `rename` also renames the group names grants and memberships refer to.
 */
export function reId(entities: Map<string, BundleEntity>, rename?: (name: string) => string,
): {entities: Map<string, BundleEntity>; idmap: Record<string, string>} {
  const idmap: Record<string, string> = {};
  for (const id of entities.keys())
    idmap[id] = randomUUID();
  for (const {json} of entities.values())
    for (const nested of nestedIds(json))
      idmap[nested] = randomUUID();

  const out = new Map<string, BundleEntity>();
  for (const [id, e] of entities) {
    const json = rewrite(e.json, idmap);
    if (rename)
      renameIn(json, rename);
    out.set(idmap[id], {type: e.type, json});
  }
  return {entities: out, idmap};
}

/** Renames each segment, so `Admin:MigTestSpace:` follows its renamed space. */
export const renameNqName = (nqName: string, rename: (name: string) => string): string =>
  nqName.split(':').map((s) => (s ? rename(s) : s)).join(':');

export function renameIn(json: any, rename: (name: string) => string): void {
  for (const key of ['name', 'friendlyName'])
    if (typeof json[key] === 'string')
      json[key] = rename(json[key]);
  if (typeof json.namespace === 'string')
    json.namespace = renameNqName(json.namespace, rename);
  // A file's path is its directory plus its name, and the server derives it back from both.
  if (typeof json.path === 'string') {
    const segments = json.path.split('/');
    segments[segments.length - 1] = rename(segments[segments.length - 1]);
    json.path = segments.join('/');
  }
  for (const grant of json._grants ?? [])
    grant.group = rename(grant.group);
  for (const member of json._members ?? [])
    if (member.kind === 'group')
      member.name = rename(member.name);
}
