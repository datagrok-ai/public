import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {COLLATOR} from '../../components/collections/catalog-table.js';

export interface SpaceEntry {
  project: DG.Project;
  parent: SpaceEntry | null;
  children: SpaceEntry[];
  /** Why nothing can be saved here; null when it can. */
  readOnly: string | null;
}

export const READ_ONLY = 'You can view this space but not save to it';
const RECENT_STORE = 'u2-spaces';
const RECENT_KEY = 'recent';
const RECENT_MAX = 5;

export class SpaceTree {
  readonly roots: SpaceEntry[];
  private readonly _byId = new Map<string, SpaceEntry>();

  constructor(roots: SpaceEntry[]) {
    this.roots = roots;
    const walk = (list: SpaceEntry[]) => {
      for (const e of list) {
        this._byId.set(e.project.id, e);
        walk(e.children);
      }
    };
    walk(roots);
  }

  get(id: string): SpaceEntry | null {
    return this._byId.get(id) ?? null;
  }

  /** The friendly names from the root down to `id`; just its own name when it is not in the tree. */
  path(space: DG.Project): string[] {
    const names: string[] = [];
    for (let e = this.get(space.id); e !== null; e = e.parent)
      names.unshift(e.project.friendlyName);
    return names.length > 0 ? names : [space.friendlyName];
  }

  /** `id` and every space under it. */
  under(id: string): SpaceEntry[] {
    const top = this.get(id);
    const out: SpaceEntry[] = [];
    const walk = (e: SpaceEntry) => {
      out.push(e);
      for (const c of e.children)
        walk(c);
    };
    if (top !== null)
      walk(top);
    return out;
  }

  /** Puts a space just created under `parent` (a root space without one) into the tree, in its sorted place. */
  add(project: DG.Project, parent: DG.Project | null): SpaceEntry {
    const up = parent === null ? null : this.get(parent.id);
    const entry: SpaceEntry = {project, parent: up, children: [], readOnly: null};
    const list = up === null ? this.roots : up.children;
    const at = list.findIndex((e) => COLLATOR.compare(e.project.friendlyName, project.friendlyName) > 0);
    list.splice(at < 0 ? list.length : at, 0, entry);
    this._byId.set(project.id, entry);
    return entry;
  }

  all(): SpaceEntry[] {
    return this.roots.flatMap((r) => this.under(r.project.id));
  }

  /** Reads the whole tree. */
  static async load(): Promise<SpaceTree> {
    // the list includes the user's own namespace, which "Not in a space" already stands for
    const own = grok.shell.user?.project?.id;
    const roots = await SpaceTree._entries((await grok.dapi.spaces.list()).filter((p) => p.id !== own), null);
    let level = roots;
    while (level.length > 0) {
      const next = await Promise.all(level.map(async (e) => {
        const kids = await grok.dapi.spaces.id(e.project.id).children.ofTypes('Project').list();
        const spaces = kids.filter((k): k is DG.Project => k instanceof DG.Project && k.isSpace);
        e.children = await SpaceTree._entries(spaces, e);
        return e.children;
      }));
      level = next.flat();
    }
    return new SpaceTree(roots);
  }

  private static async _entries(projects: DG.Project[], parent: SpaceEntry | null): Promise<SpaceEntry[]> {
    const editable = await Promise.all(projects.map((p) => grok.dapi.permissions.check(p, 'Edit').catch(() => false)));
    return projects
      .map((project, i) => ({project, parent, children: [], readOnly: editable[i] ? null : READ_ONLY}))
      .sort((a, b) => COLLATOR.compare(a.project.friendlyName, b.project.friendlyName));
  }
}

/** Creates a space called `name` under `parent`, or a root space without one. */
export function createSpace(name: string, parent: DG.Project | null): Promise<DG.Project> {
  return parent === null ? grok.dapi.spaces.createRootSpace(name) : grok.dapi.spaces.id(parent.id).addSubspace(name);
}

/** The `name` the platform gives an entity called `friendlyName`: its letters and digits, each word capitalized,
 * leading digits dropped ("EGFR SAR review" is EGFRSARReview); prop_gen `friendlyNameToName`. */
export function entityName(friendlyName: string): string {
  let out = '';
  let cap = true;
  for (const c of friendlyName.trim()) {
    if (!/[A-Za-z0-9]/.test(c))
      cap = true;
    else if (out !== '' || !/[0-9]/.test(c)) {
      out += cap ? c.toUpperCase() : c;
      cap = false;
    }
  }
  return out;
}

/** Whether `space` holds a project called `friendlyName`, compared as the platform's move does (by {@link entityName});
 * the project `except` (the one being saved, which the space may link) does not count. */
export async function spaceHasName(space: DG.Project, friendlyName: string, except?: string | null): Promise<boolean> {
  const name = entityName(friendlyName);
  if (name === '' || !await grok.dapi.spaces.id(space.id).subspaceExists(name))
    return false;
  if (!except)
    return true;
  const held = await grok.dapi.spaces.id(space.id).children.ofTypes('Project', true).list();
  return held.some((p) => p.name === name && p.id !== except);
}

/** The ids of the spaces the user saved to last, the latest first. */
export function recentSpaces(): string[] {
  try {
    const ids = JSON.parse(grok.userSettings.getValue(RECENT_STORE, RECENT_KEY) ?? '[]');
    return Array.isArray(ids) ? ids.filter((id): id is string => typeof id === 'string') : [];
  } catch {
    return [];
  }
}

export function rememberSpace(space: DG.Project): void {
  const ids = [space.id, ...recentSpaces().filter((id) => id !== space.id)].slice(0, RECENT_MAX);
  grok.userSettings.add(RECENT_STORE, RECENT_KEY, JSON.stringify(ids));
}
