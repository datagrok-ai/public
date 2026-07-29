import {randomUUID} from 'node:crypto';
import * as api from './api-client.js';

// Operation registry backing the domain tools. Previously each of these was its own MCP tool, so
// all ~32 schemas sat in the model's prompt prefix on every turn (measured: 3,847 tokens/call).
// Collapsing them behind one tool per domain keeps the same reach at a fraction of the prefix:
// the tool description carries only the operation NAMES, and a call with no `op` returns the full
// parameter schema for that domain. See docs/BENCHMARK.md "context diet".

export type ParamType = 'string' | 'number' | 'boolean' | 'array' | 'object';

export interface Param {
  type: ParamType;
  required?: boolean;
  desc: string;
}

export interface Op {
  desc: string;
  params?: Record<string, Param>;
  /** File-ish payload: never array-paged, capped at a larger size than list results. */
  raw?: boolean;
  run: (a: Record<string, any>) => Promise<unknown>;
}

export interface Domain {
  tool: string;
  blurb: string;
  ops: Record<string, Op>;
}

// Read-only by operation name, matching the runtime's verify gate (verify.ts). Naming ops so the
// prefix is meaningful is what lets the gate stay a pure name check after consolidation — the tool
// name alone can no longer say whether a call mutates anything.
const READONLY_RE = /^(whoami$|list|get|search|read|download|find)/;

export const isReadonlyOp = (op: string): boolean => READONLY_RE.test(op);

const str = (desc: string, required = false): Param => ({type: 'string', required, desc});
const bool = (desc: string): Param => ({type: 'boolean', desc});

const FILTER = str('Smart search filter, e.g. name="Chem"');
const LIMIT = {type: 'number' as const, desc: 'Max items to return (default 50)'};
const OFFSET = {type: 'number' as const, desc: 'Items to skip, for paging (default 0)'};
const PAGING = {limit: LIMIT, offset: OFFSET};

export const DOMAINS: Domain[] = [
  {
    tool: 'datagrok_functions',
    blurb: 'Datagrok functions, scripts and data queries — discover, inspect, run, and create them.',
    ops: {
      list: {
        desc: 'List functions matching an optional smart filter.',
        params: {filter: FILTER, ...PAGING},
        run: (a) => api.listFunctions(a.filter),
      },
      get: {
        desc: 'Full details of one function, including its parameter list.',
        params: {id: str('Function ID or name, e.g. "Chem:getMolProperty"', true)},
        run: (a) => api.getFunction(a.id),
      },
      call: {
        desc: 'Execute a function and return its result.',
        params: {
          name: str('Function name, e.g. "Chem:smilesToMw"', true),
          params: {type: 'object', desc: 'Arguments keyed by input name'},
        },
        run: (a) => api.callFunction(a.name, a.params),
      },
      list_scripts: {
        desc: 'List all scripts.',
        params: PAGING,
        run: () => api.listFunctions('source="script"'),
      },
      list_queries: {
        desc: 'List all data queries.',
        params: PAGING,
        run: () => api.listFunctions('source="data-query"'),
      },
      create_script: {
        desc: 'Create a new script.',
        params: {
          script: str('Script source code', true),
          name: str('Script name'),
          language: str('python | r | julia | nodejs | grok | octave (default python)'),
        },
        run: (a) => api.saveFunction({
          source: 'script', script: a.script,
          name: a.name ?? 'Untitled Script', language: a.language ?? 'python',
        }),
      },
      create_query: {
        desc: 'Create a new data query against a connection.',
        params: {
          connectionId: str('ID of the data connection', true),
          query: str('SQL query text', true),
          name: str('Query name'),
        },
        run: (a) => api.saveFunction({
          source: 'data-query', query: a.query,
          connection: {id: a.connectionId}, name: a.name ?? 'Untitled Query',
        }),
      },
    },
  },
  {
    tool: 'datagrok_files',
    blurb: 'Files in connector storage (e.g. "System:DemoFiles") — browse, read, write.',
    ops: {
      list: {
        desc: 'List files in a connector directory.',
        params: {
          connector: str('Connector name, e.g. "System:DemoFiles"', true),
          path: str('Directory path within the connector'), ...PAGING,
        },
        run: (a) => api.listFiles(a.connector, a.path),
      },
      download: {
        desc: 'Read a file\'s contents.',
        params: {connector: str('Connector name', true), path: str('File path within the connector', true)},
        raw: true,
        run: (a) => api.downloadFile(a.connector, a.path),
      },
      upload: {
        desc: 'Write content to a file.',
        params: {
          connector: str('Connector name', true), path: str('File path within the connector', true),
          content: str('File content', true),
        },
        run: (a) => api.uploadFile(a.connector, a.path, a.content),
      },
    },
  },
  {
    tool: 'datagrok_projects',
    blurb: 'Projects — the shareable unit that bundles tables, layouts, queries and connections.',
    ops: {
      list: {desc: 'List projects.', params: {filter: FILTER, ...PAGING}, run: (a) => api.listProjects(a.filter)},
      get: {desc: 'Details of one project.', params: {id: str('Project ID', true)}, run: (a) => api.getProject(a.id)},
      list_recent: {desc: 'List recently accessed projects.', params: PAGING, run: () => api.listRecentProjects()},
      search: {
        desc: 'Find a project by name.',
        params: {name: str('Project name to search for', true), ...PAGING},
        run: (a) => api.searchProject(a.name),
      },
      create: {
        desc: 'Create a new project.',
        params: {name: str('Project name', true), description: str('Project description')},
        run: (a) => api.createProject(a.name, a.description),
      },
      delete: {desc: 'Delete a project.', params: {id: str('Project ID', true)}, run: (a) => api.deleteProject(a.id)},
      attach_entity: {
        desc: 'Attach an existing SERVER-SIDE entity (table, layout, script, query, connection) to a ' +
          'project. Pass a real server id — not a client-side stub from t.getTableInfo().id.',
        params: {
          projectId: str('Project ID (UUID)', true),
          entityId: str('Entity ID (UUID) of the server-side entity', true),
          link: bool('true creates a reference; false (default) means the project owns it'),
        },
        run: (a) => api.attachEntityToProject(a.projectId, a.entityId, a.link ?? false),
      },
      share: {
        desc: 'Share a project with user groups. If the entity is not a project, its owning project ' +
          'is shared. Returns per-group success/already-shared/failure.',
        params: {
          projectId: str('Project ID (UUID) or "namespace:name"', true),
          groups: {type: 'array', required: true, desc: 'Group names, e.g. ["clinical-ops"]'},
          access: str('"View" (default) or "Edit"'),
        },
        run: (a) => api.shareProject(a.projectId, a.groups, a.access === 'Edit' ? 'Edit' : 'View'),
      },
      list_shares: {
        desc: 'Groups a project is shared with, by access level.',
        params: {projectId: str('Project ID (UUID)', true)},
        run: (a) => api.listProjectShares(a.projectId),
      },
    },
  },
  {
    tool: 'datagrok_spaces',
    blurb: 'Spaces — hierarchical containers that group entities and hold their own file storage.',
    ops: {
      list: {desc: 'List root spaces.', params: {filter: FILTER, ...PAGING}, run: (a) => api.listSpaces(a.filter)},
      get: {desc: 'Details of one space.', params: {id: str('Space ID', true)}, run: (a) => api.getSpace(a.id)},
      list_children: {
        desc: 'List a space\'s children — subspaces and contained entities.',
        params: {
          spaceId: str('Space ID', true),
          types: str('Comma-separated entity types, e.g. "Script,DataQuery"'),
          includeLinked: bool('Include linked (non-owned) children'), ...PAGING,
        },
        run: (a) => api.listSpaceChildren(a.spaceId, a.types, a.includeLinked),
      },
      create: {desc: 'Create a root space.', params: {name: str('Space name', true)}, run: (a) => api.createRootSpace(a.name)},
      delete: {desc: 'Delete a space.', params: {id: str('Space ID', true)}, run: (a) => api.deleteSpace(a.id)},
      create_subspace: {
        desc: 'Create a subspace inside a space.',
        params: {
          spaceId: str('Parent space ID', true), name: str('Subspace name', true),
          link: bool('Create as a link reference instead of an owned child'),
        },
        run: (a) => api.createSubspace(a.spaceId, a.name, a.link),
      },
      add_entity: {
        desc: 'Add an entity (script, query, connection, …) to a space.',
        params: {
          spaceId: str('Space ID', true), entityId: str('Entity ID', true),
          link: bool('Add as a link reference instead of an owned child'),
        },
        run: (a) => api.addEntityToSpace(a.spaceId, a.entityId, a.link),
      },
      remove_entity: {
        desc: 'Remove an entity from a space.',
        params: {spaceId: str('Space ID', true), entityId: str('Entity ID', true)},
        run: (a) => api.removeEntityFromSpace(a.spaceId, a.entityId),
      },
      read_file: {
        desc: 'Read a file from a space\'s storage.',
        params: {spaceId: str('Space ID', true), path: str('File path within the space', true)},
        raw: true,
        run: (a) => api.readSpaceFile(a.spaceId, a.path),
      },
      write_file: {
        desc: 'Write a file into a space\'s storage.',
        params: {
          spaceId: str('Space ID', true), path: str('File path within the space', true),
          content: str('File content', true),
        },
        run: (a) => api.writeSpaceFile(a.spaceId, a.path, a.content),
      },
      delete_file: {
        desc: 'Delete a file from a space\'s storage.',
        params: {spaceId: str('Space ID', true), path: str('File path within the space', true)},
        run: (a) => api.deleteSpaceFile(a.spaceId, a.path),
      },
    },
  },
  {
    tool: 'datagrok_platform',
    blurb: 'Who you are, and the instance\'s users, groups and data connections.',
    ops: {
      whoami: {desc: 'The authenticated user.', run: () => api.getCurrentUser()},
      list_users: {desc: 'List users.', params: {filter: FILTER, ...PAGING}, run: (a) => api.listUsers(a.filter)},
      list_groups: {desc: 'List groups.', params: {filter: FILTER, ...PAGING}, run: (a) => api.listGroups(a.filter)},
      list_connections: {
        desc: 'List data connections.',
        params: {filter: FILTER, ...PAGING},
        run: (a) => api.listConnections(a.filter),
      },
    },
  },
];

/** One-line operation menu for a domain — what goes in the tool description. */
export function opMenu(d: Domain): string {
  return Object.keys(d.ops).join(', ');
}

/** Full schema for a domain, returned when a call omits `op` or names an unknown one. */
export function catalog(d: Domain): object {
  const operations: Record<string, object> = {};
  for (const [name, op] of Object.entries(d.ops)) {
    const params: Record<string, string> = {};
    for (const [p, meta] of Object.entries(op.params ?? {}))
      params[p] = `${meta.type}${meta.required ? ', required' : ''} — ${meta.desc}`;
    operations[name] = Object.keys(params).length ? {description: op.desc, params} : {description: op.desc};
  }
  return {tool: d.tool, description: d.blurb, operations};
}

export function missingParams(op: Op, args: Record<string, any>): string[] {
  return Object.entries(op.params ?? {})
    .filter(([name, meta]) => meta.required && args[name] === undefined)
    .map(([name]) => name);
}
