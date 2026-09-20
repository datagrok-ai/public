/* The schema.json and rows the domain tests share — a two-table slice of Grit's, with one column
   of every type the memory backend maps. Not a test file: `node --test tests/*.test.js` skips it. */
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';

export const SCHEMA = {
  name: 'grit',
  tables: {
    project: {
      businessKey: ['key'], friendlyName: 'Projects',
      columns: {
        key: {type: 'string', required: true},
        name: {type: 'string', required: true},
        description: {type: 'string'},
      },
    },
    issue: {
      friendlyName: 'Issues', businessKey: ['project_id', 'number'],
      constraints: {weight_positive: {expr: 'weight >= 0', message: 'Weight is positive'}, legacy: {check: 'number > 0'}},
      permissions: ['escalate'],
      columns: {
        project_id: {type: 'ref', ref: 'project', required: true},
        number: {type: 'int', min: 1},
        title: {type: 'string', required: true, isName: true, searchable: true},
        description: {type: 'string', editor: 'textarea', searchable: true},
        done: {type: 'bool'},
        reporter: {type: 'user'},
        tags: {type: 'string_list'},
        due: {type: 'datetime'},
        weight: {type: 'float'},
        priority: {type: 'string', choices: ['low', 'high'], friendlyName: 'Priority'},
      },
    },
  },
};

export const ROWS = {
  project: [{id: 'p1', key: 'GRIT', name: 'Grit'}, {id: 'p2', key: 'DG', name: 'Datagrok'}],
  issue: [
    {id: 'i1', project_id: 'p1', number: 1, title: 'Aspirin', done: true, tags: ['a'],
      due: '2026-01-01T00:00:00Z', weight: 1.5},
    {id: 'i2', project_id: 'p1', number: 2, title: 'Ibuprofen', done: false, weight: 0.5},
    {id: 'i3', project_id: 'p2', number: 1, title: 'Naproxen', done: false, priority: 'high', weight: 2},
  ],
};

/** A fresh backend over copies of the rows, so a transaction in one test never reaches another. */
export function backend(options = {}) {
  const rows = Object.fromEntries(Object.entries(ROWS).map(([t, list]) => [t, list.map((r) => ({...r}))]));
  return new MemoryDomainBackend(SCHEMA, {rows, ...options});
}

/* A hierarchy table the way the manifest declares one: `hierarchy: true` plus exactly one ref
   column targeting the table itself. */
export const LOCATIONS = {
  name: 'stock',
  tables: {
    location: {
      hierarchy: true, singularName: 'Location',
      columns: {
        name: {type: 'string', required: true, isName: true, searchable: true},
        parent_id: {type: 'ref', ref: 'location'},
        kind: {type: 'string', choices: ['site', 'room', 'shelf']},
      },
    },
  },
};

/** site → room → shelf → box, plus an unrelated second site. */
export const TREE = {location: [
  {id: 'l1', name: 'Site', kind: 'site'},
  {id: 'l2', name: 'Room', parent_id: 'l1', kind: 'room'},
  {id: 'l3', name: 'Shelf', parent_id: 'l2', kind: 'shelf'},
  {id: 'l4', name: 'Box', parent_id: 'l3'},
  {id: 'l9', name: 'Other site', kind: 'site'},
]};

/** A fresh backend over the 3-level tree. */
export function hierarchyBackend(options = {}) {
  return new MemoryDomainBackend(LOCATIONS, {rows: {location: TREE.location.map((r) => ({...r}))}, ...options});
}

/** The `stock.location` table of a fresh {@link hierarchyBackend}. */
export function locations(options = {}) {
  return hierarchyBackend(options).tableSync('stock.location');
}
