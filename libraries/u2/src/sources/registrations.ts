/* The four platform-backed v1 sources (DD14), on the same tray `u2-state` lives on. All
   `visual: false` and `category: 'Data'`, so the palette groups them into one pane; with no
   platform behind them they render the standard broken-chip containment, which is exactly what
   the gallery and the headless tests see. */
import {FuncSource} from './func-source.js';
import {EntitySource} from './entity-source.js';
import {TableSource} from './table-source.js';
import {EntityRef} from './entity-ref.js';
import {ComponentMeta, Registry, registry as globalRegistry} from '../spec/registry.js';

/** The `grok.dapi.*` collections a spec may list — every one of them a paged entity source. */
export const COLLECTIONS = ['users', 'groups', 'projects', 'queries', 'connections', 'scripts',
  'spaces'];

const DESIGN_DATA = {
  name: 'designData', type: 'string', choices: ['live', 'schema', 'sample'],
  description: 'What the designer shows: real runs, the declared shape without running, or ' +
    'the stored `sample`.',
};

const SAMPLE = {
  name: 'sample', type: 'object',
  description: 'The preview the `sample` policy builds its outputs from; rows become a table.',
};

const METAS: ComponentMeta[] = [
  {
    tag: 'u2-func-source',
    category: 'Data',
    visual: false,
    createComponent: (props, env) => new FuncSource(props, env),
    description: 'A function, query or script run as data: its outputs are what the form binds to.',
    usage: 'The default source for anything computed or queried. Bind a param to an input and the ' +
      'source re-runs as that input changes; bind a DataFrame output\'s `currentRow` to reach ' +
      'the fields of the current record.',
    props: [
      {name: 'func', type: 'string',
        description: 'The function name, `Package:Name` for a package function.'},
      {name: 'params', type: 'object', subBindable: true,
        description: 'What it is called with; a `params.<name>` bind makes that param reactive.'},
      {name: 'auto', type: 'bool',
        description: 'Run at construction and whenever a bound param changes (default true).'},
      {name: 'debounce', type: 'int',
        description: 'Milliseconds a param change waits before the re-run (default 300).'},
      DESIGN_DATA,
      SAMPLE,
    ],
    designPreview: [{name: 'Sample row', value: 1}],
    example: {tag: 'u2-func-source', name: 'orders', props: {func: 'U2demo:demoOrders'}},
  },
  {
    tag: 'u2-entity-source',
    category: 'Data',
    visual: false,
    createComponent: (props, env) => new EntitySource(props, env),
    description: 'A paged listing of a server collection: users, groups, projects and the rest.',
    usage: 'For a picker over server entities, or a list of them. `names` feeds a choice input ' +
      'directly; `items` carries the entities themselves.',
    props: [
      {name: 'entity', type: 'string', choices: COLLECTIONS,
        description: 'The `grok.dapi` collection to list.'},
      {name: 'filter', type: 'string', bindable: true,
        description: 'A smart-search filter; changing it reloads the collection from the first page.'},
      {name: 'order', type: 'string', description: 'The field to order by.'},
      {name: 'pageSize', type: 'int', description: 'How many entities one page loads (default 20).'},
    ],
    defaults: {entity: 'users', pageSize: 20},
    example: {tag: 'u2-entity-source', name: 'people', props: {entity: 'users', pageSize: 20}},
  },
  {
    tag: 'u2-table-source',
    category: 'Data',
    visual: false,
    createComponent: (props) => new TableSource(props),
    description: 'An open workspace table, by name — the whole DataFrame binding surface.',
    usage: 'For a form over a table the user already has open. It references the table BY NAME, ' +
      'so prefer a query or a project table for a form that must keep working tomorrow.',
    props: [
      {name: 'table', type: 'string', description: 'The name of the open table.'},
    ],
    example: {tag: 'u2-table-source', name: 'demog', props: {table: 'demog'}},
  },
  {
    tag: 'u2-entity-ref',
    category: 'Data',
    visual: false,
    createComponent: (props, env) => new EntityRef(props, env),
    description: 'One server entity, looked up by id — the object a details form shows.',
    usage: 'For a form about a single entity. Bind `id` to whatever selects it and every bound ' +
      'field follows; the entity\'s own properties are binding steps.',
    props: [
      {name: 'entityType', type: 'string', choices: COLLECTIONS,
        description: 'Which collection the id belongs to.'},
      {name: 'id', type: 'string', bindable: true, description: 'The entity id to look up.'},
      DESIGN_DATA,
      SAMPLE,
    ],
    defaults: {entityType: 'users', designData: 'live'},
    designPreview: {name: 'Sample entity', friendlyName: 'Sample entity'},
    example: {tag: 'u2-entity-ref', name: 'author', props: {entityType: 'users', id: 'current'}},
  },
];

/** Called once per registry by `registerAll`, whose WeakSet is what makes a repeated import safe —
 * `Registry.register` itself refuses a tag that is already there. The manifest, the gallery and the
 * designer palette all see the sources alongside every visual tag. */
export function registerDataSources(reg: Registry = globalRegistry): void {
  for (const meta of METAS)
    reg.register(meta);
}
