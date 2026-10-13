// `DomainSource` is EMS, so its tag registers here rather than in u2's `registerDataSources`
import {ComponentMeta, Registry, registry as globalRegistry} from '@datagrok-libraries/u2';
import {DomainSource} from './domain-source.js';

const METAS: ComponentMeta[] = [
  {
    tag: 'u2-domain-source',
    category: 'Data',
    visual: false,
    createComponent: (props, env) => new DomainSource(props, env),
    description: 'A domain (EMS) table as data: its rows, the current row a form edits, the ' +
      'caller\'s access, and every pending change until `save`.',
    usage: 'The one source every domain control binds to. Bind `query` to a filter input and ' +
      '`search` to a search box to narrow the rows; bind a form\'s inputs under `currentRow` to edit ' +
      'the current row; wire Save and Discard to the `save` and `discard` functions (they go through ' +
      'the session), and `newRow` to a Create button; `draft: true` is the source of a create form. ' +
      'Every source in one spec shares the ambient session, so one Save writes them all as one ' +
      'transaction. Prefer it over a query source for anything a user edits.',
    props: [
      {name: 'table', type: 'string', description: 'The table address, `<schema>.<table>`.'},
      {name: 'query', type: 'string', bindable: true,
        description: 'A smart-filter string; changing it reloads the rows from the first page.'},
      {name: 'search', type: 'string', bindable: true,
        description: 'A case-insensitive text over the table\'s searchable columns, ANDed with the query.'},
      {name: 'pageSize', type: 'int', description: 'How many rows one page loads (default 50).'},
      {name: 'withAccess', type: 'bool',
        description: 'Fetch the per-row access columns with every row (default true); the table-level ' +
          'access is always fetched.'},
      {name: 'captions', type: 'string_list',
        description: 'Ref columns whose target names ride with the rows; every visible ref column by default.'},
      {name: 'defaults', type: 'object',
        description: 'Column values every draft row starts with — a parent\'s id on a child table.'},
      {name: 'empty', type: 'bool',
        description: 'Load no rows; drafts can still be added — a child collection under a draft parent.'},
      {name: 'draft', type: 'bool',
        description: 'Load nothing and start on one pristine draft — what a create form binds to.'},
      {name: 'deleted', type: 'string', choices: ['exclude', 'include', 'only'],
        description: 'Which rows the source answers: the live ones (default), the live and the ' +
          'soft-deleted, or the deleted alone — a trash list, read-only until its rows are restored.'},
      {name: 'live', type: 'bool', bindable: true,
        description: 'Follow the server: the source probes the table every 30 s and reloads while there ' +
          'is nothing unsaved to lose, marking itself stale while there is.'},
    ],
    defaults: {pageSize: 50, withAccess: true, empty: false, draft: false, deleted: 'exclude', live: false},
    example: {tag: 'u2-domain-source', name: 'issues', props: {table: 'grit.issue', pageSize: 50}},
  },
];

export function registerDomainSources(reg: Registry = globalRegistry): void {
  for (const meta of METAS) {
    if (reg.get(meta.tag) === undefined)
      reg.register(meta);
  }
}
