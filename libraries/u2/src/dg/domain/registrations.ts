/* The `u2-domain-*` tags (ruling 8: every domain control is spec-registered from day one). A
   control binds to its source as a whole — `bind: {source: "$.issues.source"}` — and inherits the
   source's access; registered by `registerPlatformComponents`, next to the other dg tags. */
import {Signal} from '../../core/signals.js';
import type {LiveOption} from '../../core/input-base.js';
import {Registry, registry as globalRegistry} from '../../spec/registry.js';
import type {ComponentMeta} from '../../spec/registry.js';
import {DomainSource} from '../../sources/domain-source.js';
import type {FormLayout} from '../../components/forms/form.js';
import {DomainForm} from './form.js';
import {DomainList} from './list.js';
import type {DomainListMode} from './list.js';
import {DomainPick} from './pick.js';
import {DomainGrid} from './grid.js';
import {DomainHistory} from './history.js';
import {DomainChildren} from './children.js';
import type {DomainChildrenMode} from './children.js';
import {DomainSearch} from './search.js';
import {DomainFilters} from './filters.js';
import type {DomainFiltersMode} from './filters.js';

type Props = Record<string, unknown>;

const SOURCE = {name: 'source', type: 'object', bindable: true,
  description: 'The domain source — bind `$.<source>.source`; the control inherits its access.'};

/** The source a bound `source` prop carries: the source itself, or a signal holding it. */
export function domainSourceProp(x: unknown, tag: string): DomainSource {
  const value = x instanceof Signal ? x.peek() : x;
  if (value instanceof DomainSource)
    return value;
  throw new Error(`${tag}: bind "source" to a domain source ("$.<source>.source")`);
}

function list(x: unknown): string[] | undefined {
  return Array.isArray(x) ? x.map(String) :
    typeof x === 'string' && x !== '' ? x.split(',').map((s) => s.trim()) : undefined;
}

const METAS: ComponentMeta[] = [
  {
    tag: 'u2-domain-form',
    category: 'Inputs',
    create: (props: Props) => new DomainForm(domainSourceProp(props.source, 'u2-domain-form'), {
      include: list(props.include), exclude: list(props.exclude), layout: props.layout as FormLayout | undefined,
    }),
    description: 'A form over the current row of a domain source: one editor per column the caller may ' +
      'edit, text for what they may only read, every edit tracked by the source until Save.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.issues.source`) and the form edits its current row — ' +
      'what a list over the same source selects, or the draft a `draft` source holds. Pair it with a Save button ' +
      'wired to `cmd:issues.save`. Prefer it over hand-bound inputs: the schema, the access and the ' +
      'validation come from the table — a column the caller may not see is absent, one they may not ' +
      'write (no edit on the row, no insert for a draft) is text.',
    props: [
      SOURCE,
      {name: 'include', type: 'string_list', description: 'Columns to show, in this order; every column by default.'},
      {name: 'exclude', type: 'string_list', description: 'Columns to leave out.'},
      {name: 'layout', type: 'string', choices: ['auto', 'normal', 'wide', 'tall']},
    ],
    example: {tag: 'u2-domain-form', bind: {source: '$.issues.source'}},
  },
  {
    tag: 'u2-domain-list',
    category: 'Display',
    create: (props: Props) => new DomainList(domainSourceProp(props.source, 'u2-domain-list'), {
      mode: props.mode as DomainListMode | undefined,
      itemHeight: props.itemHeight as number | undefined,
      empty: props.empty as string | undefined,
    }),
    description: 'The rows of a domain source as a virtual list — the table\'s own rendering, per-row ' +
      'actions gated by access, the next page loaded on scroll; selecting a row makes it the source\'s ' +
      'current row.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.issues.source`). `brief` shows one line per row, ' +
      '`cards` the table\'s card. Row actions follow the access: Delete needs the delete capability, taken ' +
      'per row from the `~can_*` columns the source fetches with `withAccess`. Put a `u2-domain-form` over ' +
      'the same source beside it for master–detail.',
    props: [
      SOURCE,
      {name: 'mode', type: 'string', choices: ['brief', 'cards'], description: 'One line per row, or a card.'},
      {name: 'itemHeight', type: 'int', description: 'Row height in pixels (28 brief, 96 cards).'},
      {name: 'empty', type: 'string', description: 'What an empty result says.'},
    ],
    defaults: {mode: 'brief'},
    example: {tag: 'u2-domain-list', bind: {source: '$.issues.source'}, props: {mode: 'cards'}},
  },
  {
    tag: 'u2-domain-pick',
    category: 'Inputs',
    create: (props: Props) => {
      const value = props.value;
      const bound = value instanceof Signal;
      return new DomainPick(String(props.table ?? ''), {
        label: props.label as LiveOption<string> | undefined,
        name: props.name as string | undefined,
        tooltipText: props.tooltipText as LiveOption<string> | undefined,
        enabled: props.enabled as LiveOption<boolean> | undefined,
        nullable: props.nullable as boolean | undefined,
        placeholder: props.placeholder as string | undefined,
        filter: props.filter as string | undefined,
        value: bound ? undefined : value as string | null | undefined,
        bind: bound ? value as Signal<string | null> : undefined,
      });
    },
    description: 'Picks a row of a domain table by typing its name; the value is the row id — what a ' +
      'reference column holds.',
    usage: 'Use for a `ref` column or any place a row of another table is chosen: `table` names it ' +
      '(`<schema>.<table>`), `filter` narrows the candidates with a smart filter. Bind `value` two-way to ' +
      'the column (`$.issues.currentRow.project_id`); a generated domain form does this by itself.',
    props: [
      {name: 'table', type: 'string', description: 'The table to pick from, `<schema>.<table>`.'},
      {name: 'label', type: 'string', bindable: true},
      {name: 'name', type: 'string', description: 'Stable key for forms and dumps; defaults to the label.'},
      {name: 'value', type: 'string', bindable: true, twoWay: true, description: 'The picked row\'s id.'},
      {name: 'placeholder', type: 'string'},
      {name: 'filter', type: 'string', description: 'A smart filter the candidates must match.'},
      {name: 'nullable', type: 'bool'},
      {name: 'tooltipText', type: 'string', bindable: true},
      {name: 'enabled', type: 'bool', bindable: true},
    ],
    example: {tag: 'u2-domain-pick', props: {label: 'Project', table: 'grit.project'}},
  },
  {
    tag: 'u2-domain-grid',
    category: 'Display',
    appearance: false,
    create: (props: Props) => new DomainGrid(domainSourceProp(props.source, 'u2-domain-grid')),
    description: 'The platform grid over the rows of a domain source: in-cell editing goes through the ' +
      'source\'s writer (pending cells amber, invalid red, conflicts orange), service columns hidden, the ' +
      'table\'s own column decoration; the grid\'s current row and selection are the source\'s.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.issues.source`). Prefer it over `u2-domain-list` for ' +
      'many columns or bulk edits — every cell is editable in place, under the row\'s access. Save and ' +
      'Discard are the session\'s, not the grid\'s: wire them to `cmd:issues.save` / `cmd:issues.discard`.',
    props: [SOURCE],
    example: {tag: 'u2-domain-grid', bind: {source: '$.issues.source'}},
  },
  {
    tag: 'u2-domain-history',
    category: 'Display',
    create: (props: Props) => new DomainHistory(domainSourceProp(props.source, 'u2-domain-history')),
    description: 'The audit trail of the source\'s current row, newest first: who did what when, and the ' +
      'columns an update changed as `caption: before → after`. A draft says "Not saved yet"; refreshed ' +
      'when the session saves.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.issues.source`) and put it on the entity page, beside ' +
      'or under the `u2-domain-form` over the same source — it follows the row the form edits. Only for ' +
      'tables whose schema keeps an audit trail.',
    props: [SOURCE],
    example: {tag: 'u2-domain-history', bind: {source: '$.issues.source'}},
  },
  {
    tag: 'u2-domain-children',
    category: 'Display',
    create: (props: Props) => new DomainChildren(domainSourceProp(props.source, 'u2-domain-children'), {
      tables: list(props.tables), mode: props.mode as DomainChildrenMode | undefined,
    }),
    description: 'One tab per table that refers to the source\'s rows, over its current row: the child ' +
      'rows queried by the foreign key, New pre-filled with it, every child in the same session — a draft ' +
      'parent and its child drafts save as one transaction.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.projects.source`) on the entity page, under the ' +
      '`u2-domain-form`. `tables` narrows the tabs; `grid` (default) edits the children in place, `list` ' +
      'pairs a list with a form. Save and Discard are the session\'s — one Save lands the parent and its ' +
      'children.',
    props: [
      SOURCE,
      {name: 'tables', type: 'string_list', description: 'The child tables to show; every one by default.'},
      {name: 'mode', type: 'string', choices: ['grid', 'list'],
        description: 'The platform grid, or a list beside a form.'},
    ],
    defaults: {mode: 'grid'},
    example: {tag: 'u2-domain-children', bind: {source: '$.projects.source'}},
  },
  {
    tag: 'u2-domain-search',
    category: 'Inputs',
    create: (props: Props) => new DomainSearch(domainSourceProp(props.source, 'u2-domain-search'), {
      placeholder: props.placeholder as string | undefined,
      debounceMs: props.debounceMs as number | undefined,
    }),
    description: 'A search box over a domain source: the text goes to the table\'s searchable columns, ' +
      'AND-ed with the source\'s query, after a short pause or on Enter.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.issues.source`) and put it in the ribbon beside ' +
      '`u2-domain-filters`. It writes the source\'s `search`; Escape clears it. The columns it searches are ' +
      'the schema\'s `searchable` ones (the name column by default) — no `like` chains in the app.',
    props: [
      SOURCE,
      {name: 'placeholder', type: 'string'},
      {name: 'debounceMs', type: 'int', description: 'The pause before the search runs (300 by default).'},
    ],
    example: {tag: 'u2-domain-search', bind: {source: '$.issues.source'}},
  },
  {
    tag: 'u2-domain-filters',
    category: 'Inputs',
    create: (props: Props) => new DomainFilters(domainSourceProp(props.source, 'u2-domain-filters'), {
      mode: props.mode as DomainFiltersMode | undefined,
      placeholder: props.placeholder as string | undefined,
    }),
    description: 'The filter query box (or the condition builder) over a domain source\'s query, two-way, ' +
      'with the table\'s columns and values completed.',
    usage: 'Bind `source` to a `u2-domain-source` (`$.issues.source`). `query` is the one-line smart-filter ' +
      'box with completion, `builder` a row per condition. A change while the session holds unsaved edits ' +
      'asks to save or discard them first; the query round-trips to the URL (`?q=`) through the app\'s path.',
    props: [
      SOURCE,
      {name: 'mode', type: 'string', choices: ['query', 'builder']},
      {name: 'placeholder', type: 'string'},
    ],
    defaults: {mode: 'query'},
    example: {tag: 'u2-domain-filters', bind: {source: '$.issues.source'}},
  },
];

export function registerDomainComponents(reg: Registry = globalRegistry): void {
  for (const meta of METAS) {
    if (reg.get(meta.tag) === undefined)
      reg.register(meta);
  }
}
