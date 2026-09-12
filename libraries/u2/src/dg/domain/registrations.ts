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
];

export function registerDomainComponents(reg: Registry = globalRegistry): void {
  for (const meta of METAS) {
    if (reg.get(meta.tag) === undefined)
      reg.register(meta);
  }
}
