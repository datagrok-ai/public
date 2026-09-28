/* Compile-only: the typed domain surface (`DomainTable<TRow>` and what hangs off it) against an
   app's own row interface, and the untyped default beside it. Type-checked by `npm run build`
   (tsconfig includes this folder); never run. */
import {domains, DomainTable, DomainApp} from '../../src/dg/index.js';
import type {DomainRibbon} from '../../src/dg/index.js';
import type {DomainSource, RowView} from '../../src/index.js';

interface IssueRow {
  id: string;
  title: string;
  done: boolean;
  weight?: number;
}

async function typed(): Promise<void> {
  const issues: DomainTable<IssueRow> = await domains.table<IssueRow>('grit.issue');
  const src: DomainSource<IssueRow> = issues.source({query: 'done = false', defaults: {done: false}});
  const draft = issues.draft({title: 'x'});
  const row: RowView<IssueRow> = src.newRow({weight: 2});
  const title: string = row.title;
  const weight: number | undefined = row.weight;
  const other: unknown = row.anything;
  const current: RowView<IssueRow> | null = draft.currentRow.value;
  issues.actions.add({name: 'Close', when: (r) => !r.done, run: (r) => void r.title});
  issues.actions.add({name: 'Escalate', run: (r: IssueRow) => void r.id});
  issues.validators.add('title', (value, r) => r.done && value === '' ? 'A closed issue needs a title' : null);
  issues.validators.add('weight', (_value, r: IssueRow) => r.weight === undefined ? null : 'x');
  issues.renderer = {caption: (r) => r.title, card: (r: IssueRow) => document.createElement(r.done ? 'del' : 'b')};
  // @ts-expect-error a column holds its declared type
  issues.source({defaults: {done: 'yes'}});
  // @ts-expect-error a draft's values too
  issues.draft({title: 1});
  // @ts-expect-error the row is typed in an action
  issues.actions.add({name: 'x', run: (r) => r.title.foo()});
  void [title, weight, other, current];
}

/** A typed table goes wherever an untyped one is expected — the controls have no parameter. */
async function untyped(): Promise<void> {
  const issues = await domains.table<IssueRow>('grit.issue');
  const plain: DomainTable = issues;
  const src: DomainSource = issues.source();
  plain.validators.add('anything', (_value, r) => r.anything === null ? 'x' : null);
  plain.actions.add({name: 'x', run: (r) => void r.whatever});
  const values = plain.source({defaults: {whatever: 1}}).newRow({more: true});
  const anything: unknown = values.more;
  domains.form(src);
  domains.list(issues.source(), {mode: 'cards'});
  domains.grid(src);
  domains.search(src);
  domains.filters(src);
  domains.history(src);
  domains.children(src);
  domains.saveButton(src);
  domains.discardButton(src);
  domains.newButton(src, {title: 'x'});
  class Sub extends DomainApp {
    shortcuts = {'Ctrl+Shift+C': 'Close'};
    ribbon(): DomainRibbon {
      return {...super.ribbon(), presets: [this.presets(['Mine', 'reporter = $me'])]};
    }
  }
  issues.app({app: Sub, shortcuts: {'Delete': 'Delete'}});
  domains.app({table: issues, base: '/x'});
  void anything;
}

void [typed, untyped];
