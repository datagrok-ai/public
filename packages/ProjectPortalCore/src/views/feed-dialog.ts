/** The activity feed in the app (F-21): one list renderer shared by the program page's Activity
 * section and the My feed dialog (everything the user follows). Newest first, filtered by kind,
 * paged; an entry over a record traces its provenance, and in My feed the program opens. */
import * as grok from 'datagrok-api/grok';
import {button, divH, divV, Scope, span} from '@datagrok-libraries/u2';
import {countByKind, describeEntry, FEED_KIND_LABELS, FEED_KINDS, FeedEntry, FeedKind, filterFeed, pageFeed} from '../model/feed';
import {feedUserNames, myFeedOf} from '../feed';
import {choice, labelled, linkTo, named, openDialog, PortalDialog} from './u2-controls';

/** Entries per page, on the program page and in My feed. */
export const FEED_PAGE_SIZE = 8;

/** A labelled u2 choice in a filter row; its wrapper carries the name. [onChange] runs after the change has
 * settled, so it may redraw the page that holds the choice. */
export function feedSelect(name: string, label: string, options: {code: string, name: string}[], value: string,
  onChange: (v: string) => void): HTMLElement {
  const input = choice(`${name}-choice`, undefined, options, value, {inline: true,
    onChanged: (v) => queueMicrotask(() => onChange(v))});
  return labelled(name, label, input);
}

/** The kind filter: All, then each kind present, with its count. */
export function kindOptions(entries: FeedEntry[]): {code: string, name: string}[] {
  const counts = countByKind(entries);
  return [{code: '', name: `All (${entries.length})`},
    ...FEED_KINDS.filter((k) => counts.has(k)).map((k) => ({code: k, name: `${FEED_KIND_LABELS[k]} (${counts.get(k)})`}))];
}

export interface FeedListOptions {
  /** Already filtered, newest first. */
  entries: FeedEntry[];
  names: Map<string, string>;
  page: number;
  /** data-u2-name prefix: `feed` on the program page, `my-feed` in the dialog. */
  prefix: string;
  withProgram?: boolean;
  onPage: (page: number) => void;
  /** Traces the record an entry is about. */
  onTrace?: (e: FeedEntry) => void;
  /** Opens the entry's program. */
  onProgram?: (e: FeedEntry) => void;
}

/** One page of entries and the pager. */
export function feedList(o: FeedListOptions): HTMLElement[] {
  if (!o.entries.length)
    return [named(span('Nothing yet', 'pp-muted'), `${o.prefix}-empty`)];
  const p = pageFeed(o.entries, o.page, FEED_PAGE_SIZE);
  const items = p.items.map((e, i) => {
    const parts: HTMLElement[] = [];
    if (o.withProgram && o.onProgram) {
      const link = linkTo(e.program, `${o.prefix}-program-${i}`, () => o.onProgram!(e), {title: 'Open the program'});
      parts.push(link);
    }
    parts.push(span(describeEntry(e, e.who ? o.names.get(e.who) ?? null : null)));
    if (e.ref && e.ref.kind !== 'program' && o.onTrace) {
      const trace = linkTo('Provenance', `${o.prefix}-trace-${i}`, () => o.onTrace!(e),
        {cls: 'pp-gap', title: 'What this rests on, and what rests on it'});
      parts.push(trace);
    }
    return named(divH(parts, 'pp-unit'), `${o.prefix}-entry-${i}`);
  });
  const newer = named(button('Newer', () => o.onPage(p.page - 1)), `${o.prefix}-newer-button`);
  const older = named(button('Older', () => o.onPage(p.page + 1)), `${o.prefix}-older-button`);
  newer.disabled = p.page === 0;
  older.disabled = p.page >= p.pages - 1;
  return [
    named(divV(items, 'pp-units'), `${o.prefix}-entries`),
    ...(p.pages > 1 ? [divH([newer, named(span(`Page ${p.page + 1} of ${p.pages}`, 'pp-muted'), `${o.prefix}-page`), older],
      'pp-actions')] : []),
  ];
}

export interface MyFeedOptions {
  openProgram: (portfolioCode: string, programId: string) => void | Promise<void>;
  trace: (kind: string, id: string) => void;
}

/** My feed: the activity of every program (or series) the current user follows. */
export async function openMyFeedDialog(opts: MyFeedOptions): Promise<PortalDialog> {
  const body = named(divV([span('Reading what you follow…', 'pp-muted')], 'pp-feed'), 'my-feed-content');
  const dlg = openDialog({title: 'My feed', name: 'my-feed-dialog', width: 680}, () => ({content: body}));
  // each draw owns its kind filter; the last one goes with the dialog
  let scope: Scope | null = null;
  dlg.own(() => scope?.dispose());
  try {
    const {entries, programs} = await myFeedOf();
    const names = await feedUserNames(entries);
    let kind = '';
    let page = 0;
    const draw = () => {
      scope?.dispose();
      scope = new Scope();
      const shown = filterFeed(entries, {kinds: kind ? [kind as FeedKind] : null});
      Scope.runWith(scope, () => body.replaceChildren(
        named(span(programs ? `${programs} program${programs === 1 ? '' : 's'} followed; newest first` :
          'You follow nothing yet: Follow a program or a series on its page', 'pp-muted'), 'my-feed-summary'),
        ...(entries.length ? [divH([feedSelect('my-feed-kind-input', 'Show', kindOptions(entries), kind, (v) => {
          kind = v;
          page = 0;
          draw();
        })], 'pp-actions')] : []),
        ...feedList({entries: shown, names, page, prefix: 'my-feed', withProgram: true,
          onPage: (n) => {
            page = n;
            draw();
          },
          onTrace: (e) => opts.trace(e.ref!.kind, e.ref!.id),
          onProgram: (e) => {
            dlg.close();
            void Promise.resolve(opts.openProgram(e.portfolio ?? '', e.programId))
              .catch((x) => grok.shell.error(x?.message ?? String(x)));
          }}),
      ));
    };
    draw();
  } catch (e: any) {
    body.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'my-feed-error'));
  }
  return dlg;
}
