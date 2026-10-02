/** "My programs" on Home: a tile per pinned program the user can see, across portfolios. A tile
 * says where the program stands on its track, how ready it is for its next gate, what needs a
 * look and whether new data landed; a click opens its page in the Project Portal. The widget reads
 * again whenever Home comes back, so a pin made elsewhere shows without a reload. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, divH, divV, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {pinnedTiles} from '../pinned-programs';
import {describeTile, NOTHING_PINNED, ProgramTile, tileTitle} from '../model/pinning';
import {installStyles} from './styles';
import {openPortal} from './portal-view';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

const statusVariant: Record<string, BadgeVariant> = {
  active: 'success', on_hold: 'warning', terminated: 'error', completed: 'accent',
  transitioned: 'accent', missing_from_feed: 'error',
};

/** Opens the program's page: in a Project Portal view already open, else by its address. */
export async function openProgramPage(t: Pick<ProgramTile, 'portfolioCode' | 'programId' | 'address'>): Promise<void> {
  const open = openPortal();
  if (open) {
    grok.shell.v = open.view;
    await open.portal.showProgram(t.portfolioCode, t.programId);
    return;
  }
  grok.shell.route(t.address);
}

function tile(t: ProgramTile): HTMLElement {
  // the track as segments: done, current, ahead
  const track = divH(t.track.map((s) => {
    const seg = span('', `pp-tile-seg pp-tile-seg-${s.state}`);
    seg.title = s.code;
    return seg;
  }), 'pp-tile-track');
  const readiness = t.readiness === null ? span(t.readinessLine, 'pp-muted') :
    named(badge(t.ready ? 'ready' : `${Math.round(t.readiness * 100)}%`,
      {variant: t.ready ? 'success' : t.readiness === 0 ? 'error' : 'warning'}), 'tile-readiness');
  const el = named(divV([
    divH([span(tileTitle(t), 'pp-card-title'),
      ...(t.pinnedFor.length && !t.pinnedByMe ? [badge(`for ${t.pinnedFor[0]}`, {variant: 'accent'})] : [])], 'pp-card-meta'),
    span(t.portfolioName, 'pp-muted'),
    span(`${t.stage}${t.steps ? ` · ${t.step || '–'} of ${t.steps}` : ''}`, 'pp-tile-stage'),
    ...(t.steps ? [track] : []),
    divH([
      badge(t.status, {variant: statusVariant[t.status] ?? 'default'}),
      ...(t.next ? [span(`→ ${t.next.split(' · ')[0]}`, 'pp-muted')] : []),
      readiness,
      ...(t.attention ? [named(badge(`${t.attention} need attention`, {variant: 'warning'}), 'tile-attention')] : []),
      ...(t.newData ? [named(badge('new data', {variant: 'accent'}), 'tile-new-data')] : []),
    ], 'pp-card-meta'),
  ], 'pp-card pp-tile'), `pinned-tile-${slug(t.code)}`);
  el.title = [describeTile(t), t.readinessLine, ...t.gaps.slice(0, 3).map((g) => `• ${g}`),
    ...(t.newData ? [t.newData] : []), ...(t.pinnedFor.length ? [`Pinned for ${t.pinnedFor.join(', ')}`] : [])].join('\n');
  // after the click is done with: Home takes the click as its own and would stay the current view
  el.onclick = (e) => {
    e.stopPropagation();
    setTimeout(() => void openProgramPage(t).catch((err) => grok.shell.error(err?.message ?? String(err))), 50);
  };
  return el;
}

export class MyProgramsWidget extends DG.Widget {
  private readonly body = named(divV([], 'pp-tiles'), 'my-programs');
  private loading = false;

  constructor() {
    // two columns of the Home grid: room for tiles side by side
    super(divV([], 'pp-my-programs'));
    this.root.classList.add('widget-wide');
    installStyles();
    this.root.append(this.body);
    void this.refresh();
    // Home is kept while other views open; read again each time it comes back
    this.subs.push(grok.events.onCurrentViewChanged.subscribe(() => setTimeout(() => {
      if (this.root.isConnected)
        void this.refresh();
    }, 50)));
  }

  async refresh(): Promise<void> {
    if (this.loading)
      return;
    this.loading = true;
    try {
      const tiles = await pinnedTiles();
      this.body.replaceChildren(...(tiles.length ? tiles.map(tile) :
        [named(span(NOTHING_PINNED, 'pp-muted'), 'no-pinned-programs')]));
    } catch (e: any) {
      this.body.replaceChildren(named(span(`Pinned programs could not be read: ${e?.message ?? e}`, 'pp-reason'), 'pinned-error'));
    } finally {
      this.loading = false;
    }
  }
}
