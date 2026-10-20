/** A program's card: where it stands on its track, how ready it is for its next gate, what needs a
 * look and whether new data landed. One card for every place a program shows: a My programs tile on
 * Home, the context panel (the program's record), the Workspace's preview. */
import * as grok from 'datagrok-api/grok';
import {badge, divH, divV, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {describeTile, ProgramTile, tileTitle} from '../model/pinning';
import {portalIcon} from './icons';
import {openPortal} from './portal-view';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

export const statusVariant: Record<string, BadgeVariant> = {
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

/** Opens the program's portfolio board with the program marked and its card in the context panel. */
export async function showOnBoard(t: Pick<ProgramTile, 'portfolioCode' | 'programId' | 'address'>): Promise<void> {
  const open = openPortal();
  if (open) {
    grok.shell.v = open.view;
    await open.portal.previewProgram(t.portfolioCode, t.programId);
    return;
  }
  grok.shell.route(t.address.slice(0, t.address.lastIndexOf('/')));
}

/** The card of [t], named [name]; [extra] goes at the end of its title line (a favorite star). */
export function programCard(t: ProgramTile, name: string, extra: HTMLElement[] = []): HTMLElement {
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
    divH([portalIcon('program'), span(tileTitle(t), 'pp-card-title pp-grow'),
      ...(t.pinnedFor.length && !t.pinnedByMe ? [badge(`for ${t.pinnedFor[0]}`, {variant: 'accent'})] : []), ...extra], 'pp-card-meta'),
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
  ], 'pp-card pp-tile'), name);
  el.title = [describeTile(t), t.readinessLine, ...t.gaps.slice(0, 3).map((g) => `• ${g}`),
    ...(t.newData ? [t.newData] : []), ...(t.pinnedFor.length ? [`Pinned for ${t.pinnedFor.join(', ')}`] : [])].join('\n');
  return el;
}
