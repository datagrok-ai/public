/** Pinned programs: a program in favorites is its platform record (`portal.program`) among the
 * user's favorites or a group's. The favorite is the only record of it, so it shows in the
 * Workspace on Home under the program's name, and a favorite added or removed there, or with the
 * star, shows here and the other way round. Before 0.3.0 a pin was a favorite on the program's
 * workspace folder; such a favorite still counts, and removing the program removes both. A
 * favorite program also gets a tile in the "My programs" widget on Home: where it stands, how ready
 * it is for its next gate, and what needs a look.
 *
 * Pure: the caller reads the favorites, the programs' facts and the follows. */
import {isReady} from './readiness';

/** The app's address; a program's page is `<portfolio>/<program>` under it. */
export const PORTAL_ADDRESS = '/apps/ProjectPortalCore/Portal';

export const programAddress = (portfolioCode: string, programCode: string): string =>
  `${PORTAL_ADDRESS}/${encodeURIComponent(portfolioCode)}/${encodeURIComponent(programCode)}`;

/** What the widget says when nothing is pinned. */
export const NOTHING_PINNED = 'No programs in your favorites. Open a program in the Project Portal and click ' +
  '"Add to favorites": it shows here, and in your Workspace on Home.';

/** The entity name of a program's platform record: `portal.program:<portfolio id>|<code>`. */
export const PROGRAM_RECORD = 'portal.program';

/** Whether an entity named [name] is a program's record (what a favorite program is). */
export const isProgramRecordName = (name: string | null | undefined): boolean =>
  String(name ?? '').toLowerCase().startsWith(`${PROGRAM_RECORD}:`);

/** Why [p] cannot be added to favorites, null when it can: the favorite is the program's record,
 * which the platform must find as an entity. No workspace folder is needed. */
export function whyNotPinnable(p: {recordFound: boolean}): string | null {
  return p.recordFound ? null : 'Its record is not found as a platform entity, so it cannot be a favorite';
}

/** Why the user may not pin for a group, empty when they may: the portfolio's approvers. */
export function checkGroupPin(c: {isApprover: boolean}): string[] {
  return c.isApprover ? [] : ['Only the portfolio\'s approvers pin a program for a group'];
}

/** Who pinned a program: the user, and the groups that have it among their favorites. */
export interface PinHolders {
  mine: boolean;
  groups: string[];
}

/** A favorite as read: its id, and whether it is a program's record or a folder (a pin made
 * before 0.3.0); null for anything else. */
export interface FavoriteRef {
  id: string;
  kind: 'program' | 'folder' | null;
}

/** The pins by entity id (a program record's id is the program's id; a folder's is its space id),
 * from the user's favorites and those of their groups. Other favorites are no pins. */
export function pinsOf(personal: FavoriteRef[], byGroup: {group: string, entities: FavoriteRef[]}[]): Map<string, PinHolders> {
  const pins = new Map<string, PinHolders>();
  const at = (id: string) => pins.get(id) ?? pins.set(id, {mine: false, groups: []}).get(id)!;
  for (const e of personal.filter((x) => x.kind))
    at(e.id).mine = true;
  for (const g of byGroup) {
    for (const e of g.entities.filter((x) => x.kind)) {
      const h = at(e.id);
      if (!h.groups.includes(g.group))
        h.groups.push(g.group);
    }
  }
  return pins;
}

/** Who pinned [p]: its record's favorites and its folder's together; undefined when neither is. */
export function pinOf(pins: Map<string, PinHolders>, p: {programId: string, spaceId?: string | null}): PinHolders | undefined {
  const held = [pins.get(p.programId), p.spaceId ? pins.get(p.spaceId) : undefined].filter((h): h is PinHolders => !!h);
  if (!held.length)
    return undefined;
  return {mine: held.some((h) => h.mine), groups: [...new Set(held.flatMap((h) => h.groups))]};
}

/** One pinned program as read: the plain facts a tile is made of. */
export interface TileInput {
  programId: string;
  code: string;
  name: string;
  spaceId: string | null;
  portfolioCode: string;
  portfolioName: string;
  status: string | null;
  /** The program's stage model in order, and its stage there; null when unstaged. */
  stages: {code: string, name: string}[];
  stage: string | null;
  next: {code: string, name: string} | null;
  /** score is null when the next gate has no required criteria. */
  readiness: {score: number | null, requiredMet: number, requiredTotal: number, gaps: string[], dataError?: string} | null;
  /** What trips the portfolio's attention rules; null when it declares none or the data could not be read. */
  attention: number | null;
  /** The user's active follows of the program and its work units. */
  follows: {newsOn: string | null, seenOn: string | null, lastNews: string | null}[];
}

export interface ProgramTile {
  programId: string;
  code: string;
  name: string;
  portfolioCode: string;
  portfolioName: string;
  status: string;
  /** `KIN · Lead optimization`, or Unstaged. */
  stage: string;
  /** 1-based place on the track and the track's length; 0 when unstaged. */
  step: number;
  steps: number;
  /** The stage codes, the current one marked done/current/ahead, for the track. */
  track: {code: string, state: 'done' | 'current' | 'ahead'}[];
  next: string | null;
  /** 0..1, null when there is no next gate or it has no criteria. */
  readiness: number | null;
  ready: boolean;
  /** `3 of 4 required met`, or why there is nothing to measure. */
  readinessLine: string;
  gaps: string[];
  attention: number | null;
  followed: boolean;
  /** The last news line not seen yet; null when nothing new or not followed. */
  newData: string | null;
  pinnedByMe: boolean;
  pinnedFor: string[];
  address: string;
}

const iso = (v: string | null) => v ?? '';

/** The news a program's follows hold that the user has not seen: landed after their last visit. */
export function unseenNews(follows: TileInput['follows']): string | null {
  const unseen = follows.filter((f) => f.newsOn && f.lastNews && iso(f.newsOn) > iso(f.seenOn))
    .sort((a, b) => iso(b.newsOn).localeCompare(iso(a.newsOn)));
  return unseen[0]?.lastNews ?? null;
}

/** The tiles of the pinned programs, by portfolio name and then code. A program no favorite holds
 * (neither its record nor its folder) gets no tile, so the tiles follow the favorites exactly. */
export function assembleTiles(inputs: TileInput[], pins: Map<string, PinHolders>,
  opts: {unpinned?: boolean} = {}): ProgramTile[] {
  const tiles: ProgramTile[] = [];
  for (const x of inputs) {
    const found = pinOf(pins, x);
    const pinned = !!found && (found.mine || found.groups.length > 0);
    // a program's card shows wherever the program does; My programs keeps the pinned ones
    if (!pinned && !opts.unpinned)
      continue;
    const held = found ?? {mine: false, groups: []};
    const at = x.stage ? x.stages.findIndex((s) => s.code === x.stage) : -1;
    const r = x.readiness;
    tiles.push({
      programId: x.programId, code: x.code, name: x.name, portfolioCode: x.portfolioCode, portfolioName: x.portfolioName,
      status: x.status ?? 'active',
      stage: at < 0 ? 'Unstaged' : `${x.stages[at].code} · ${x.stages[at].name}`,
      step: at + 1, steps: x.stages.length,
      track: x.stages.map((s, i) => ({code: s.code, state: at < 0 || i > at ? 'ahead' : i === at ? 'current' : 'done'})),
      next: x.next ? `${x.next.code} · ${x.next.name}` : null,
      readiness: r?.score ?? null,
      ready: isReady(r),
      readinessLine: !x.next ? 'No further gate' : !r ? 'Readiness not computed' :
        r.score === null ? 'No criteria defined' : `${r.requiredMet} of ${r.requiredTotal} required met`,
      gaps: r?.gaps ?? [],
      attention: x.attention,
      followed: x.follows.length > 0,
      newData: unseenNews(x.follows),
      pinnedByMe: held.mine,
      pinnedFor: [...held.groups].sort(),
      address: programAddress(x.portfolioCode, x.code),
    });
  }
  return tiles.sort((a, b) => a.portfolioName.localeCompare(b.portfolioName) || a.code.localeCompare(b.code));
}

/** `EXS-KIN1 · Kinase`; the code alone when the name is the code. */
export const tileTitle = (t: Pick<ProgramTile, 'code' | 'name'>): string =>
  t.name && t.name !== t.code ? `${t.code} · ${t.name}` : t.code;

/** A tile as one line: what the widget's tooltip and the pinnedPrograms table say. */
export function describeTile(t: ProgramTile): string {
  return [`${tileTitle(t)} (${t.portfolioCode})`, t.stage + (t.steps ? ` [${t.step}/${t.steps}]` : ''), t.status,
    t.next ? `${t.next}: ${t.readiness === null ? t.readinessLine : `${Math.round(t.readiness * 100)}%${t.ready ? ' ready' : ''}`}` :
      t.readinessLine,
    ...(t.attention ? [`${t.attention} need attention`] : []),
    ...(t.newData ? ['new data'] : [])].join(' · ');
}
