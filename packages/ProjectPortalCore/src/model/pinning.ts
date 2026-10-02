/** Pinned programs: pinning a program pins its workspace folder. The platform favorite on
 * the folder is the only record, so a pin made or removed in the Workspace, or with the star,
 * shows here and the other way round. A pinned program also gets a tile in the "My programs"
 * widget on Home: where it stands, how ready it is for its next gate, and what needs a look.
 *
 * Pure: the caller reads the favorites, the programs' facts and the follows. */

/** The app's address; a program's page is `<portfolio>/<program>` under it. */
export const PORTAL_ADDRESS = '/apps/ProjectPortalCore/Portal';

export const programAddress = (portfolioCode: string, programCode: string): string =>
  `${PORTAL_ADDRESS}/${encodeURIComponent(portfolioCode)}/${encodeURIComponent(programCode)}`;

/** What the widget says when nothing is pinned. */
export const NOTHING_PINNED = 'No programs pinned. Open a program in the Project Portal and click "Pin to my Workspace": ' +
  'it shows here, and its folder in your Workspace.';

/** Why [p] cannot be pinned, null when it can: only a folder can be a favorite, so a program
 * without a built workspace folder has nothing to pin. */
export function whyNotPinnable(p: {spaceId: string | null | undefined, spaceLive: boolean}): string | null {
  if (!p.spaceId)
    return 'It has no workspace folder yet; pinning pins the folder. An approver builds it (buildWorkspace)';
  if (!p.spaceLive)
    return 'Its workspace folder is gone; rebuild the portfolio\'s workspace to pin it';
  return null;
}

/** Why the user may not pin for a group, empty when they may: the portfolio's approvers. */
export function checkGroupPin(c: {isApprover: boolean}): string[] {
  return c.isApprover ? [] : ['Only the portfolio\'s approvers pin a program for a group'];
}

/** Who pinned a folder: the user, and the groups that have it among their favorites. */
export interface PinHolders {
  mine: boolean;
  groups: string[];
}

/** The pinned folders by id, from the user's favorites and those of their groups, any entity. */
export function pinsOf(personal: {id: string, isSpace: boolean}[], byGroup: {group: string, entities: {id: string, isSpace: boolean}[]}[]):
  Map<string, PinHolders> {
  const pins = new Map<string, PinHolders>();
  const at = (id: string) => pins.get(id) ?? pins.set(id, {mine: false, groups: []}).get(id)!;
  for (const e of personal.filter((x) => x.isSpace))
    at(e.id).mine = true;
  for (const g of byGroup) {
    for (const e of g.entities.filter((x) => x.isSpace)) {
      const h = at(e.id);
      if (!h.groups.includes(g.group))
        h.groups.push(g.group);
    }
  }
  return pins;
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
  readiness: {score: number, requiredMet: number, requiredTotal: number, gaps: string[], dataError?: string} | null;
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
  /** 0..1, null when there is no next gate. */
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

/** The tiles of the pinned programs, by portfolio name and then code. A program whose folder no
 * pin holds gets no tile, so the tiles follow the favorites exactly. */
export function assembleTiles(inputs: TileInput[], pins: Map<string, PinHolders>): ProgramTile[] {
  const tiles: ProgramTile[] = [];
  for (const x of inputs) {
    const held = x.spaceId ? pins.get(x.spaceId) : undefined;
    if (!held || (!held.mine && !held.groups.length))
      continue;
    const at = x.stage ? x.stages.findIndex((s) => s.code === x.stage) : -1;
    const r = x.readiness;
    tiles.push({
      programId: x.programId, code: x.code, name: x.name, portfolioCode: x.portfolioCode, portfolioName: x.portfolioName,
      status: x.status ?? 'active',
      stage: at < 0 ? 'Unstaged' : `${x.stages[at].code} · ${x.stages[at].name}`,
      step: at + 1, steps: x.stages.length,
      track: x.stages.map((s, i) => ({code: s.code, state: at < 0 || i > at ? 'ahead' : i === at ? 'current' : 'done'})),
      next: x.next ? `${x.next.code} · ${x.next.name}` : null,
      readiness: r ? r.score : null,
      ready: !!r && r.score === 1,
      readinessLine: !x.next ? 'No further gate' : !r ? 'Readiness not computed' :
        r.requiredTotal === 0 ? 'No required criteria' : `${r.requiredMet} of ${r.requiredTotal} required met`,
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
