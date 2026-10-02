/** The Project Portal app: portfolios on the left, the selected portfolio's own stage model as
 * a board in the middle, the selected program's gate on the right. Every portfolio draws its
 * columns from its own stages, so three portfolios show three different boards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, button, Control, divV, divH, signal, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';
import {portalDb, PortfolioRow, ProgramRow, StageModelRow, StageRow, StageTransitionRow} from '../generated/db';
import {ContextCriterion, programContext, ProgramContext, programData, readinessFor, recordDecision, setCriterion, stagesOf,
  verdictFor, enabledWorkUnits, workUnitsOf, compoundsOf} from '../portal';
import {COMPOUND_ROLES} from '../model/compounds';
import {WORK_UNIT_LABELS, WorkUnitKind} from '../model/work-units';
import {CriterionReadiness, Readiness} from '../model/readiness';
import {DataRow} from '../model/conditions';
import {Decision, nextStage} from '../model/gates';
import {openPortfolioTable} from './portfolio-table-view';
import {existingSpaces} from '../workspace';
import {CriterionRuns, criterionRuns, SavedRunItem, savedRunsOf} from '../model-runs';
import {describeSelection} from '../model/run-selection';
import {openChooseRunDialog} from './choose-run-dialog';
import {datasetsOf, openProgramTable} from '../program-data';
import {MAIN_DATASET} from '../model/program-data';
import {showModelCard} from './model-card';
import {openSaveRunDialog} from './save-run-dialog';
import {attentionFor, attentionRulesOf, openAttentionItem} from '../attention';
import {AttentionSeverity, ProgramAttention} from '../model/attention';
import {recordedSync} from '../sync-runs';
import {openSyncHistory} from './sync-history';
import {openExportDialog} from './export-dialog';
import {ProgramRecommendation, recommendFor, runModelFor} from '../recommend';
import {describeMissing} from '../model/recommend';
import {follow, myFollows, unfollow, whatsNew, WorkUnitRef} from '../new-data';
import {FollowRow} from '../generated/db';
import {ADVANCED_BOARD_BUTTONS, BoardButton, MORE_SECTIONS, PageSection, Persona, PERSONA_LABELS, PERSONAS, PORTAL_MODES,
  PersonaConfig, PortalMode, QuickActionCode, QuickCard, quickCards, sectionShown} from '../model/personas';
import {personaIn, runFunctionCard, savedMode, saveMode, savePersona} from '../personas';
import {missingAssayFrame, recommendationFrame} from '../recommend';
import {attentionFrame} from '../attention';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const statusVariant: Record<string, BadgeVariant> = {
  active: 'success', on_hold: 'warning', terminated: 'error', completed: 'accent',
  transitioned: 'accent', missing_from_feed: 'error',
};

/** Automation names are derived from codes: 'SM-DISC' → 'sm-disc', 'Program A' → 'program-a'. */
export const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

const APP_PATH = '/apps/ProjectPortalCore/Portal';

export class PortalView {
  readonly root = named(divH([], 'pp-root'), 'portal');
  private readonly portfolioList = named(divV([], 'pp-portfolios'), 'portfolio-list');
  private readonly board = named(divV([], 'pp-board'), 'stage-board');
  private readonly detail = named(divV([], 'pp-detail'), 'program-detail');
  readonly status = signal('Loading portfolios…');
  /** The app's address: the open portfolio's code is its last segment, so a board can be linked. */
  readonly path = signal(APP_PATH);

  private portfolios: PortfolioRow[] = [];
  private models = new Map<string, StageModelRow>();
  private portfolio: PortfolioRow | null = null;
  private stages: StageRow[] = [];
  private programs: ProgramRow[] = [];
  private selected: string | null = null;
  /** The stage whose criteria the gate panel shows; null = the program's next stage. */
  private viewedStage: string | null = null;
  /** Where the last run saved from this page went, by program id. */
  private savedLine = new Map<string, string>();
  /** The current user's follows of the selected program and its work units. */
  private follows: FollowRow[] = [];
  /** What landed since the user's last visit, read once per open; undefined = not read yet. */
  private news: {since: string | null, lines: string[]} | null | undefined;
  /** The selected program's dataset, read once per open; undefined = not read yet. */
  private data: {rows: DataRow[] | null, error?: string} | undefined;
  /** Simple or Advanced, remembered per user. */
  private mode: PortalMode = savedMode();
  /** The persona whose quick actions head the program page, in the open portfolio. */
  private persona: {persona: Persona, chosen: boolean} | null = null;
  /** In Simple mode, whether More is open on the program page. */
  private moreOpen = false;
  /** Draws the program page again from what was last read: a change of mode or persona reads nothing. */
  private redraw: (() => void) | null = null;
  /** Ribbon buttons Simple mode hides; the app view hands them over. */
  readonly boardButtons = new Map<BoardButton, HTMLElement>();

  constructor() {
    this.root.append(this.portfolioList, this.board, this.detail);
  }

  /** The program page; the portfolio table docks it beside its grid. */
  get detailRoot(): HTMLElement {
    return this.detail;
  }

  /** The open portfolio's code, or null before one is open. */
  get portfolioCode(): string | null {
    return this.portfolio?.code ?? null;
  }

  /** Opens [programId] of [portfolioCode], switching portfolios first when needed. */
  async showProgram(portfolioCode: string, programId: string): Promise<void> {
    if (this.portfolio?.code !== portfolioCode)
      await this.openPortfolio(portfolioCode);
    await this.openProgram(programId);
  }

  async load(portfolioCode?: string): Promise<void> {
    const [portfolios, models] = await Promise.all([portalDb.portfolios.query({}), portalDb.stageModels.query({})]);
    this.portfolios = portfolios.sort((a, b) => a.name.localeCompare(b.name));
    this.models = new Map(models.map((m) => [m.id, m]));
    const counts = await Promise.all(this.portfolios.map((p) =>
      portalDb.programs.count({filter: `portfolio_id = "${p.id}"`}).catch(() => 0)));
    this.renderPortfolios(counts);
    const first = this.portfolios.find((p) => p.code === portfolioCode) ?? this.portfolios[0];
    if (first)
      await this.openPortfolio(first.code);
    else
      this.status.value = 'No portfolios are loaded; run ProjectPortalCore:loadPortfolios';
  }

  private renderPortfolios(counts: number[]): void {
    this.portfolioList.replaceChildren(span('Portfolios', 'pp-heading'));
    this.portfolios.forEach((p, i) => {
      const model = this.models.get(p.stage_model_id);
      const card = named(divV([
        span(p.name, 'pp-portfolio-name'),
        span(`${counts[i]} ${(p.program_noun ?? 'program').toLowerCase()}s · ${model?.name ?? ''}`, 'pp-muted'),
      ], 'pp-portfolio'), `portfolio-${slug(p.code)}`);
      card.onclick = () => void this.openPortfolio(p.code);
      if (this.portfolio?.id === p.id)
        card.classList.add('pp-selected');
      this.portfolioList.append(card);
    });
  }

  async openPortfolio(code: string): Promise<void> {
    const p = this.portfolios.find((x) => x.code === code);
    if (!p)
      return;
    this.portfolio = p;
    this.selected = null;
    this.redraw = null;
    this.persona = await personaIn(p).catch(() => ({persona: 'lead' as Persona, chosen: false}));
    this.path.value = `${APP_PATH}/${p.code}`;
    for (const el of Array.from(this.portfolioList.querySelectorAll('.pp-portfolio')))
      el.classList.toggle('pp-selected', (el as HTMLElement).dataset.u2Name === `portfolio-${slug(code)}`);
    [this.stages, this.programs] = await Promise.all([
      stagesOf(p.stage_model_id),
      portalDb.programs.query({filter: `portfolio_id = "${p.id}"`}),
    ]);
    this.renderBoard();
    this.detail.replaceChildren(span(`Select a ${(p.program_noun ?? 'program').toLowerCase()}`, 'pp-muted'));
    const model = this.models.get(p.stage_model_id);
    this.status.value = `${p.name}: ${this.programs.length} ${(p.program_noun ?? 'program').toLowerCase()}s on ` +
      `${model?.name ?? '?'} ${model?.model_version ?? ''} (${this.stages.length} stages); stages owned by the ` +
      `${p.stage_authority === 'feed' ? 'registry' : 'portal'}`;
  }

  private renderBoard(): void {
    const p = this.portfolio!;
    const model = this.models.get(p.stage_model_id);
    const header = named(divH([
      span(p.name, 'pp-title'),
      badge(p.setting ?? ''),
      badge(`${model?.code ?? ''}`, {variant: 'accent'}),
      badge(p.stage_authority === 'feed' ? 'stages from registry' : 'stages decided here',
        {variant: p.stage_authority === 'feed' ? 'warning' : 'success'}),
      this.preferences(),
    ], 'pp-board-header'), 'board-header');
    const columns = divH([], 'pp-columns');
    const unstaged = this.programs.filter((x) => !x.stage_id);
    if (unstaged.length > 0)
      columns.append(this.column('unstaged', 'Unstaged', 'Not placed in the model yet', unstaged));
    for (const s of this.stages) {
      columns.append(this.column(slug(s.code), `${s.code} · ${s.name}`, s.phase ?? '',
        this.programs.filter((x) => x.stage_id === s.id)));
    }
    this.board.replaceChildren(header, columns);
  }

  private column(key: string, title: string, subtitle: string, programs: ProgramRow[]): HTMLElement {
    const cards = programs.sort((a, b) => a.name.localeCompare(b.name)).map((x) => this.card(x));
    // an empty stage is a narrow strip, so a 16-gate model keeps its populated stages in view
    const cls = programs.length === 0 ? 'pp-column pp-column-empty' : 'pp-column';
    const column = named(divV([
      divV([span(title, 'pp-column-title'), span(subtitle, 'pp-muted')], 'pp-column-header'),
      named(span(String(programs.length), 'pp-count'), `count-${key}`),
      ...cards,
    ], cls), `stage-${key}`);
    column.title = subtitle ? `${title} — ${subtitle}` : title;
    return column;
  }

  private card(x: ProgramRow): HTMLElement {
    const el = named(divV([
      span(x.name, 'pp-card-title'),
      divH([
        badge(x.status ?? 'active', {variant: statusVariant[x.status ?? 'active'] ?? 'default'}),
        ...(x.modality ? [span(x.modality, 'pp-muted')] : []),
        ...(x.predecessor_id ? [badge('continued', {variant: 'accent'})] : []),
      ], 'pp-card-meta'),
    ], 'pp-card'), `program-${slug(x.code)}`);
    if (x.id === this.selected)
      el.classList.add('pp-selected');
    el.onclick = () => void this.openProgram(x.id);
    return el;
  }

  async openProgram(id: string, viewedStage: string | null = null): Promise<void> {
    if (id !== this.selected) {
      this.data = undefined;
      this.news = undefined;
    }
    this.selected = id;
    this.viewedStage = viewedStage;
    for (const el of Array.from(this.board.querySelectorAll('.pp-card')))
      el.classList.remove('pp-selected');
    const program = this.programs.find((x) => x.id === id);
    if (program)
      this.board.querySelector(`[data-u2-name="program-${slug(program.code)}"]`)?.classList.add('pp-selected');
    const ctx = await programContext(id);
    const transitions = (await portalDb.stageTransitions.query({filter: `program_id = "${id}"`}))
      .sort((a, b) => String(b.created_on).localeCompare(String(a.created_on)));
    const successors = await portalDb.programs.query({filter: `predecessor_id = "${id}"`}).catch(() => []);
    const readiness = await this.readinessOf(ctx, viewedStage);
    const units = await Promise.all(enabledWorkUnits(ctx.portfolio).map(async (kind) =>
      ({kind, rows: await workUnitsOf(kind, [id]).catch(() => [])})));
    const compounds = await compoundsOf(id).catch(() => []);
    const spaceId = ctx.program.space_id;
    const spaceLive = !!spaceId && (await existingSpaces([spaceId]).catch(() => new Set<string>())).has(spaceId);
    const runs = await savedRunsOf(id).catch(() => []);
    // criteria of the shown stage that rest on a run: their candidates and choices
    const choices = new Map<string, CriterionRuns>();
    for (const c of readiness?.criteria ?? []) {
      if (ctx.criteria.find((x) => x.id === c.id)?.artifact === 'model_run') {
        const r = await criterionRuns(id, c.id).catch(() => null);
        if (r)
          choices.set(c.id, r);
      }
    }
    const attention = await this.attentionOf(ctx);
    const recommendation = await recommendFor(ctx, {readiness, rows: this.data?.rows, dataError: this.data?.error})
      .catch((e) => {
        grok.log.warning(`Project Portal: no model recommendation for ${ctx.program.code}: ${e?.message ?? e}`);
        return null;
      });
    // what landed since the last visit, for a user who follows the program or one of its series
    const follows = await myFollows(id).catch(() => [] as FollowRow[]);
    // read once per open of a program: a click on the stage track is not a new visit
    if (this.news === undefined && follows.length) {
      if (this.data === undefined) {
        try {
          this.data = {rows: await programData(ctx.portfolio, ctx.program.code)};
        } catch (e: any) {
          this.data = {rows: null, error: e?.message ?? String(e)};
        }
      }
      this.news = await whatsNew(ctx.portfolio, ctx.program, this.data.rows, follows).catch(() => null);
    }
    const news = follows.length ? this.news ?? null : null;
    this.follows = follows;
    this.redraw = () => this.renderDetail(ctx, transitions, successors, readiness, units, compounds, spaceLive, runs,
      choices, attention, recommendation, news);
    this.redraw();
  }

  /** Mode and persona, on the board header. */
  private preferences(): HTMLElement {
    const mode = this.select('mode', 'Mode', PORTAL_MODES.map((m) => ({code: m, name: m === 'simple' ? 'Simple' : 'Advanced'})),
      this.mode);
    mode.title = 'Simple keeps analytics, runs and history behind More; remembered for you';
    mode.onchange = () => this.setMode(mode.value as PortalMode);
    const persona = this.select('persona', 'Persona', PERSONAS.map((p) => ({code: p, name: PERSONA_LABELS[p]})),
      this.persona?.persona ?? '');
    persona.title = this.persona?.chosen ? 'Your persona: it picks the quick actions; remembered for you' :
      'Picked from your groups in this portfolio; choose one to keep it';
    persona.onchange = () => this.setPersona(persona.value as Persona);
    return named(divH([mode.parentElement!, persona.parentElement!], 'pp-preferences'), 'preferences');
  }

  /** A labelled native select; its wrapper carries the name, as in the dialogs. */
  private select(name: string, label: string, options: {code: string, name: string}[], value: string): HTMLSelectElement {
    const el = document.createElement('select');
    el.className = 'pp-select pp-select-small';
    el.append(...options.map((o) => Object.assign(document.createElement('option'), {value: o.code, textContent: o.name})));
    el.value = value;
    named(divH([span(label, 'pp-muted'), el], 'pp-preference'), `${name}-input`);
    return el;
  }

  get currentMode(): PortalMode {
    return this.mode;
  }

  setMode(mode: PortalMode): void {
    this.mode = mode;
    this.moreOpen = false;
    saveMode(mode);
    this.applyBoardButtons();
    this.redraw?.();
  }

  setPersona(persona: Persona): void {
    this.persona = {persona, chosen: true};
    savePersona(persona);
    this.redraw?.();
  }

  /** Simple mode leaves the keepers' buttons (sync history, export, all portfolios) off the ribbon. */
  applyBoardButtons(): void {
    for (const [key, el] of this.boardButtons)
      el.style.display = this.mode === 'simple' && ADVANCED_BOARD_BUTTONS.includes(key) ? 'none' : '';
  }

  /** The portfolio's attention rules over the program's data, read once per open; null without rules. */
  private async attentionOf(ctx: ProgramContext): Promise<{result: ProgramAttention | null, error?: string} | null> {
    if (attentionRulesOf(ctx.portfolio).length === 0)
      return null;
    if (this.data === undefined) {
      try {
        this.data = {rows: await programData(ctx.portfolio, ctx.program.code)};
      } catch (e: any) {
        this.data = {rows: null, error: e?.message ?? String(e)};
      }
    }
    return this.data.error ? {result: null, error: this.data.error} : {result: attentionFor(ctx.portfolio, this.data.rows)};
  }

  /** Readiness for the viewed stage: the context already holds the next stage's. */
  private async readinessOf(ctx: ProgramContext, stageId: string | null): Promise<Readiness | null> {
    if (stageId === null || stageId === ctx.readiness?.stageId)
      return ctx.readiness;
    if (ctx.criteria.some((c) => c.stageId === stageId && c.condition) && this.data === undefined) {
      try {
        this.data = {rows: await programData(ctx.portfolio, ctx.program.code)};
      } catch (e: any) {
        this.data = {rows: null, error: e?.message ?? String(e)};
      }
    }
    return readinessFor(ctx, stageId, this.data?.rows ?? null, this.data?.error);
  }

  private renderDetail(ctx: ProgramContext, transitions: StageTransitionRow[], successors: ProgramRow[],
    readiness: Readiness | null, units: {kind: WorkUnitKind, rows: Record<string, any>[]}[] = [],
    compounds: {link: any, compound: any}[] = [], spaceLive = false, runs: SavedRunItem[] = [],
    choices = new Map<string, CriterionRuns>(), attention: {result: ProgramAttention | null, error?: string} | null = null,
    recommendation: ProgramRecommendation | null = null, news: {since: string | null, lines: string[]} | null = null): void {
    const {program, stages} = ctx;
    const current = stages.find((s) => s.id === program.stage_id) ?? null;
    const next = nextStage(stages.map((s) => ({...s, isTerminal: s.is_terminal})), program.stage_id ?? null);
    const stageName = (id?: string | null) => stages.find((s) => s.id === id)?.name ?? '—';

    const viewed = stages.find((s) => s.id === (readiness?.stageId ?? null)) ?? null;
    const track = named(divH(stages.map((s) => {
      const cls = s.id === program.stage_id ? 'pp-step pp-step-current' :
        current && s.ordinal < current.ordinal ? 'pp-step pp-step-done' : 'pp-step';
      const step = named(span(s.code, cls + (s.id === viewed?.id ? ' pp-step-viewed' : '')), `step-${slug(s.code)}`);
      step.title = `${s.name}: show what entering it takes`;
      step.onclick = () => void this.openProgram(program.id, s.id === next?.id ? null : s.id);
      return step;
    }), 'pp-track'), 'stage-track');

    const editable = ctx.approver && ctx.portfolio.stage_authority === 'portal';
    const criteriaList = named(divV((readiness?.criteria ?? []).map((c) => this.criterionRow(program.id, c, editable,
      choices.get(c.id) ? {runs: choices.get(c.id)!, criterion: ctx.criteria.find((x) => x.id === c.id)!} : null)),
    'pp-criteria'), 'gate-criteria');
    const score = readiness ? named(divH([
      span(`${readiness.requiredMet} of ${readiness.requiredTotal} required met or waived`),
      badge(`${Math.round(readiness.score * 100)}%`,
        {variant: readiness.score === 1 ? 'success' : readiness.score === 0 ? 'error' : 'warning'}),
    ], 'pp-readiness'), 'readiness-score') : null;

    const rationale = named(document.createElement('textarea'), 'rationale-input');
    rationale.placeholder = 'Rationale';
    rationale.rows = 2;
    const verdict = verdictFor(ctx, 'advance');
    const reasons = named(divV(verdict.reasons.map((r) => span(r, 'pp-reason'))), 'gate-reasons');
    const decide = (d: Decision) => async () => {
      try {
        const done = await recordDecision(program.id, d, null, rationale.value);
        grok.shell.info(`${program.name}: ${d} recorded${done.workspace ? `\n${done.workspace}` : ''}`);
        await this.openPortfolio(this.portfolio!.code);
        await this.openProgram(program.id);
      } catch (e: any) {
        grok.shell.error(e?.message ?? String(e));
      }
    };
    const advance = named(button(next ? `Advance to ${next.name}` : 'Advance', decide('advance'), {primary: true}),
      'advance-button');
    advance.disabled = !verdict.allowed;
    const hold = named(button('Hold', decide('hold')), 'hold-button');
    const terminate = named(button('Terminate', decide('terminate')), 'terminate-button');
    hold.disabled = terminate.disabled = !ctx.approver || ctx.portfolio.stage_authority === 'feed';

    const history = named(divV(transitions.map((t) => {
      const row = span(
        `${(t.decided_on ?? t.created_on)?.format?.('YYYY-MM-DD') ?? ''} · ${t.decision}: ${stageName(t.from_stage_id)} → ` +
        `${stageName(t.to_stage_id)}${t.source === 'data' ? ' (by data)' : ''}${t.rationale ? ` — ${t.rationale}` : ''}`,
        'pp-history-row');
      // what the data said when it advanced the program
      if (t.evidence)
        row.title = t.evidence;
      return row;
    })), 'transition-history');

    const lineage: HTMLElement[] = [];
    if (program.predecessor_id)
      lineage.push(named(span(`Continues a program of another portfolio`, 'pp-lineage'), 'lineage-predecessor'));
    for (const s of successors)
      lineage.push(named(span(`Continued as ${s.code} in another portfolio`, 'pp-lineage'), 'lineage-successor'));

    const workspace = spaceLive ?
      named(button('Open workspace', () => void this.openWorkspace(program.space_id!)), 'open-workspace-button') :
      named(span(program.space_id ? 'The workspace folder is gone: rebuild it with buildWorkspace' :
        'No workspace yet: build it with buildWorkspace', 'pp-muted'), 'no-workspace');

    // the program's governed data, one button per dataset the portfolio declares
    const data = datasetsOf(ctx.portfolio).map((d) => {
      const b = named(button(d.code === MAIN_DATASET ? 'Open data' : `Open ${d.name}`, () =>
        void openProgramTable(ctx.portfolio.code, program.code, d.code)
          .catch((e) => grok.shell.error(e?.message ?? String(e)))),
      d.code === MAIN_DATASET ? 'open-data-button' : `open-dataset-${slug(d.code)}`);
      b.title = d.code === MAIN_DATASET ? `${d.name}: this ${(ctx.portfolio.program_noun ?? 'program').toLowerCase()}'s rows as a table` :
        d.description ?? `${d.name} of this ${(ctx.portfolio.program_noun ?? 'program').toLowerCase()} as a table`;
      return b;
    });

    const followed = this.follows.some((f) => !f.work_unit_code);
    const followButton = named(button(followed ? 'Unfollow' : 'Follow', () => void this.toggleFollow(ctx, null)),
      'follow-button');
    followButton.title = followed ? 'Stop being notified of new data in this program' :
      'Be notified, in the platform\'s notifications, when new results land in this program\'s data';
    const newsLine = news ? named(divV(news.lines.length ?
      [span(`What's new since your last visit${news.since ? ` (${news.since.slice(0, 16).replace('T', ' ')} UTC)` : ''}`, 'pp-muted'),
        ...news.lines.map((l) => span(l))] :
      [span(`Nothing new since your last visit${news.since ? ` (${news.since.slice(0, 16).replace('T', ' ')} UTC)` : ''}`, 'pp-muted')],
    'pp-units'), 'whats-new') : null;

    const openSave = () => void openSaveRunDialog({
      context: {portfolio: ctx.portfolio.code, program: program.code},
      onSaved: async (_r, line) => {
        this.savedLine.set(program.id, line);
        if (this.selected === program.id)
          await this.openProgram(program.id, this.viewedStage);
      },
    }).catch((e) => grok.shell.error(e?.message ?? String(e)));
    const saveRun = named(button('Save a run…', openSave), 'save-a-run-button');
    saveRun.title = 'Save a run of a catalogued model, made in this session, to this program';

    const gateHeading = named(span(!viewed ? 'No further stage' : viewed.id === next?.id ? `Gate to ${viewed.name}` :
      `What entering ${viewed.name} takes`, 'pp-heading'), 'gate-heading');
    const quick = this.quickActions(ctx, {next, spaceLive, choices, attention, recommendation, followed,
      openSave, gate: () => {
        gateHeading.scrollIntoView({block: 'start'});
        rationale.focus();
      }});

    // Simple mode keeps the heavier sections behind More
    const shown = (x: PageSection) => sectionShown(this.mode, x, this.moreOpen);
    const more = this.mode === 'simple' ? named(button(this.moreOpen ? 'Less' : `More (${MORE_SECTIONS.length})`, () => {
      this.moreOpen = !this.moreOpen;
      this.redraw?.();
    }), 'more-button') : null;
    if (more)
      more.title = this.moreOpen ? 'Back to the simple page' : 'Missing assays and suggested models, attention, saved runs, history';

    this.detail.replaceChildren(
      named(span(program.name, 'pp-title'), 'program-title'),
      divH([badge(program.status ?? 'active', {variant: statusVariant[program.status ?? 'active'] ?? 'default'}),
        span(`${ctx.portfolio.name} · ${program.code}`, 'pp-muted')]),
      named(divH([workspace, followButton], 'pp-actions'), 'program-actions'),
      ...(newsLine ? [newsLine] : []),
      ...quick,
      ...(data.length ? [named(divH(data, 'pp-actions pp-data-actions'), 'program-data-actions')] : []),
      ...(program.source_stage ? [span(`Registry stage: ${program.source_stage}`, 'pp-muted')] : []),
      ...lineage,
      track,
      gateHeading,
      ...(score && readiness!.criteria.length ? [score] : []),
      readiness?.criteria.length ? criteriaList : span('No criteria declared', 'pp-muted'),
      ...(readiness?.dataError ? [named(span(`Program data could not be read: ${readiness.dataError}`, 'pp-reason'),
        'data-error')] : []),
      reasons,
      rationale,
      divH([advance, hold, terminate], 'pp-actions'),
      ...(units.length ? [span('Work units', 'pp-heading'), this.workUnits(units, ctx)] : []),
      ...(compounds.length ? [span('Compounds', 'pp-heading'), this.compounds(compounds)] : []),
      ...(more ? [divH([more], 'pp-actions')] : []),
      ...(recommendation && shown('recommendation') ? this.recommendation(recommendation) : []),
      ...(attention && shown('attention') ? this.attention(ctx.portfolio.code, program.code, attention) : []),
      ...(shown('saved_runs') ? [
        named(span(`Saved runs (${runs.length})`, 'pp-heading'), 'saved-runs-heading'),
        divH([saveRun], 'pp-actions'),
        ...(this.savedLine.has(program.id) ? [named(span(this.savedLine.get(program.id)!, 'pp-muted'), 'save-run-result')] : []),
        runs.length ? this.savedRuns(runs) : named(span('No runs saved to this program yet', 'pp-muted'), 'no-saved-runs'),
      ] : []),
      ...(shown('history') ? [
        span('History', 'pp-heading'),
        transitions.length ? history : span('No transitions yet', 'pp-muted'),
      ] : []),
    );
  }

  /** The persona's quick-action cards for this program: each opens something the page already
   * does. A card that cannot be used now says why on hover. */
  private quickActions(ctx: ProgramContext, s: {next: {id: string, name: string} | null, spaceLive: boolean,
    choices: Map<string, CriterionRuns>, attention: {result: ProgramAttention | null, error?: string} | null,
    recommendation: ProgramRecommendation | null, followed: boolean, openSave: () => void, gate: () => void}): HTMLElement[] {
    if (!this.persona)
      return [];
    const {program, portfolio} = ctx;
    const noun = (portfolio.program_noun ?? 'program').toLowerCase();
    const datasets = datasetsOf(portfolio);
    const main = datasets.find((d) => d.code === MAIN_DATASET);
    const other = datasets.find((d) => d.code !== MAIN_DATASET);
    const run = [...s.choices.entries()][0];
    const r = s.recommendation;
    const a = s.attention;
    const unavailable: Partial<Record<QuickActionCode, string>> = {};
    if (!main)
      unavailable.open_data = 'The portfolio declares no program dataset';
    if (!other)
      unavailable.open_other_data = 'The portfolio declares no further dataset';
    if (!a)
      unavailable.needs_attention = 'The portfolio declares no attention rules';
    else if (!a.result)
      unavailable.needs_attention = `Program data could not be read: ${a.error}`;
    else if (a.result.count === 0)
      unavailable.needs_attention = 'Nothing needs attention';
    if (!r || !r.missing.length)
      unavailable.missing_assays = `Nothing missing for ${s.next?.name ?? 'the gate'}`;
    if (!r || !r.models.length)
      unavailable.suggested_models = 'No catalogued model fits this program';
    if (!run)
      unavailable.choose_run = 'No criterion of this gate rests on a run';
    else if (run[1].refused.length)
      unavailable.choose_run = run[1].refused.join('; ');
    if (!s.spaceLive)
      unavailable.open_workspace = 'No workspace yet';
    if (!s.next)
      unavailable.gate_decision = 'No further stage';
    const labels: Partial<Record<QuickActionCode, {label?: string, description?: string}>> = {
      open_data: {label: main ? `Open ${main.name}` : undefined, description: `This ${noun}'s rows as a table`},
      open_other_data: other ? {label: `Open ${other.name}`, description: other.description ?? `${other.name} of this ${noun}`} : {},
      needs_attention: a?.result ? {label: `Needs attention (${a.result.count})`} : {},
      missing_assays: r?.missing.length ? {label: `Missing assays (${r.missing.length})`} : {},
      suggested_models: r?.models.length ? {label: `Suggested models (${r.models.length})`,
        description: r.models.slice(0, 2).map((m) => m.name).join(', ')} : {},
      follow: s.followed ? {label: 'Unfollow', description: 'Stop hearing about new results'} : {},
      gate_decision: s.next ? {label: `Gate to ${s.next.name}`, description: ctx.approver && portfolio.stage_authority === 'portal' ?
        'Readiness, then advance, hold or terminate' : 'Readiness; the approvers decide'} : {},
      choose_run: run ? {description: ctx.criteria.find((c) => c.id === run[0])?.name} : {},
    };
    const failed = (e: any) => grok.shell.error(e?.message ?? String(e));
    const act: Record<QuickActionCode, () => void> = {
      open_data: () => void openProgramTable(portfolio.code, program.code, MAIN_DATASET).catch(failed),
      open_other_data: () => void openProgramTable(portfolio.code, program.code, other!.code).catch(failed),
      needs_attention: () => grok.shell.addTableView(attentionFrame(a!.result!, attentionRulesOf(portfolio),
        `${program.code} · Needs attention`)),
      missing_assays: () => grok.shell.addTableView(missingAssayFrame(r!)),
      suggested_models: () => grok.shell.addTableView(recommendationFrame(r!)),
      save_run: s.openSave,
      choose_run: () => {
        const c = ctx.criteria.find((x) => x.id === run![0])!;
        void openChooseRunDialog({programId: program.id, criterion: {id: c.id, code: c.code, name: c.name,
          runFilter: c.runFilter}, runs: run![1], onSelected: () => this.openProgram(program.id, this.viewedStage)}).catch(failed);
      },
      follow: () => void this.toggleFollow(ctx, null),
      open_workspace: () => void this.openWorkspace(program.space_id!),
      gate_decision: s.gate,
      portfolio_table: () => void openPortfolioTable(portfolio.code),
    };
    const cards = quickCards(this.persona.persona, portfolio.personas as PersonaConfig | undefined, unavailable, labels);
    const stage = ctx.stages.find((x) => x.id === program.stage_id)?.code ?? null;
    const card = (c: QuickCard) => {
      const el = named(divV([span(c.label, 'pp-quick-label'), span(c.description, 'pp-muted')],
        c.unavailable ? 'pp-quick pp-quick-off' : 'pp-quick'), `quick-${slug(c.key)}`);
      el.title = c.unavailable ?? c.description;
      if (!c.unavailable) {
        el.onclick = c.action ? act[c.action] : () => void runFunctionCard(c.function!,
          {portfolio: portfolio.code, program: program.code, stage}).catch(failed);
      }
      return el;
    };
    return [
      named(span(`Quick actions · ${PERSONA_LABELS[this.persona.persona]}`, 'pp-heading'), 'quick-actions-heading'),
      named(divH(cards.map(card), 'pp-quick-cards'), 'quick-actions'),
    ];
  }

  /** Missing assays for the shown gate (what would close it), then the catalogued models that fit
   * the program, best first. A model that predicts a missing assay says so; Run opens its dialog
   * filled in from the program. */
  private recommendation(r: ProgramRecommendation): HTMLElement[] {
    const out: HTMLElement[] = [];
    const modelLink = (code: string, name: string, prefix = 'recommended-card') => {
      const link = named(span(name, 'pp-link'), `${prefix}-${slug(code)}`);
      link.title = 'Open the model card';
      link.onclick = () => void showModelCard(code);
      return link;
    };
    if (r.missing.length) {
      out.push(named(span(`Missing assays${r.stage ? ` for ${r.stage.name}` : ''} (${r.missing.length})`, 'pp-heading'),
        'missing-assays-heading'));
      out.push(named(divV(r.missing.map((m) => {
        const by = m.column ? r.predictors.get(m.column) ?? [] : [];
        return named(divV([
          span(describeMissing(m)),
          ...(by.length ? [divH([span('Predicted meanwhile by', 'pp-muted'), ...by.map((x) => modelLink(x.code, x.name, 'predictor-card'))],
            'pp-run-meta')] : []),
        ], 'pp-run'), `missing-assay-${slug(m.criterion)}${m.column ? `-${slug(m.column)}` : ''}`);
      }), 'pp-runs'), 'missing-assays'));
    }
    out.push(named(span(`Suggested models (${r.models.length})`, 'pp-heading'), 'suggested-models-heading'));
    if (!r.models.length) {
      out.push(named(span('No catalogued model fits this program\'s stage, modality and work units', 'pp-muted'),
        'no-suggested-models'));
      return out;
    }
    out.push(named(divV(r.models.map((m) => {
      const run = named(button('Run', () => void runModelFor(m.code, {portfolio: r.portfolio.code, program: r.program.code,
        series: r.series}).catch((e) => grok.shell.error(e?.message ?? String(e)))), `run-model-${slug(m.code)}`);
      run.disabled = !m.installed;
      run.title = m.installed ? 'Open the model\'s dialog, filled in from this program' : 'The model\'s function is not installed here';
      return named(divV([
        divH([modelLink(m.code, m.name), ...(m.approval ? [badge(m.approval.replace('_', ' '),
          {variant: m.approval === 'approved' ? 'success' : 'default'})] : []), run], 'pp-run-meta'),
        ...(m.reasons.length ? [span(m.reasons.join(' · '), 'pp-muted')] : []),
      ], 'pp-run'), `suggested-model-${slug(m.code)}`);
    }), 'pp-runs'), 'suggested-models'));
    return out;
  }

  /** Needs attention: a count per rule that trips, then each flagged compound (or lot) with the
   * rules it trips and why. A click opens the program's data with that row current. */
  private attention(portfolioCode: string, programCode: string,
    a: {result: ProgramAttention | null, error?: string}): HTMLElement[] {
    const heading = named(span(`Needs attention${a.result ? ` (${a.result.count})` : ''}`, 'pp-heading'), 'attention-heading');
    if (!a.result)
      return [heading, named(span(`Program data could not be read: ${a.error}`, 'pp-reason'), 'attention-error')];
    const r = a.result;
    if (r.count === 0)
      return [heading, named(span('Nothing needs attention', 'pp-muted'), 'attention-none')];
    const variant = (s: AttentionSeverity): BadgeVariant => s === 'critical' ? 'error' : s === 'warning' ? 'warning' : 'default';
    const nameOf = new Map(r.rules.map((x) => [x.code, x.name]));
    const flagged = r.items.map((x) => x.item);
    const counts = named(divH(r.rules.filter((x) => x.items.length).map((x) =>
      named(badge(`${x.items.length} ${x.name}`, {variant: variant(x.severity)}), `attention-rule-${slug(x.code)}`)),
    'pp-card-meta'), 'attention-rules');
    const items = r.items.map((it) => {
      const link = named(span(it.item, 'pp-link'), `attention-open-${slug(it.item)}`);
      link.title = 'Open the program data at this row';
      link.onclick = () => void openAttentionItem(portfolioCode, programCode, it.item, flagged)
        .catch((e) => grok.shell.error(e?.message ?? String(e)));
      return named(divV([
        divH([link, ...it.rules.map((code) => badge(nameOf.get(code) ?? code,
          {variant: variant(r.rules.find((x) => x.code === code)!.severity)}))], 'pp-unit'),
        ...it.messages.map((m) => span(m, 'pp-muted')),
      ], 'pp-attention-item'), `attention-item-${slug(it.item)}`);
    });
    return [heading, named(divV([counts, ...items], 'pp-units'), 'attention')];
  }

  /** The program's campaigns, series, studies and process stages, as the portfolio enables them. */
  private workUnits(units: {kind: WorkUnitKind, rows: Record<string, any>[]}[], ctx?: ProgramContext): HTMLElement {
    // a series has rows of its own in the program dataset, so it can be followed on its own
    const followSeries = (code: string): HTMLElement[] => {
      if (!ctx)
        return [];
      const on = this.follows.some((f) => f.work_unit_kind === 'series' && f.work_unit_code === code);
      const b = named(span(on ? 'Unfollow' : 'Follow', 'pp-link'), `follow-series-${slug(code)}`);
      b.title = on ? 'Stop being notified of new data in this series' : 'Be notified when new results land in this series';
      b.onclick = () => void this.toggleFollow(ctx, {kind: 'series', code});
      return [b];
    };
    const unitVariant = (s?: string): BadgeVariant => ['active', 'running', 'in_development'].includes(s ?? '') ?
      'success' : ['closed', 'completed', 'locked', 'transferred'].includes(s ?? '') ? 'accent' :
        ['dropped', 'cancelled', 'paused', 'deprioritized'].includes(s ?? '') ? 'warning' : 'default';
    return named(divV(units.map(({kind, rows}) => {
      const sorted = [...rows].sort((a, b) => (a.process_type ?? '').localeCompare(b.process_type ?? '') ||
        (a.ordinal ?? 0) - (b.ordinal ?? 0) || a.code.localeCompare(b.code));
      return named(divV([
        span(`${WORK_UNIT_LABELS[kind]} (${rows.length})`, 'pp-unit-kind'),
        ...(rows.length === 0 ? [span('None yet', 'pp-muted')] : sorted.map((u) => named(divH([
          span(u.name === u.code ? u.code : `${u.code} · ${u.name}`),
          ...(u.status ? [badge(String(u.status).replace('_', ' '), {variant: unitVariant(u.status)})] : []),
          ...(u.maturity ? [badge(u.maturity, {variant: 'accent'})] : []),
          ...(u.outcome && u.outcome !== 'open' ? [badge(u.outcome, {variant: u.outcome === 'succeeded' ? 'success' : 'warning'})] : []),
          ...(u.type ? [span(String(u.type).replace('_', ' '), 'pp-muted')] : []),
          ...(u.process_type ? [span(u.process_type, 'pp-muted')] : []),
          ...(u.missing_from_feed ? [badge('missing from feed', {variant: 'error'})] : []),
          ...(kind === 'series' ? followSeries(u.code) : []),
        ], 'pp-unit'), `unit-${slug(kind)}-${slug(u.code)}`))),
      ], 'pp-unit-group'), `units-${slug(kind)}`);
    }), 'pp-units'), 'work-units');
  }

  /** Compounds by role: a count per role, and the ones that matter by name (lead and beyond).
   * Rows, never folders (FR-26). */
  private compounds(items: {link: any, compound: any}[]): HTMLElement {
    const counts = COMPOUND_ROLES.map((r) => [r, items.filter((x) => x.link.role === r).length] as const).filter(([, n]) => n > 0);
    const shown = items.filter((x) => ['candidate', 'lead', 'backup', 'tool'].includes(x.link.role))
      .sort((a, b) => COMPOUND_ROLES.indexOf(b.link.role) - COMPOUND_ROLES.indexOf(a.link.role) ||
        a.compound.code.localeCompare(b.compound.code)).slice(0, 12);
    return named(divV([
      named(divH(counts.map(([r, n]) => badge(`${n} ${r}`, {variant: r === 'candidate' || r === 'lead' ? 'success' :
        r === 'discontinued' ? 'warning' : 'default'})), 'pp-card-meta'), 'compound-roles'),
      ...shown.map((x) => named(divH([span(x.compound.code), badge(x.link.role, {variant: 'accent'}),
        ...(x.compound.name && x.compound.name !== x.compound.code ? [span(x.compound.name, 'pp-muted')] : []),
        ...(x.link.missing_from_feed ? [badge('missing from feed', {variant: 'error'})] : []),
      ], 'pp-unit'), `compound-${slug(x.compound.code)}`)),
    ], 'pp-units'), 'program-compounds');
  }

  /** The model runs saved to the program, newest first: title, then the model (its card opens on
   * a click), the stage at run, unit operation, run type, purpose, author's date. */
  private savedRuns(items: SavedRunItem[]): HTMLElement {
    return named(divV(items.map(({run, model, stage, unitOperation, runType}) => {
      const modelLink = named(span(model?.name ?? run.model_name, model ? 'pp-link' : 'pp-muted'),
        `model-link-${slug(run.model_name)}`);
      if (model) {
        modelLink.title = 'Open the model card';
        modelLink.onclick = () => void showModelCard(model.code);
      } else
        modelLink.title = 'Not in the model catalogue';
      return named(divV([
        span(run.title ?? run.run_id ?? run.model_name),
        divH([
          modelLink,
          ...(stage ? [badge(stage, {variant: 'accent'})] : []),
          ...(unitOperation ? [span(unitOperation, 'pp-muted')] : []),
          ...(runType ? [span(runType, 'pp-muted')] : []),
          ...(run.purpose ? [badge(run.purpose)] : []),
          ...(run.is_reference ? [badge('reference', {variant: 'success'})] : []),
          span((run.run_started ?? run.created_on)?.format?.('YYYY-MM-DD') ?? '', 'pp-muted'),
        ], 'pp-run-meta'),
      ], 'pp-run'), `saved-run-${slug(run.title ?? run.run_id ?? run.id)}`);
    }), 'pp-runs'), 'saved-runs');
  }

  /** One criterion of the gate: its state, and for a data condition what the data says. A computed
   * condition can only be waived; "Met" is for what a person records. */
  private criterionRow(programId: string, c: CriterionReadiness, editable: boolean,
    run: {runs: CriterionRuns, criterion: ContextCriterion} | null = null): HTMLElement {
    const variant: BadgeVariant = c.state === 'met' ? 'success' : c.state === 'waived' ? 'accent' :
      c.state === 'not_met' ? 'error' : c.state === 'in_progress' ? 'warning' : 'default';
    const row = named(divH([
      badge(c.state.replace('_', ' '), {variant}),
      span(`${c.name}${c.required ? '' : ' (optional)'}`),
      span(c.basis === 'evaluated' ? 'from data' : c.kind.replace('_', ' '), 'pp-muted'),
    ], 'pp-criterion'), `criterion-${slug(c.code)}`);
    if (editable) {
      if (c.basis === 'recorded' && c.state !== 'met' && c.state !== 'waived')
        row.append(named(button('Met', () => void this.mark(programId, c.id, 'met')), `mark-met-${slug(c.code)}`));
      if (c.state !== 'waived' && c.state !== 'met')
        row.append(named(button('Waive', () => void this.mark(programId, c.id, 'waived')), `mark-waived-${slug(c.code)}`));
    }
    const runLines = run ? [this.criterionRun(programId, c, run.runs, run.criterion)] : [];
    if (!c.result && !c.gap && !runLines.length)
      return row;
    const lines = c.result ? [span(c.result.summary, 'pp-muted')] : c.gap ? [span(c.gap, 'pp-muted')] : [];
    if (c.result && !c.result.met && c.result.gaps.length)
      lines.push(span(c.gap.slice(c.result.summary.length + 2), 'pp-gap'));
    return divV([row, named(divV([...lines, ...runLines], 'pp-criterion-detail'), `criterion-detail-${slug(c.code)}`)],
      'pp-criterion-block');
  }

  /** A criterion that rests on a run: the run chosen, the ones it superseded, and Choose run. */
  private criterionRun(programId: string, c: CriterionReadiness, runs: CriterionRuns, criterion: ContextCriterion): HTMLElement {
    const choose = named(button('Choose run…', () => void openChooseRunDialog({programId,
      criterion: {id: c.id, code: c.code, name: c.name, runFilter: criterion.runFilter}, runs,
      onSelected: () => this.openProgram(programId, this.viewedStage)})
      .catch((e) => grok.shell.error(e?.message ?? String(e)))), `choose-run-${slug(c.code)}`);
    choose.disabled = runs.refused.length > 0;
    choose.title = runs.refused.join('; ') || 'Compare the program\'s runs for this criterion and select one';
    return named(divV([
      named(divH([
        runs.selected ? span(`Run: ${describeSelection(runs.selected)}`) :
          span(`No run selected · ${runs.candidates.length} candidate${runs.candidates.length === 1 ? '' : 's'}`, 'pp-muted'),
        choose,
      ], 'pp-criterion-run'), `selected-run-${slug(c.code)}`),
      ...(runs.superseded.length ? [named(divV(runs.superseded.map((e) =>
        span(`Superseded: ${describeSelection(e)}`, 'pp-gap'))), `run-history-${slug(c.code)}`)] : []),
    ]), `criterion-run-${slug(c.code)}`);
  }

  /** Follows or unfollows the program, or one of its series, then redraws the page. */
  private async toggleFollow(ctx: ProgramContext, unit: WorkUnitRef | null): Promise<void> {
    try {
      const on = this.follows.some((f) => unit ? f.work_unit_kind === unit.kind && f.work_unit_code === unit.code : !f.work_unit_code);
      if (on)
        await unfollow(ctx.program, unit);
      else
        await follow(ctx.portfolio, ctx.program, unit);
      this.news = undefined;
      grok.shell.info(`${unit ? `${ctx.program.code} ${unit.code}` : ctx.program.name}: ${on ? 'no longer followed' :
        'followed; new results will show in your notifications'}`);
      await this.openProgram(ctx.program.id, this.viewedStage);
    } catch (e: any) {
      grok.shell.error(e?.message ?? String(e));
    }
  }

  /** Opens the program's folder the way the Browse tree does. A folder deleted by hand is reported,
   * not recreated: the builder is the one place that creates folders. */
  private async openWorkspace(spaceId: string): Promise<void> {
    const space = await grok.dapi.spaces.find(spaceId).catch(() => null);
    if (!space) {
      grok.shell.warning('This program\'s workspace folder is gone; rebuild the portfolio\'s workspace');
      return;
    }
    grok.shell.addView(DG.SpaceView.forProject(space));
  }

  private async mark(programId: string, criterionId: string, state: 'met' | 'waived'): Promise<void> {
    await setCriterion(programId, criterionId, state);
    await this.openProgram(programId, this.viewedStage);
  }

  /** The open portfolio as a portfolio file. */
  async exportPortfolio(): Promise<void> {
    if (this.portfolio)
      await openExportDialog(this.portfolio.code).catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  /** The open portfolio's schedule and sync runs. */
  async syncHistory(): Promise<void> {
    if (this.portfolio)
      await openSyncHistory(this.portfolio).catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  async sync(): Promise<void> {
    if (!this.portfolio)
      return;
    const code = this.portfolio.code;
    const pi = DG.TaskBarProgressIndicator.create(`Syncing ${this.portfolio.name}`);
    try {
      // recorded in the portfolio's sync history, like a scheduled run
      const r = await recordedSync(code, 'manual');
      if (r.error)
        throw new Error(r.error);
      grok.shell.info(r.report + (r.historyError ? `\n${r.historyError}` : ''));
      await this.load(code);
    } catch (e: any) {
      grok.shell.error(e?.message ?? String(e));
    } finally {
      pi.close();
    }
  }
}

export async function portalAppView(portfolioCode?: string): Promise<DG.ViewBase> {
  const portal = new PortalView();
  const content = new Control(portal.root);
  const sync = named(button('Sync', () => void portal.sync()), 'sync-button');
  sync.title = 'Run this portfolio\'s registry feeds';
  const history = named(button('Sync history', () => void portal.syncHistory()), 'sync-history-button');
  history.title = 'This portfolio\'s schedule, and what each sync changed';
  const exportButton = named(button('Export', () => void portal.exportPortfolio()), 'export-button');
  exportButton.title = 'This portfolio as a portfolio file, to diff or load on another instance';
  const table = named(button('Table', () => void openPortfolioTable(portal.portfolioCode)), 'table-button');
  table.title = 'This portfolio as a table, with filters and analytics';
  const all = named(button('All portfolios', () => void openPortfolioTable(null)), 'all-portfolios-button');
  all.title = 'Every portfolio you can see as one table';
  const view = appView({name: 'Project Portal', content, ribbon: [[sync, history, exportButton], [table, all]], status: portal.status,
    path: portal.path});
  portal.boardButtons.set('sync_history', history);
  portal.boardButtons.set('export', exportButton);
  portal.boardButtons.set('all_portfolios', all);
  portal.applyBoardButtons();
  await portal.load(portfolioCode);
  (view as any).portal = portal;
  return view;
}
