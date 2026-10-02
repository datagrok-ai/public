/** The Project Portal app: portfolios on the left, the selected portfolio's own stage model as
 * a board in the middle, the selected program's gate on the right. Every portfolio draws its
 * columns from its own stages, so three portfolios show three different boards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, button, Control, divV, divH, signal, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';
import {portalDb, PortfolioRow, ProgramRow, StageModelRow, StageRow, StageTransitionRow} from '../generated/db';
import {programContext, ProgramContext, programData, readinessFor, recordDecision, setCriterion, stagesOf,
  summarizePlan, syncPortfolio, syncWorkUnits, verdictFor, enabledWorkUnits, workUnitsOf, compoundsOf, syncCompounds} from '../portal';
import {COMPOUND_ROLES, summarizeCompoundPlan} from '../model/compounds';
import {summarizeWorkUnitPlan, WORK_UNIT_LABELS, WorkUnitKind} from '../model/work-units';
import {CriterionReadiness, Readiness} from '../model/readiness';
import {DataRow} from '../model/conditions';
import {Decision, nextStage} from '../model/gates';
import {openPortfolioTable} from './portfolio-table-view';
import {existingSpaces} from '../workspace';

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
  /** The selected program's dataset, read once per open; undefined = not read yet. */
  private data: {rows: DataRow[] | null, error?: string} | undefined;

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
    if (id !== this.selected)
      this.data = undefined;
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
    this.renderDetail(ctx, transitions, successors, readiness, units, compounds, spaceLive);
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
    compounds: {link: any, compound: any}[] = [], spaceLive = false): void {
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
    const criteriaList = named(divV((readiness?.criteria ?? []).map((c) => this.criterionRow(program.id, c, editable)),
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

    const history = named(divV(transitions.map((t) => span(
      `${(t.decided_on ?? t.created_on)?.format?.('YYYY-MM-DD') ?? ''} · ${t.decision}: ${stageName(t.from_stage_id)} → ` +
      `${stageName(t.to_stage_id)}${t.rationale ? ` — ${t.rationale}` : ''}`, 'pp-history-row'))), 'transition-history');

    const lineage: HTMLElement[] = [];
    if (program.predecessor_id)
      lineage.push(named(span(`Continues a program of another portfolio`, 'pp-lineage'), 'lineage-predecessor'));
    for (const s of successors)
      lineage.push(named(span(`Continued as ${s.code} in another portfolio`, 'pp-lineage'), 'lineage-successor'));

    const workspace = spaceLive ?
      named(button('Open workspace', () => void this.openWorkspace(program.space_id!)), 'open-workspace-button') :
      named(span(program.space_id ? 'The workspace folder is gone: rebuild it with buildWorkspace' :
        'No workspace yet: build it with buildWorkspace', 'pp-muted'), 'no-workspace');

    this.detail.replaceChildren(
      named(span(program.name, 'pp-title'), 'program-title'),
      divH([badge(program.status ?? 'active', {variant: statusVariant[program.status ?? 'active'] ?? 'default'}),
        span(`${ctx.portfolio.name} · ${program.code}`, 'pp-muted')]),
      workspace,
      ...(program.source_stage ? [span(`Registry stage: ${program.source_stage}`, 'pp-muted')] : []),
      ...lineage,
      track,
      named(span(!viewed ? 'No further stage' : viewed.id === next?.id ? `Gate to ${viewed.name}` :
        `What entering ${viewed.name} takes`, 'pp-heading'), 'gate-heading'),
      ...(score && readiness!.criteria.length ? [score] : []),
      readiness?.criteria.length ? criteriaList : span('No criteria declared', 'pp-muted'),
      ...(readiness?.dataError ? [named(span(`Program data could not be read: ${readiness.dataError}`, 'pp-reason'),
        'data-error')] : []),
      reasons,
      rationale,
      divH([advance, hold, terminate], 'pp-actions'),
      ...(units.length ? [span('Work units', 'pp-heading'), this.workUnits(units)] : []),
      ...(compounds.length ? [span('Compounds', 'pp-heading'), this.compounds(compounds)] : []),
      span('History', 'pp-heading'),
      transitions.length ? history : span('No transitions yet', 'pp-muted'),
    );
  }

  /** The program's campaigns, series, studies and process stages, as the portfolio enables them. */
  private workUnits(units: {kind: WorkUnitKind, rows: Record<string, any>[]}[]): HTMLElement {
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

  /** One criterion of the gate: its state, and for a data condition what the data says. A computed
   * condition can only be waived; "Met" is for what a person records. */
  private criterionRow(programId: string, c: CriterionReadiness, editable: boolean): HTMLElement {
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
    if (!c.result && !c.gap)
      return row;
    const lines = [c.result ? span(c.result.summary, 'pp-muted') : span(c.gap, 'pp-muted')];
    if (c.result && !c.result.met && c.result.gaps.length)
      lines.push(span(c.gap.slice(c.result.summary.length + 2), 'pp-gap'));
    return divV([row, named(divV(lines, 'pp-criterion-detail'), `criterion-detail-${slug(c.code)}`)], 'pp-criterion-block');
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

  async sync(): Promise<void> {
    if (!this.portfolio)
      return;
    const code = this.portfolio.code;
    const pi = DG.TaskBarProgressIndicator.create(`Syncing ${this.portfolio.name}`);
    try {
      const plans = await syncPortfolio(code);
      const units = await syncWorkUnits(code);
      const compounds = await syncCompounds(code);
      grok.shell.info([...plans.map(summarizePlan), ...units.map(summarizeWorkUnitPlan),
        ...(compounds ? [summarizeCompoundPlan(compounds)] : [])].join('\n'));
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
  const table = named(button('Table', () => void openPortfolioTable(portal.portfolioCode)), 'table-button');
  table.title = 'This portfolio as a table, with filters and analytics';
  const all = named(button('All portfolios', () => void openPortfolioTable(null)), 'all-portfolios-button');
  all.title = 'Every portfolio you can see as one table';
  const view = appView({name: 'Project Portal', content, ribbon: [[sync], [table, all]], status: portal.status,
    path: portal.path});
  await portal.load(portfolioCode);
  (view as any).portal = portal;
  return view;
}
