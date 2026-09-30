/** The Project Portal app: portfolios on the left, the selected portfolio's own stage model as
 * a board in the middle, the selected program's gate on the right. Every portfolio draws its
 * columns from its own stages, so three portfolios show three different boards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, button, Control, divV, divH, signal, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';
import {portalDb, PortfolioRow, ProgramRow, StageModelRow, StageRow, StageTransitionRow} from '../generated/db';
import {programContext, ProgramContext, recordDecision, setCriterion, stagesOf, summarizePlan,
  syncPortfolio, verdictFor} from '../portal';
import {Decision, nextStage} from '../model/gates';

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

export class PortalView {
  readonly root = named(divH([], 'pp-root'), 'portal');
  private readonly portfolioList = named(divV([], 'pp-portfolios'), 'portfolio-list');
  private readonly board = named(divV([], 'pp-board'), 'stage-board');
  private readonly detail = named(divV([], 'pp-detail'), 'program-detail');
  readonly status = signal('Loading portfolios…');

  private portfolios: PortfolioRow[] = [];
  private models = new Map<string, StageModelRow>();
  private portfolio: PortfolioRow | null = null;
  private stages: StageRow[] = [];
  private programs: ProgramRow[] = [];
  private selected: string | null = null;

  constructor() {
    this.root.append(this.portfolioList, this.board, this.detail);
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

  async openProgram(id: string): Promise<void> {
    this.selected = id;
    for (const el of Array.from(this.board.querySelectorAll('.pp-card')))
      el.classList.remove('pp-selected');
    const program = this.programs.find((x) => x.id === id);
    if (program)
      this.board.querySelector(`[data-u2-name="program-${slug(program.code)}"]`)?.classList.add('pp-selected');
    const ctx = await programContext(id);
    const transitions = (await portalDb.stageTransitions.query({filter: `program_id = "${id}"`}))
      .sort((a, b) => String(b.created_on).localeCompare(String(a.created_on)));
    const successors = await portalDb.programs.query({filter: `predecessor_id = "${id}"`}).catch(() => []);
    this.renderDetail(ctx, transitions, successors);
  }

  private renderDetail(ctx: ProgramContext, transitions: StageTransitionRow[], successors: ProgramRow[]): void {
    const {program, stages} = ctx;
    const current = stages.find((s) => s.id === program.stage_id) ?? null;
    const next = nextStage(stages.map((s) => ({...s, isTerminal: s.is_terminal})), program.stage_id ?? null);
    const stageName = (id?: string | null) => stages.find((s) => s.id === id)?.name ?? '—';

    const track = named(divH(stages.map((s) => {
      const cls = s.id === program.stage_id ? 'pp-step pp-step-current' :
        current && s.ordinal < current.ordinal ? 'pp-step pp-step-done' : 'pp-step';
      const step = span(s.code, cls);
      step.title = s.name;
      return step;
    }), 'pp-track'), 'stage-track');

    const criteria = next ? ctx.criteria.filter((c) => c.stageId === next.id) : [];
    const criteriaList = named(divV(criteria.map((c) => {
      const state = ctx.states[c.id] ?? 'not_started';
      const row = named(divH([
        badge(state, {variant: state === 'met' ? 'success' : state === 'waived' ? 'accent' : state === 'not_met' ? 'error' : state === 'in_progress' ? 'warning' : 'default'}),
        span(`${c.name}${c.required ? '' : ' (optional)'}`),
        span(c.kind.replace('_', ' '), 'pp-muted'),
      ], 'pp-criterion'), `criterion-${slug(c.code)}`);
      if (ctx.approver && ctx.portfolio.stage_authority === 'portal') {
        row.append(named(button('Met', () => void this.mark(program.id, c.id, 'met')), `mark-met-${slug(c.code)}`),
          named(button('Waive', () => void this.mark(program.id, c.id, 'waived')), `mark-waived-${slug(c.code)}`));
      }
      return row;
    }), 'pp-criteria'), 'gate-criteria');

    const rationale = named(document.createElement('textarea'), 'rationale-input');
    rationale.placeholder = 'Rationale';
    rationale.rows = 2;
    const verdict = verdictFor(ctx, 'advance');
    const reasons = named(divV(verdict.reasons.map((r) => span(r, 'pp-reason'))), 'gate-reasons');
    const decide = (d: Decision) => async () => {
      try {
        await recordDecision(program.id, d, null, rationale.value);
        grok.shell.info(`${program.name}: ${d} recorded`);
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

    this.detail.replaceChildren(
      named(span(program.name, 'pp-title'), 'program-title'),
      divH([badge(program.status ?? 'active', {variant: statusVariant[program.status ?? 'active'] ?? 'default'}),
        span(`${ctx.portfolio.name} · ${program.code}`, 'pp-muted')]),
      ...(program.source_stage ? [span(`Registry stage: ${program.source_stage}`, 'pp-muted')] : []),
      ...lineage,
      track,
      span(next ? `Gate to ${next.name}` : 'No further stage', 'pp-heading'),
      criteria.length ? criteriaList : span('No criteria declared', 'pp-muted'),
      reasons,
      rationale,
      divH([advance, hold, terminate], 'pp-actions'),
      span('History', 'pp-heading'),
      transitions.length ? history : span('No transitions yet', 'pp-muted'),
    );
  }

  private async mark(programId: string, criterionId: string, state: 'met' | 'waived'): Promise<void> {
    await setCriterion(programId, criterionId, state);
    await this.openProgram(programId);
  }

  async sync(): Promise<void> {
    if (!this.portfolio)
      return;
    const code = this.portfolio.code;
    const pi = DG.TaskBarProgressIndicator.create(`Syncing ${this.portfolio.name}`);
    try {
      const plans = await syncPortfolio(code);
      grok.shell.info(plans.map(summarizePlan).join('\n'));
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
  const view = appView({name: 'Project Portal', content, ribbon: [[sync]], status: portal.status,
    path: '/apps/ProjectPortalCore/Portal'});
  await portal.load(portfolioCode);
  (view as any).portal = portal;
  return view;
}
