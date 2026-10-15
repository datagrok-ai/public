/** The Project Portal app: the landing page (find a program, the portfolios), a portfolio's own
 * stage model as a board, and a program's page. Every portfolio draws its columns from its own
 * stages, so three portfolios show three different boards. A click on a board card shows the
 * program's summary in the context panel; opening the program puts its page in the main area. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, button, Card, Control, divV, divH, Scope, signal, span} from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';
import {portalDb, PortfolioRow, ProgramRow, StageModelRow, StageRow, StageTransitionRow} from '../generated/db';
import {ContextCriterion, programContext, ProgramContext, programData, readinessFor, setCriterion, stagesOf,
  verdictFor, enabledWorkUnits, workUnitsOf, compoundsOf} from '../portal';
import {COMPOUND_ROLES} from '../model/compounds';
import {WORK_UNIT_LABELS, WorkUnitKind} from '../model/work-units';
import {CriterionReadiness, Readiness, readinessLabel} from '../model/readiness';
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
import {portalIcon} from './icons';
import {ProgramRecommendation, recommendFor, runModelFor} from '../recommend';
import {describeMissing} from '../model/recommend';
import {follow, myFollows, unfollow, whatsNew, WorkUnitRef} from '../new-data';
import {FollowingRow} from '../generated/db';
import {ADVANCED_BOARD_BUTTONS, BoardButton, MORE_SECTIONS, PageSection, Persona, PERSONA_LABELS, PERSONAS, PORTAL_MODES,
  PersonaConfig, PortalMode, QuickActionCode, QuickCard, quickCards, PagePanel, PANEL_LABELS, panelLayout} from '../model/personas';
import {personaIn, runFunctionCard, savedMode, saveMode, savePersona} from '../personas';
import {missingAssayFrame, recommendationFrame} from '../recommend';
import {attentionFrame} from '../attention';
import {isApprover} from '../portal';
import {showsBoardButton} from '../model/keeper-actions';
import {namedRolesOf, resolvePrincipal, setProgramRole} from '../program-roles';
import {PROGRAM_ROLE_LABELS, PROGRAM_ROLES, ProgramRole} from '../model/program-roles';
import {ProgramRoleRow} from '../generated/db';
import {EvidenceLine, evidenceOf, suggestedEvidence} from '../evidence';
import {openSnapshotDialog} from './snapshot-dialog';
import {snapshotsOf, summaryOf, takeSnapshot} from '../snapshots';
import {describeSnapshot} from '../model/snapshots';
import {ProgramSnapshotRow} from '../generated/db';
import {EVIDENCE_LABELS, EvidenceItem, evidenceKey, linkEvidence, validateEvidence} from '../model/evidence';
import {documentCardsOf, documentsOf, setDocumentStatus, vocabularyOf} from '../documents';
import {describeDocument, DocumentCard, DocumentStatus} from '../model/documents';
import {DocumentRow} from '../generated/db';
import {openDocumentFile, openLinkDocumentDialog} from './document-dialog';
import {evidenceChangesOf, evidenceTableOf, generatedOf, generateEvidenceTableFor, storeChecks, summarizeGenerated} from '../evidence-tables';
import {describeChanges, EvidenceChange} from '../model/evidence-tables';
import {programAccess} from '../program-roles';
import {openProvenanceDialog} from './provenance-dialog';
import {openCompareRunsDialog} from './compare-runs-dialog';
import {observationsOf} from '../observations';
import {checkObservationAccess, ObservationView, subjectLabel} from '../model/observations';
import {openObservationDialog, openRetractDialog} from './observation-dialog';
import {DesignSetItem, designRefusal, designSetsOf, selectionSetsOf} from '../design-sets';
import {describeDesignSet, selectionCompounds, summarizeDesignSet} from '../model/design-sets';
import {openDesignSetDialog, openNewDesignSetDialog, openNewSelectionSetDialog} from './design-set-dialog';
import {SelectionSetRow} from '../generated/db';
import {DmtaStatus, dmtaFrame, dmtaOf} from '../dmta';
import {describeCycle} from '../model/dmta';
import {programFeedOf, feedUserNames, seriesCompounds} from '../feed';
import {FeedEntry, FeedKind, filterFeed} from '../model/feed';
import {feedList, feedSelect, kindOptions, openMyFeedDialog} from './feed-dialog';
import {openPredictedVsMeasured, predictionsOf, scoreCompounds} from '../predictions';
import {describeFit, describePrediction, Prediction, predictedVsMeasured, withPredictions} from '../model/predictions';
import {ApprovalRow} from '../generated/db';
import {decideOrPropose, requestsOf, submitDocument, userNames} from '../approvals';
import {approvalOf, describeApproval, needsApproval, reviewsOf} from '../model/approvals';
import {openWaitingDialog, requestItem} from './approvals-dialog';
import {openWorkspaceDialog} from './workspace-dialog';
import {myPins} from '../pinning';
import {PinHolders, pinOf} from '../model/pinning';
import {newestFirstBy} from '../model/time';
import {filterValue, ROW_LIMIT} from '../rows';
import {findablePrograms, recentPrograms, rememberProgram} from '../find';
import {recordGateDecision} from '../decision-records';
import {DECISION_OUTCOMES, DecisionOutcome, describeRecord, OUTCOME_LABELS} from '../model/decision-records';
import {openDecisionSearchDialog} from './decision-search-dialog';
import {FindProgram, findPrograms, matchedBy, MIN_QUERY, resolveRecents} from '../model/find';
import {accessFrame, portfolioAccess} from '../data-access';
import {accessNotice} from '../model/data-access';
import {discussionOf, postToProgram, ProgramDiscussion} from '../discussion';
import {
  Anchor, decisionsCiting, plainMentions, SmeContext, smeContextFrom, Thread, threadTitle, threadUrl,
} from '../model/discussion';
import {openAskSmeDialog, openDiscussDialog} from './discussion-dialog';
import {area, check, choice, choiceValue, composer, labelled, linkTo, table, text} from './u2-controls';
import {programAddress} from '../model/pinning';
import {showProgramInPanel} from './program-panel';

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

const viewedStageName = (ctx: ProgramContext, stageId: string | null) => ctx.stages.find((s) => s.id === stageId)?.name ?? 'the gate';

/** The portfolio the app opens on: the one asked for if the user can see it, else the only one;
 * null with several, so the landing page asks which. */
export function startingPortfolio(codes: string[], asked?: string | null): string | null {
  if (asked && codes.includes(asked))
    return asked;
  return codes.length === 1 ? codes[0] : null;
}

/** A group's name in the possessive: "Viewers'", "Team's". */
export const possessive = (name: string): string => name.endsWith('s') ? `${name}'` : `${name}'s`;

export class PortalView {
  readonly root = named(divH([], 'pp-root'), 'portal');
  /** The landing page where there are several portfolios: find a program, the recent ones, and a
   * card per portfolio. A portfolio's page lists no other portfolio; the Browse tree switches. */
  private readonly portfolioList = named(divV([], 'pp-portfolios'), 'portfolio-list');
  private readonly board = named(divV([], 'pp-board'), 'stage-board');
  /** The program page, in the main area in place of the board; the context panel stays free for
   * what is selected on it (a run, a document, a compound). */
  private readonly detail = named(divV([], 'pp-detail'), 'program-detail');
  private readonly page = named(divV([this.detail], 'pp-program-page'), 'program-page');
  /** Find a program across portfolios, and the ones opened last: the top of the landing column. */
  private readonly searchInput = text('program-search-input', undefined, {search: true,
    placeholder: 'Find a program: code, name, alias',
    tooltipText: 'Searches every portfolio you can see, by code, name, target or a registry alias'});
  private readonly searchResults = named(divV([], 'pp-search-results'), 'search-results');
  private readonly recentList = named(divV([], 'pp-recent'), 'recent-programs');
  /** The programs the user can see, every portfolio, read once per load. */
  private findable: FindProgram[] = [];
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
  private follows: FollowingRow[] = [];
  /** What landed since the user's last visit, read once per open; undefined = not read yet. */
  private news: {since: string | null, lines: string[]} | null | undefined;
  /** The selected program's dataset, read once per open; undefined = not read yet. */
  private data: {rows: DataRow[] | null, error?: string} | undefined;
  /** What a click's work that failed says: the platform's error balloon, never an unhandled rejection. */
  private readonly failed = (e: any): void => {
    grok.shell.error(e?.message ?? String(e));
  };
  /** The program [data] was read for: a read that lands after another program was opened is not kept. */
  private dataFor: string | null = null;
  /** Simple or Advanced, remembered per user. */
  private mode: PortalMode = savedMode();
  /** The persona whose quick actions head the program page, in the open portfolio. */
  private persona: {persona: Persona, chosen: boolean} | null = null;
  /** In Simple mode, whether More is open on the program page. */
  private moreOpen = false;
  /** Draws the program page again from what was last read: a change of mode or persona reads nothing. */
  private redraw: (() => void) | null = null;
  /** The selected program's roles, and whether the user may change them (a portfolio approver). */
  private roles: {rows: (ProgramRoleRow & {name: string})[], manage: boolean} = {rows: [], manage: false};
  /** What the last role change on this page did, by program id. */
  private roleLine = new Map<string, string>();
  /** The evidence the next decision on [programId] carries: suggested items the decider unticked,
   * and links added by hand. Kept across redraws of the same program. */
  private picked: {programId: string, off: Set<string>, links: EvidenceItem[], snapshot?: boolean} = {programId: '', off: new Set(), links: []};
  /** The selected program's activity, read on each open, and how the Activity section is filtered
   * and paged; the filter is kept while the same program stays open. */
  private feed: {programId: string, entries: FeedEntry[], names: Map<string, string>, error?: string, kind: string,
    series: string, page: number, compounds: Map<string, string[]>} =
    {programId: '', entries: [], names: new Map(), kind: '', series: '', page: 0, compounds: new Map()};
  /** The selected program's stored predictions, the latest per compound, column and model. */
  private predictions: Prediction[] = [];
  /** What the last scoring on this page did, by program id. */
  private scoreLine = new Map<string, string>();
  /** The selected program's approval requests (F-30), newest first, and the names of the people in them. */
  private approvals: {rows: ApprovalRow[], names: Map<string, string>} = {rows: [], names: new Map()};
  /** What the last proposal or review on this page did, by program id. */
  private approvalLine = new Map<string, string>();
  /** The selected program's observations (F-26), who wrote them, and whether the user may add one. */
  private observations: {list: ObservationView[], names: Map<string, string>, refused: string[], approver: boolean} =
    {list: [], names: new Map(), refused: [], approver: false};
  /** The selected program's discussion (F-51): the platform chat on its folder, as threads. */
  private discussion: {d: ProgramDiscussion | null, error?: string, all: boolean} = {d: null, all: false};
  /** Owns the u2 controls of the program page as last drawn; disposed when it is drawn again or cleared. */
  private detailScope: Scope | null = null;
  /** Counts the program page's loads: only the latest draws, so a slow load begun before a change
   * (a follow from the context panel) cannot draw over the reload that came after it. */
  private pageLoads = 0;
  /** Owns the u2 controls of the board as last drawn: its cards, mode and persona. */
  private boardScope: Scope | null = null;
  /** The program page last drawn: Discuss and Ask SME beside a criterion, run or document use it. */
  private pageCtx: ProgramContext | null = null;
  /** The selected program's documents by id, for the context a question about one carries. */
  private pageDocs = new Map<string, DocumentRow>();
  /** The selected program's design and selection sets (F-23), and why the user may not add one. */
  private designs: {sets: DesignSetItem[], selections: SelectionSetRow[], refused: string[]} = {sets: [], selections: [], refused: []};
  /** The selected program's DMTA status (F-24); null where the portfolio is not discovery. */
  private dmta: {status: DmtaStatus | null, error?: string} = {status: null};
  /** The pins in the user's Workspace, by entity id (program records, and folders pinned before
   * 0.3.0): read with the board and on each program open. */
  private pins = new Map<string, PinHolders>();
  /** The open portfolio's data access notice, for its approvers and admins; null when all is well. */
  private accessNoticeEl: HTMLElement | null = null;
  /** Ribbon buttons Simple mode hides; the app view hands them over. */
  readonly boardButtons = new Map<BoardButton, HTMLElement>();
  /** Whether the user is one of the open portfolio's approvers (admins are): the keeper's buttons. */
  keeper = false;

  constructor() {
    this.root.append(this.portfolioList, this.board, this.page);
    this.page.style.display = 'none';
    this.searchInput.root.classList.add('pp-grow');
    let timer: any = null;
    this.searchInput.effect(() => {
      this.searchInput.value.value;
      clearTimeout(timer);
      timer = setTimeout(() => this.renderSearch(), 150);
    });
    this.searchInput.root.addEventListener('keydown', (e) => {
      if (e.key === 'Enter')
        (this.searchResults.querySelector('.pp-search-hit') as HTMLElement | null)?.click();
      else if (e.key === 'Escape') {
        this.searchInput.value.value = '';
        clearTimeout(timer);
        this.renderSearch();
      }
    });
  }

  /** The program page, which opening a program puts in the main area. */
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

  /** Shows [programId] of [portfolioCode] at a glance in the context panel, switching portfolios first when needed. */
  async previewProgram(portfolioCode: string, programId: string): Promise<void> {
    if (this.portfolio?.code !== portfolioCode)
      await this.openPortfolio(portfolioCode);
    else
      this.toBoard();
    await this.selectProgram(programId);
  }

  async load(portfolioCode?: string): Promise<void> {
    const [portfolios, models] = await Promise.all([portalDb.portfolios.query({limit: ROW_LIMIT}), portalDb.stageModels.query({limit: ROW_LIMIT})]);
    this.portfolios = portfolios.sort((a, b) => a.name.localeCompare(b.name));
    this.models = new Map(models.map((m) => [m.id, m]));
    const counts = await Promise.all(this.portfolios.map((p) =>
      portalDb.programs.count({filter: `portfolio_id = ${filterValue(p.id)}`}).catch(() => 0)));
    this.findable = await findablePrograms().catch(() => []);
    this.renderPortfolios(counts);
    const first = startingPortfolio(this.portfolios.map((p) => p.code), portfolioCode);
    if (first)
      await this.openPortfolio(first);
    else if (this.portfolios.length)
      this.showLanding();
    else
      this.status.value = 'No portfolios are loaded; run ProjectPortalCore:loadPortfolios';
  }

  /** Whether there are portfolios to choose between, so a portfolio's page links back to the choice. */
  get severalPortfolios(): boolean {
    return this.portfolios.length > 1;
  }

  /** The landing page: no portfolio open, a card for each. */
  showLanding(): void {
    this.portfolio = null;
    this.selected = null;
    this.redraw = null;
    this.clearDetailScope(true);
    this.clearBoard();
    this.keeper = false;
    this.portfolioList.style.display = '';
    this.board.style.display = 'none';
    this.page.style.display = 'none';
    this.applyBoardButtons();
    this.path.value = APP_PATH;
    this.renderRecent();
    this.status.value = `${this.portfolios.length} portfolios; choose one`;
  }

  /** The search's results under the box: each opens its program, in its portfolio. */
  private renderSearch(): void {
    const q = this.searchInput.value.value;
    const hits = findPrograms(q, this.findable);
    this.searchResults.replaceChildren(...(q.trim().length < MIN_QUERY ? [] : hits.length ? hits.map((h) => {
      const why = matchedBy(h);
      return this.programLink(h.program, 'search-hit', why ? `${h.program.portfolioName} · ${why}` : h.program.portfolioName);
    }) : [named(span(`Nothing matches "${q.trim()}" in the portfolios you can see`, 'pp-muted'),
      'search-empty')]));
  }

  /** The programs the user opened last that they can still see. */
  private renderRecent(): void {
    const recent = resolveRecents(recentPrograms(), this.findable);
    this.recentList.replaceChildren(...(recent.length ? [span('Recent', 'pp-heading'),
      ...recent.map((p) => this.programLink(p, 'recent', p.portfolioName))] : []));
  }

  /** A line that opens [p] in its portfolio; named `<kind>-<portfolio>-<code>`. */
  private programLink(p: FindProgram, kind: string, sub: string): HTMLElement {
    const el = named(divV([divH([portalIcon('program'), span(p.name, 'pp-card-title')], 'pp-icon-line'), span(sub, 'pp-muted')],
      `pp-search-hit pp-${kind}`),
      `${kind}-${slug(p.portfolio)}-${slug(p.code)}`);
    el.title = `${p.code} in ${p.portfolioName}`;
    el.onclick = () => void this.showProgram(p.portfolio, p.id).catch(this.failed);
    return el;
  }

  private renderPortfolios(counts: number[]): void {
    this.renderRecent();
    this.portfolioList.replaceChildren(named(divV([divH([this.searchInput.root]), this.searchResults], 'pp-search-box'), 'program-search'),
      this.recentList, span('Portfolios', 'pp-heading'));
    const cards = named(divH([], 'pp-portfolio-cards'), 'portfolio-cards');
    this.portfolioList.append(cards);
    this.portfolios.forEach((p, i) => {
      const model = this.models.get(p.stage_model_id);
      const card = new Card({icon: portalIcon('portfolio'), title: p.name,
        subtitle: `${counts[i]} ${(p.program_noun ?? 'program').toLowerCase()}s · ${model?.name ?? ''}`,
        onClick: () => void this.openPortfolio(p.code).catch(this.failed)});
      card.root.classList.add('pp-portfolio');
      cards.append(named(card.root, `portfolio-${slug(p.code)}`));
    });
  }

  async openPortfolio(code: string): Promise<void> {
    const p = this.portfolios.find((x) => x.code === code);
    if (!p)
      return;
    // another portfolio's board would stay on show, dead, while this one is read
    if (this.portfolio && this.portfolio !== p)
      this.clearBoard();
    this.portfolio = p;
    this.selected = null;
    this.redraw = null;
    this.clearDetailScope(true);
    const persona = await personaIn(p).catch(() => ({persona: 'lead' as Persona, chosen: false}));
    // another portfolio was opened meanwhile
    if (this.portfolio !== p)
      return;
    this.persona = persona;
    this.keeper = await isApprover(p).catch(() => false);
    if (this.portfolio !== p)
      return;
    this.applyBoardButtons();
    this.path.value = `${APP_PATH}/${p.code}`;
    this.portfolioList.style.display = 'none';
    this.page.style.display = 'none';
    this.board.style.display = '';
    const [stages, programs, pins] = await Promise.all([
      stagesOf(p.stage_model_id),
      portalDb.programs.query({filter: `portfolio_id = ${filterValue(p.id)}`, limit: ROW_LIMIT}),
      myPins().catch(() => new Map<string, PinHolders>()),
    ]);
    if (this.portfolio !== p)
      return;
    [this.stages, this.programs, this.pins] = [stages, programs, pins];
    this.accessNoticeEl = null;
    this.renderBoard();
    void this.checkDataAccess(p);
    this.detail.replaceChildren(span(`Select a ${(p.program_noun ?? 'program').toLowerCase()}`, 'pp-muted'));
    const model = this.models.get(p.stage_model_id);
    this.status.value = `${p.name}: ${this.programs.length} ${(p.program_noun ?? 'program').toLowerCase()}s on ` +
      `${model?.name ?? '?'} ${model?.model_version ?? ''} (${this.stages.length} stages); stages owned by the ` +
      `${p.stage_authority === 'feed' ? 'registry' : 'portal'}`;
  }

  /** The board, drawn in a scope of its own (its cards and the mode and persona inputs): it is drawn again
   * on each open of a portfolio and on a change of persona. */
  private renderBoard(): void {
    this.boardScope?.dispose();
    const scope = this.boardScope = new Scope();
    Scope.runWith(scope, () => this.drawBoard());
  }

  private drawBoard(): void {
    const p = this.portfolio!;
    const model = this.models.get(p.stage_model_id);
    const back: HTMLElement[] = [];
    if (this.severalPortfolios) {
      const toLanding = linkTo('Portfolios ›', 'portfolios-link', () => this.showLanding(),
        {title: 'Back to the choice of portfolio'});
      back.push(toLanding);
    }
    const header = named(divH([
      ...back,
      divH([portalIcon('portfolio', 'Portfolio'), span(p.name, 'pp-title')], 'pp-icon-line'),
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
    this.board.replaceChildren(header, ...(this.accessNoticeEl ? [this.accessNoticeEl] : []), columns);
  }

  /** For the portfolio's approvers and admins, who can share a connection or ask its owner to: a
   * notice when a group cannot query a connection its datasets use. Plain viewers and contributors
   * see nothing; they could not act on it. Read-only, and quiet when the check fails. */
  private async checkDataAccess(p: PortfolioRow): Promise<void> {
    try {
      if (!await isApprover(p))
        return;
      const checks = await portfolioAccess(p);
      const notice = accessNotice(checks);
      if (this.portfolio !== p || !notice)
        return;
      this.accessNoticeEl = named(divH([
        badge('data access', {variant: notice.missing ? 'error' : 'warning'}),
        span(notice.text, 'pp-access-text'),
        named(button('Details', () => setTimeout(() => grok.shell.addTableView(accessFrame(checks)))), 'data-access-details'),
      ], 'pp-access-notice'), 'data-access-notice');
      this.accessNoticeEl.title = 'Share the connection with the group with View and Query (the connection\'s owner or an admin)';
      this.renderBoard();
    } catch (e: any) {
      grok.log.warning(`Project Portal: data access of ${p.code} not checked: ${e?.message ?? e}`);
    }
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
    const card = new Card({title: x.name, onClick: () => void this.selectProgram(x.id).catch(this.failed), children: [
      divH([
        badge(x.status ?? 'active', {variant: statusVariant[x.status ?? 'active'] ?? 'default'}),
        ...(x.modality ? [span(x.modality, 'pp-muted')] : []),
        ...(x.predecessor_id ? [badge('continued', {variant: 'accent'})] : []),
        ...(pinOf(this.pins, {programId: x.id, spaceId: x.space_id}) ? [named(span('pinned', 'pp-pinned-mark'), `pinned-${slug(x.code)}`)] : []),
      ], 'pp-card-meta'),
    ]});
    // the open program is marked here, not with u2's selected state: that one only comes with selectable,
    // where a click toggles it
    card.root.classList.add('pp-card');
    card.root.title = `${x.code}: click for a summary, double-click to open`;
    card.root.addEventListener('dblclick', () => void this.openProgram(x.id).catch(this.failed));
    if (x.id === this.selected)
      card.root.classList.add('pp-selected');
    return named(card.root, `program-${slug(x.code)}`);
  }

  async openProgram(id: string, viewedStage: string | null = null): Promise<void> {
    // a program selected on the board is not open yet: its page is drawn afresh
    const opened = id !== this.selected || this.page.style.display === 'none';
    if (opened) {
      this.showPage();
      void showProgramInPanel(id).catch(this.failed);
      this.data = undefined;
      this.news = undefined;
    }
    if (this.picked.programId !== id)
      this.picked = {programId: id, off: new Set(), links: []};
    this.selected = id;
    this.viewedStage = viewedStage;
    const load = ++this.pageLoads;
    const program = this.programs.find((x) => x.id === id);
    this.markSelected(program ?? null);
    if (program) {
      // .../Portal/<portfolio>/<program>: the program's page can be linked
      if (this.portfolio) {
        this.path.value = `${APP_PATH}/${this.portfolio.code}/${program.code}`;
        if (opened) {
          rememberProgram(this.portfolio.code, program.id);
          this.renderRecent();
        }
      }
    }
    const ctx = await programContext(id);
    const transitions = (await portalDb.stageTransitions.query({filter: `program_id = ${filterValue(id)}`, limit: ROW_LIMIT}))
      .sort(newestFirstBy((x) => x.created_on));
    const successors = await portalDb.programs.query({filter: `predecessor_id = ${filterValue(id)}`, limit: ROW_LIMIT}).catch(() => []);
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
    const known = this.dataFor === id ? this.data : undefined;
    const recommendation = await recommendFor(ctx, {readiness, rows: known?.rows, dataError: known?.error})
      .catch((e) => {
        grok.log.warning(`Project Portal: no model recommendation for ${ctx.program.code}: ${e?.message ?? e}`);
        return null;
      });
    // what landed since the last visit, for a user who follows the program or one of its series
    const follows = await myFollows(id).catch(() => [] as FollowingRow[]);
    // read once per open of a program: a click on the stage track is not a new visit
    if (this.news === undefined && follows.length) {
      const data = await this.programDataOf(ctx);
      const news = await whatsNew(ctx.portfolio, ctx.program, data.rows, follows).catch(() => null);
      if (this.selected === id)
        this.news = news;
    }
    const news = follows.length ? this.news ?? null : null;
    const roles = {rows: await namedRolesOf(id).catch(() => []), manage: await isApprover(ctx.portfolio).catch(() => false)};
    // pinned or not is the platform's favorite, which the Workspace may have changed since
    const pins = await myPins().catch(() => this.pins);
    const snapshots = await snapshotsOf(id).catch(() => [] as ProgramSnapshotRow[]);
    // documents: the cards of the current and next gate's document criteria, and the rest
    const docRows = await documentsOf(id).catch(() => [] as DocumentRow[]);
    // generated evidence tables whose evidence changed since written (F-41); what is found is stored
    const changed = docRows.some((d) => generatedOf(d)) ? await evidenceChangesOf(ctx, docRows).catch(() => null) : null;
    if (changed)
      await storeChecks(docRows, changed).catch(() => undefined);
    const documents = {rows: docRows, cards: await documentCardsOf(ctx, docRows).catch(() => [] as DocumentCard[]),
      contributor: (await programAccess(ctx.portfolio, id).catch(() => null))?.contributor ?? false,
      changes: new Map([...changed ?? new Map()].map(([k, v]) => [k, v.changes]))};
    // the evidence behind each decision, and what the next one would rest on
    const evidence = await evidenceOf(id).catch(() => new Map<string, EvidenceLine[]>());
    const suggested = ctx.approver && ctx.portfolio.stage_authority === 'portal' ? await suggestedEvidence(id, ctx.readiness,
      ctx.stages.find((s) => s.id === ctx.readiness?.stageId)?.name ?? '').catch(() => []) : [];
    // stored predictions, and the program data to set them against
    const predictions = await predictionsOf(id).then((p) => p.latest).catch(() => [] as Prediction[]);
    if (predictions.length)
      await this.programDataOf(ctx);
    // the program's activity, from its own records
    let feedEntries: FeedEntry[] = [];
    let feedNames = new Map<string, string>();
    let feedError: string | undefined;
    try {
      feedEntries = await programFeedOf(ctx.portfolio, ctx.program);
      feedNames = await feedUserNames(feedEntries);
    } catch (e: any) {
      feedError = e?.message ?? String(e);
    }
    // the approval workflow's requests about this program, when the portfolio asks for approvals
    const requests = await requestsOf(id).catch(() => [] as ApprovalRow[]);
    const approvalNames = await userNames(requests.flatMap((r) =>
      [r.proposed_by, ...reviewsOf(r.reviews).map((x) => x.by)])).catch(() => new Map<string, string>());
    // observations (F-26): curated conclusions on the program, its units and compounds
    const observations = await observationsOf(id).catch(() => [] as ObservationView[]);
    const obsAccess = await programAccess(ctx.portfolio, id).catch(() => null);
    // design and selection sets (F-23), where the portfolio has series or the program has sets
    const designSets = await designSetsOf(id).catch(() => [] as DesignSetItem[]);
    const selections = await selectionSetsOf(id).catch(() => [] as SelectionSetRow[]);
    const designRefused = await designRefusal(ctx.program);
    // DMTA status (F-24), discovery
    let dmta: {status: DmtaStatus | null, error?: string} = {status: null};
    if (ctx.portfolio.setting === 'discovery') {
      try {
        dmta = {status: await dmtaOf(ctx.portfolio, ctx.program)};
      } catch (e: any) {
        dmta = {status: null, error: e?.message ?? String(e)};
      }
    }
    const obsNames = await userNames(observations.flatMap((o) => [o.by, o.retractedBy])).catch(() => new Map<string, string>());
    // the discussion (F-51): the platform chat on the program's folder
    const discussion = await discussionOf(ctx.program).then((d) => ({d}))
      .catch((e: any) => ({d: null, error: e?.message ?? String(e)}));
    // another program was opened meanwhile, or this one read again: the page is not this load's to draw
    if (this.selected !== id || load !== this.pageLoads)
      return;
    this.designs = {sets: designSets, selections, refused: designRefused};
    this.dmta = dmta;
    this.discussion = {...discussion, all: opened ? false : this.discussion.all};
    this.observations = {list: observations, names: obsNames, approver: obsAccess?.approver ?? false,
      refused: checkObservationAccess({isApprover: obsAccess?.approver ?? false, isContributor: obsAccess?.contributor ?? false})};
    this.follows = follows;
    this.roles = roles;
    this.pins = pins;
    this.predictions = predictions;
    if (this.feed.programId !== id)
      this.feed = {programId: id, entries: [], names: new Map(), kind: '', series: '', page: 0, compounds: new Map()};
    this.feed.entries = feedEntries;
    this.feed.names = feedNames;
    this.feed.error = feedError;
    this.approvals = {rows: requests, names: approvalNames};
    this.redraw = () => {
      this.clearDetailScope();
      const scope = this.detailScope = new Scope();
      Scope.runWith(scope, () => this.renderDetail(ctx, transitions, successors, readiness, units, compounds, spaceLive,
        runs, choices, attention, recommendation, news, {lines: evidence, suggested}, snapshots, documents));
    };
    this.redraw();
  }

  /** Disposes the program page's controls; with [empty], also takes the page down, so no link or input
   * whose listeners went with the scope stays on show. A redraw replaces the page itself. */
  private clearDetailScope(empty = false): void {
    this.detailScope?.dispose();
    this.detailScope = null;
    if (empty)
      this.detail.replaceChildren();
  }

  /** The same for the board: its cards, links and inputs go with its scope. */
  private clearBoard(): void {
    this.boardScope?.dispose();
    this.boardScope = null;
    this.board.replaceChildren();
  }

  /** Mode and persona, on the board header. */
  private preferences(): HTMLElement {
    // a change settles first: setPersona draws the board, and with it these inputs, again
    const mode = choice('mode-choice', undefined, PORTAL_MODES.map((m) => ({code: m, name: m === 'simple' ? 'Simple' : 'Advanced'})),
      this.mode, {inline: true, onChanged: (v) => queueMicrotask(() => this.setMode(v as PortalMode))});
    const persona = choice('persona-choice', undefined, PERSONAS.map((p) => ({code: p, name: PERSONA_LABELS[p]})),
      this.persona?.persona ?? '', {inline: true, onChanged: (v) => queueMicrotask(() => this.setPersona(v as Persona))});
    const modeRow = labelled('mode-input', 'Mode', mode);
    modeRow.title = 'Simple keeps analytics, runs and history behind More; remembered for you';
    const personaRow = labelled('persona-input', 'Persona', persona);
    personaRow.title = this.persona?.chosen ? 'Your persona: it picks the quick actions; remembered for you' :
      'Picked from your groups in this portfolio; choose one to keep it';
    return named(divH([modeRow, personaRow], 'pp-preferences'), 'preferences');
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

  /** Puts the program page in the main area in place of the board. The board's summary leaves the
   * context panel: the page says all it said. */
  private showPage(): void {
    this.portfolioList.style.display = 'none';
    this.board.style.display = 'none';
    this.page.style.display = '';
    this.page.scrollTop = 0;
  }

  /** Whether [programId]'s page is on screen here. */
  showsPage(programId: string): boolean {
    return this.selected === programId && this.page.style.display !== 'none' && this.root.isConnected;
  }

  /** Whether [portfolioCode]'s board is on screen here. */
  showsBoard(portfolioCode: string): boolean {
    return this.portfolio?.code === portfolioCode && this.board.style.display !== 'none' && this.root.isConnected;
  }

  /** Reads [programId]'s page again when it is on screen, and its board card: after a change made from
   * the context panel (a follow, a favorite, its folder). */
  async refreshProgram(programId: string): Promise<void> {
    if (!this.portfolio || !this.programs.some((x) => x.id === programId))
      return;
    await this.refreshBoard();
    if (this.showsPage(programId))
      await this.openProgram(programId, this.viewedStage);
  }

  /** Back from a program's page to its portfolio's board, the program still marked there. */
  showBoard(): void {
    if (!this.portfolio)
      return;
    this.toBoard();
    const program = this.programs.find((x) => x.id === this.selected);
    if (program)
      void this.selectProgram(program.id).catch(this.failed);
  }

  /** Takes the program page down, if it is up, and puts the open portfolio's board in its place. */
  private toBoard(): void {
    if (!this.portfolio || this.page.style.display === 'none')
      return;
    this.redraw = null;
    this.clearDetailScope(true);
    this.page.style.display = 'none';
    // mode and persona may have changed on the page: the board's header shows them again
    this.renderBoard();
    this.board.style.display = '';
    this.path.value = `${APP_PATH}/${this.portfolio.code}`;
  }

  /** Marks [program]'s card on the board, and only it. */
  private markSelected(program: ProgramRow | null): void {
    for (const el of Array.from(this.board.querySelectorAll('.pp-card')))
      el.classList.remove('pp-selected');
    if (program)
      this.board.querySelector(`[data-u2-name="program-${slug(program.code)}"]`)?.classList.add('pp-selected');
  }

  /** A click on a board card: the program's record in the context panel, its card and Open program.
   * The board stays; a double-click or Open program opens the page. */
  async selectProgram(id: string): Promise<void> {
    const program = this.programs.find((x) => x.id === id);
    if (!program || !this.portfolio)
      return;
    this.selected = id;
    this.markSelected(program);
    await showProgramInPanel(id);
  }

  /** Reads the open portfolio's programs and pins again and redraws its board, wherever the user is:
   * after a decision on a program's page, the page stays. */
  private async refreshBoard(): Promise<void> {
    const p = this.portfolio;
    if (!p)
      return;
    const [programs, pins] = await Promise.all([
      portalDb.programs.query({filter: `portfolio_id = ${filterValue(p.id)}`, limit: ROW_LIMIT}),
      myPins().catch(() => this.pins),
    ]);
    if (this.portfolio !== p)
      return;
    [this.programs, this.pins] = [programs, pins];
    this.renderBoard();
  }

  /** Simple mode leaves the advanced buttons (sync history, export, all portfolios) off the ribbon;
   * Sync, Workspace and Export show only to the open portfolio's approvers (and admins). The
   * functions behind them check the role too. */
  applyBoardButtons(): void {
    for (const [key, el] of this.boardButtons) {
      // the landing page has no portfolio for Sync, Workspace and the like to act on
      const forPortfolio = key !== 'all_portfolios';
      el.style.display = (!forPortfolio || this.portfolio) && showsBoardButton(key, {simple: this.mode === 'simple',
        advanced: ADVANCED_BOARD_BUTTONS, isApprover: this.keeper}) ? '' : 'none';
    }
  }

  /** The portfolio's attention rules over the program's data, read once per open; null without rules. */
  private async attentionOf(ctx: ProgramContext): Promise<{result: ProgramAttention | null, error?: string} | null> {
    if (attentionRulesOf(ctx.portfolio).length === 0)
      return null;
    const data = await this.programDataOf(ctx);
    return data.error ? {result: null, error: data.error} : {result: attentionFor(ctx.portfolio, data.rows)};
  }

  /** Readiness for the viewed stage: the context already holds the next stage's. */
  private async readinessOf(ctx: ProgramContext, stageId: string | null): Promise<Readiness | null> {
    if (stageId === null || stageId === ctx.readiness?.stageId)
      return ctx.readiness;
    const data = ctx.criteria.some((c) => c.stageId === stageId && c.condition) ? await this.programDataOf(ctx) :
      this.dataFor === ctx.program.id ? this.data : undefined;
    return readinessFor(ctx, stageId, data?.rows ?? null, data?.error);
  }

  /** The program dataset of the program [ctx] is about, read once per open. Kept only while that
   * program is still the one shown. */
  private async programDataOf(ctx: ProgramContext): Promise<{rows: DataRow[] | null, error?: string}> {
    const id = ctx.program.id;
    if (this.data !== undefined && this.dataFor === id)
      return this.data;
    let data: {rows: DataRow[] | null, error?: string};
    try {
      data = {rows: await programData(ctx.portfolio, ctx.program.code)};
    } catch (e: any) {
      data = {rows: null, error: e?.message ?? String(e)};
    }
    if (this.selected === id) {
      this.data = data;
      this.dataFor = id;
    }
    return data;
  }

  private renderDetail(ctx: ProgramContext, transitions: StageTransitionRow[], successors: ProgramRow[],
    readiness: Readiness | null, units: {kind: WorkUnitKind, rows: Record<string, any>[]}[] = [],
    compounds: {link: any, compound: any}[] = [], spaceLive = false, runs: SavedRunItem[] = [],
    choices = new Map<string, CriterionRuns>(), attention: {result: ProgramAttention | null, error?: string} | null = null,
    recommendation: ProgramRecommendation | null = null, news: {since: string | null, lines: string[]} | null = null,
    evidence: {lines: Map<string, EvidenceLine[]>, suggested: EvidenceItem[]} = {lines: new Map(), suggested: []},
    snapshots: ProgramSnapshotRow[] = [],
    documents: {rows: DocumentRow[], cards: DocumentCard[], contributor: boolean, changes?: Map<string, EvidenceChange[]>} =
    {rows: [], cards: [], contributor: false}): void {
    this.pageCtx = ctx;
    this.pageDocs = new Map(documents.rows.map((r) => [r.id, r]));
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
      step.onclick = () => void this.openProgram(program.id, s.id === next?.id ? null : s.id).catch(this.failed);
      return step;
    }), 'pp-track'), 'stage-track');

    const editable = ctx.approver && ctx.portfolio.stage_authority === 'portal';
    const criteriaList = named(divV((readiness?.criteria ?? []).map((c) => this.criterionRow(program.id, c, editable,
      choices.get(c.id) ? {runs: choices.get(c.id)!, criterion: ctx.criteria.find((x) => x.id === c.id)!} : null)),
    'pp-criteria'), 'gate-criteria');
    // a gate without required criteria is not ready: "no criteria", never 100%
    const score = readiness ? named(divH([
      span(readiness.empty ? 'No required criteria defined' :
        `${readiness.requiredMet} of ${readiness.requiredTotal} required met or waived`),
      badge(readinessLabel(readiness), {variant: readiness.score === null ? 'default' : readiness.score === 1 ? 'success' :
        readiness.score === 0 ? 'error' : 'warning'}),
    ], 'pp-readiness'), 'readiness-score') : null;

    const rationale = area('rationale-input', undefined, {placeholder: 'Rationale', rows: 2});
    rationale.root.classList.add('pp-wide');
    const verdict = verdictFor(ctx, 'advance');
    const reasons = named(divV(verdict.reasons.map((r) => span(r, 'pp-reason'))), 'gate-reasons');
    const approval = approvalOf(ctx.portfolio.approval);
    const waiting = this.approvals.rows.find((r) => r.kind === 'decision' && r.status === 'pending') ?? null;
    const decide = (d: Decision) => async () => {
      try {
        const items = [...evidence.suggested, ...this.picked.links].filter((e) => !this.picked.off.has(evidenceKey(e)));
        // in a portfolio that asks for approval of this decision, it is proposed instead (F-30)
        const out = await decideOrPropose(program.id, d, rationale.value.value, {source: 'approver', items,
          snapshot: !!this.picked.snapshot});
        this.picked = {programId: program.id, off: new Set(), links: []};
        if (out.proposed) {
          this.approvalLine.set(program.id, out.proposed.line);
          grok.shell.info([out.proposed.line, ...(out.proposed.notified ? [`Approvers: ${out.proposed.notified}`] : []),
            ...(out.proposed.unplaced ? [`Visibility: ${out.proposed.unplaced}`] : [])].join('\n'));
          await this.openProgram(program.id);
          return;
        }
        const done = out.recorded!;
        const ev = done.evidence;
        grok.shell.info(`${program.name}: ${d} recorded${ev?.written ? ` with ${ev.written} piece(s) of evidence` : ''}` +
          `${ev?.errors.length ? `\n${ev.errors.join('\n')}` : ''}${done.workspace ? `\n${done.workspace}` : ''}`);
        await this.refreshBoard();
        await this.openProgram(program.id);
      } catch (e: any) {
        grok.shell.error(e?.message ?? String(e));
      }
    };
    const proposes = (d: Decision) => needsApproval(approval, 'decision', d);
    const advance = named(button(`${proposes('advance') ? 'Propose advance' : 'Advance'}${next ? ` to ${next.name}` : ''}`,
      decide('advance'), {primary: true}), 'advance-button');
    advance.disabled = !verdict.allowed || !!waiting;
    const hold = named(button(proposes('hold') ? 'Propose hold' : 'Hold', decide('hold')), 'hold-button');
    const terminate = named(button(proposes('terminate') ? 'Propose termination' : 'Terminate', decide('terminate')),
      'terminate-button');
    hold.disabled = terminate.disabled = !ctx.approver || ctx.portfolio.stage_authority === 'feed' || !!waiting;
    for (const [b, d] of [[advance, 'advance'], [hold, 'hold'], [terminate, 'terminate']] as [HTMLButtonElement, Decision][]) {
      if (waiting)
        b.title = `"${waiting.title}" is waiting for approval`;
      else if (proposes(d))
        b.title = `${ctx.portfolio.code} asks for approval: this proposes it, and another approver approves it`;
    }

    const approvedBy = new Map(this.approvals.rows.filter((r) => r.stage_transition_id).map((r) => [r.stage_transition_id!,
      ` (proposed by ${r.source === 'data' ? 'data' : this.approvals.names.get(r.proposed_by ?? '') ?? 'someone'}, approved by ` +
      `${reviewsOf(r.reviews).filter((x) => x.verdict === 'approve').map((x) => this.approvals.names.get(x.by) ?? 'someone').join(', ')})`]));
    const history = named(divV(transitions.map((t, i) => {
      const when = (t.decided_on ?? t.created_on)?.format?.('YYYY-MM-DD') ?? '';
      const row = span(t.decision === 'record' ? `${when} · ${describeRecord(t.outcome, stageName(t.to_stage_id), t.rationale)}` :
        `${when} · ${t.decision}: ${stageName(t.from_stage_id)} → ` +
        `${stageName(t.to_stage_id)}${t.source === 'data' ? ' (by data)' : ''}${approvedBy.get(t.id) ?? ''}` +
        `${t.rationale ? ` — ${t.rationale}` : ''}`,
        'pp-history-row');
      // what the data said when it advanced the program
      if (t.evidence)
        row.title = t.evidence;
      const lines = evidence.lines.get(t.id) ?? [];
      const head = divH([row, this.provenanceLink('decision', t.id, `provenance-decision-${i}`)], 'pp-unit');
      if (!lines.length)
        return named(head, `decision-${i}`);
      return named(divV([head, named(divV(lines.map((l, j) => divH([this.evidenceLink(ctx, l, `evidence-${i}-${slug(l.row.kind)}-${j}`),
        this.provenanceLink('evidence', l.row.id, `provenance-evidence-${i}-${j}`)], 'pp-unit')),
        'pp-evidence'), `decision-evidence-${i}`)]), `decision-${i}`);
    })), 'transition-history');
    // where the registry owns the stage, a decision is recorded here with its evidence, and nothing moves
    let record: HTMLElement | null = null;
    if (ctx.portfolio.stage_authority === 'feed' && next) {
      const outcome = choice('record-outcome-choice', undefined, DECISION_OUTCOMES.map((o) => ({code: o, name: OUTCOME_LABELS[o]})),
        'go', {inline: true});
      const rec = named(button(`Record decision for ${next.name}`, async () => {
        try {
          const items = [...evidence.suggested, ...this.picked.links].filter((e) => !this.picked.off.has(evidenceKey(e)));
          const out = await recordGateDecision(program.id, {outcome: choiceValue(outcome) as DecisionOutcome, rationale: rationale.value.value,
            items, snapshot: !!this.picked.snapshot});
          this.picked = {programId: program.id, off: new Set(), links: []};
          grok.shell.info([out.line, ...out.evidence.errors].join('\n'));
          await this.openProgram(program.id, this.viewedStage);
        } catch (e: any) {
          grok.shell.error(e?.message ?? String(e));
        }
      }), 'record-decision-button');
      rec.disabled = !ctx.approver;
      rec.title = ctx.approver ? 'Keep the gate review\'s outcome, rationale and evidence with the program. The stage does not move: ' +
        'the registry sets it' : 'Only the portfolio\'s approvers record gate decisions';
      record = named(divV([divH([labelled('record-outcome-input', 'Outcome', outcome), rec], 'pp-actions'),
        named(span('Recording keeps the decision and its evidence here; the stage is set by the registry', 'pp-muted'), 'record-note')]),
      'record-decision');
    }
    const picker = evidence.suggested.length || ctx.approver ?
      this.evidencePicker(program.id, evidence.suggested, viewedStageName(ctx, ctx.readiness?.stageId ?? null)) : null;

    const lineage: HTMLElement[] = [];
    if (program.predecessor_id)
      lineage.push(named(span(`Continues a program of another portfolio`, 'pp-lineage'), 'lineage-predecessor'));
    for (const s of successors)
      lineage.push(named(span(`Continued as ${s.code} in another portfolio`, 'pp-lineage'), 'lineage-successor'));


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
        rationale.root.querySelector('textarea')?.focus();
      }});

    // Simple mode keeps the heavier sections behind More
    // the page's panels (F-31): the persona's composition from the portfolio, else the mode's
    const panels: Record<PagePanel, () => HTMLElement[]> = {
      quick_actions: () => quick,
      data: () => data.length ? [named(divH(data, 'pp-actions pp-data-actions'), 'program-data-actions')] : [],
      gate: () => [
        ...(program.source_stage ? [span(`Registry stage: ${program.source_stage}`, 'pp-muted')] : []),
        ...lineage,
        track,
        gateHeading,
        ...(score && readiness!.criteria.length ? [score] : []),
        readiness?.criteria.length ? criteriaList : span('No criteria declared', 'pp-muted'),
        ...(readiness?.dataError ? [named(span(`Program data could not be read: ${readiness.dataError}`, 'pp-reason'),
          'data-error')] : []),
        reasons,
        rationale.root,
        ...(picker ? [picker] : []),
        divH([advance, hold, terminate], 'pp-actions'),
        ...(record ? [record] : []),
      ],
      approvals: () => this.approvalSection(ctx, approval),
      documents: () => documents.cards.length || documents.rows.length ? this.documents(ctx, documents) : [],
      work_units: () => units.length ? [span('Work units', 'pp-heading'), this.workUnits(units, ctx)] : [],
      compounds: () => compounds.length ? [span('Compounds', 'pp-heading'), this.compounds(compounds)] : [],
      observations: () => this.observationsSection(ctx),
      discussion: () => this.discussionSection(ctx, transitions, evidence.lines),
      design_sets: () => this.designSection(ctx, units, compounds),
      dmta: () => this.dmtaSection(ctx),
      recommendation: () => recommendation ? this.recommendation(recommendation) : [],
      predictions: () => ctx.portfolio.setting === 'discovery' || this.predictions.length ? this.predictionsSection(ctx) : [],
      attention: () => attention ? this.attention(ctx.portfolio.code, program.code, attention) : [],
      saved_runs: () => [
        named(span(`Saved runs (${runs.length})`, 'pp-heading'), 'saved-runs-heading'),
        divH([saveRun], 'pp-actions'),
        ...(this.savedLine.has(program.id) ? [named(span(this.savedLine.get(program.id)!, 'pp-muted'), 'save-run-result')] : []),
        runs.length ? this.savedRuns(runs) : named(span('No runs saved to this program yet', 'pp-muted'), 'no-saved-runs'),
      ],
      roles: () => this.programRoles(ctx),
      snapshots: () => this.snapshots(ctx, snapshots),
      feed: () => this.activity(ctx, units.find((u) => u.kind === 'series')?.rows.map((r) => String(r.code)) ?? []),
      history: () => [
        span('History', 'pp-heading'),
        transitions.length ? history : span('No transitions yet', 'pp-muted'),
      ],
    };
    const layout = panelLayout(this.persona?.persona ?? null, ctx.portfolio.personas as PersonaConfig | undefined,
      this.mode, this.moreOpen);
    // a composed page counts what is behind More and has something to show
    const behind = layout.composed ? layout.more.filter((x) => panels[x]().length) : [];
    const more = layout.composed ? (behind.length ? named(button(this.moreOpen ? 'Less' : `More (${behind.length})`, () => {
      this.moreOpen = !this.moreOpen;
      this.redraw?.();
    }), 'more-button') : null) : this.mode === 'simple' ? named(button(this.moreOpen ? 'Less' : `More (${MORE_SECTIONS.length})`, () => {
      this.moreOpen = !this.moreOpen;
      this.redraw?.();
    }), 'more-button') : null;
    if (more) {
      more.title = this.moreOpen ? (layout.composed ? 'Back to this persona\'s panels' : 'Back to the simple page') :
        layout.composed ? behind.map((x) => PANEL_LABELS[x]).join(', ') :
          'Missing assays and suggested models, attention, saved runs, roles, snapshots, activity, history';
    }
    // More sits after the panels shown first, before the ones it opens
    const firstBehind = layout.shown.findIndex((x) => layout.more.includes(x));
    const body: HTMLElement[] = [];
    layout.shown.forEach((x, i) => {
      if (i === firstBehind && more)
        body.push(divH([more], 'pp-actions'));
      // a panel is kept whole: the page flows its panels into columns where it is wide enough
      const els = panels[x]();
      if (els.length)
        body.push(named(divV(els, 'pp-panel'), `section-${slug(x)}`));
    });
    if (more && firstBehind < 0)
      body.push(divH([more], 'pp-actions'));

    // back to the board, and from there to the choice of portfolio
    const crumbs = named(divH([
      ...(this.severalPortfolios ? [linkTo('Portfolios ›', 'portfolios-link', () => this.showLanding(),
        {title: 'Back to the choice of portfolio'})] : []),
      linkTo(`${ctx.portfolio.name} ›`, 'board-link', () => this.showBoard(), {title: `Back to the ${ctx.portfolio.name} board`}),
    ], 'pp-crumbs'), 'program-crumbs');
    this.detail.replaceChildren(
      crumbs,
      named(divH([
        divH([portalIcon('program', 'Program'), named(span(program.name, 'pp-title'), 'program-title')], 'pp-icon-line'),
        badge(program.status ?? 'active', {variant: statusVariant[program.status ?? 'active'] ?? 'default'}),
        span([...(program.name !== program.code ? [program.code] : []), current?.name ?? 'Unstaged'].join(' · '), 'pp-muted'),
        // the persona composes this page: it is chosen here as on the board
        this.preferences(),
      ], 'pp-page-header'), 'program-header'),
      ...(newsLine ? [newsLine] : []),
      ...(layout.composed ? [named(span(`Panels for ${PERSONA_LABELS[this.persona!.persona]}, from ${ctx.portfolio.code}`, 'pp-muted'),
        'panels-composed')] : []),
      named(divV(body, 'pp-page-body'), 'program-body'),
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
      follow: () => void this.toggleFollow(ctx, null).catch(this.failed),
      open_workspace: () => void this.openWorkspace(program.space_id!).catch(this.failed),
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
      const link = linkTo(name, `${prefix}-${slug(code)}`, () => void showModelCard(code),
        {title: 'Open the model card'});
      return link;
    };
    if (r.missing.length) {
      out.push(named(span(`Missing assays${r.stage ? ` for ${r.stage.name}` : ''} (${r.missing.length})`, 'pp-heading'),
        'missing-assays-heading'));
      const stored = withPredictions(r.missing, this.predictions);
      out.push(named(divV(r.missing.map((m, mi) => {
        const by = m.column ? r.predictors.get(m.column) ?? [] : [];
        const predicted = stored[mi].predicted;
        return named(divV([
          span(describeMissing(m)),
          ...(predicted.length ? [named(span(`Predicted: ${predicted.slice(0, 5).map((p) => `${p.compound} ${describePrediction(p)}`).join('; ')}` +
            (predicted.length > 5 ? `; ${predicted.length - 5} more` : ''), 'pp-muted'),
          `missing-assay-${slug(m.criterion)}-${slug(m.column!)}-predicted`)] : []),
          ...(by.length ? [divH([span('Predicted meanwhile by', 'pp-muted'), ...by.map((x) => modelLink(x.code, x.name, 'predictor-card'))],
            'pp-run-meta')] : []),
        ], 'pp-run'), `missing-assay-${slug(m.criterion)}${m.column ? `-${slug(m.column)}` : ''}`);
      }), 'pp-runs'), 'missing-assays'));
    }
    out.push(named(span(`Suggested models (${r.models.length})`, 'pp-heading'), 'suggested-models-heading'));
    if (!r.models.length) {
      out.push(named(span(`No model in the ${r.portfolio.code} catalogue fits this program's stage, modality and work units`, 'pp-muted'),
        'no-suggested-models'));
      return out;
    }
    out.push(named(divV(r.models.map((m) => {
      const run = named(button('Run', () => void runModelFor(m.code, {portfolio: r.portfolio.code, program: r.program.code,
        series: r.series}).catch((e) => grok.shell.error(e?.message ?? String(e)))), `run-model-${slug(m.code)}`);
      run.disabled = !m.installed;
      run.title = m.installed ? 'Open the model\'s dialog, filled in from this program' : 'The model\'s function is not installed here';
      const scores = m.scores ?? [];
      const score = scores.length && r.portfolio.setting === 'discovery' ? named(button('Score', () => void this.score(r, m.code).catch(this.failed)),
        `score-model-${slug(m.code)}`) : null;
      if (score) {
        score.disabled = !m.installed;
        score.title = `Run once over the program's compounds and store the predictions of ${scores.join(', ')}`;
      }
      return named(divV([
        divH([modelLink(m.code, m.name), ...(m.approval ? [badge(m.approval.replace('_', ' '),
          {variant: m.approval === 'approved' ? 'success' : 'default'})] : []), run, ...(score ? [score] : [])], 'pp-run-meta'),
        ...(m.question ? [named(span(m.question, 'pp-question'), `suggested-model-question-${slug(m.code)}`)] : []),
        ...(m.reasons.length ? [span(m.reasons.join(' · '), 'pp-muted')] : []),
      ], 'pp-run'), `suggested-model-${slug(m.code)}`);
    }), 'pp-runs'), 'suggested-models'));
    return out;
  }

  /** Scores the program's compounds with a model, then reads the page again. */
  private async score(r: ProgramRecommendation, model: string): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create(`Scoring ${r.program.code}`);
    try {
      const done = await scoreCompounds(r.portfolio.code, r.program.code, model);
      this.scoreLine.set(r.program.id, done.summary.split('\n')[0]);
    } catch (e: any) {
      this.scoreLine.set(r.program.id, e?.message ?? String(e));
    } finally {
      pi.close();
    }
    if (this.mode === 'simple')
      this.moreOpen = true;
    await this.openProgram(r.program.id, this.viewedStage);
  }

  /** Predictions: per model and predicted column, how they compare with what was measured, the
   * compounds side by side, and predicted against measured as a scatter plot. */
  private predictionsSection(ctx: ProgramContext): HTMLElement[] {
    const {program, portfolio} = ctx;
    const fits = predictedVsMeasured(this.predictions, this.data?.rows ?? []);
    const out: HTMLElement[] = [named(span(`Predictions (${this.predictions.length})`, 'pp-heading'), 'predictions-heading')];
    if (this.scoreLine.has(program.id))
      out.push(named(span(this.scoreLine.get(program.id)!, 'pp-muted'), 'score-result'));
    if (!fits.length) {
      out.push(named(span('No predictions stored: Score a suggested model that predicts a column of the program data', 'pp-muted'),
        'no-predictions'));
      return out;
    }
    const fmt = (v: number | null) => v === null ? '—' : Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(1);
    for (const f of fits) {
      const key = `${slug(f.model)}-${slug(f.predicts)}`;
      const scatter = named(button('Predicted vs measured', () => void openPredictedVsMeasured(portfolio.code, program.code, f.model, f.predicts)
        .catch((e) => grok.shell.error(e?.message ?? String(e)))), `predicted-vs-measured-${key}`);
      scatter.title = 'The predictions and the measured values as a table, with a scatter plot';
      const runId = this.predictions.find((p) => p.model === f.model && p.predicts === f.predicts && p.modelRunId)?.modelRunId;
      type Pair = typeof f.pairs[number];
      const t = table<Pair>(`predictions-table-${key}`, [
        {header: 'Compound', render: (pr) => pr.compound},
        {header: `Predicted${f.units ? ` (${f.units})` : ''}`, render: (pr) => fmt(pr.predicted)},
        {header: 'Measured', render: (pr) => fmt(pr.measured)},
      ], f.pairs.slice(0, 12)).root;
      out.push(named(divV([
        divH([span(`${f.modelName} ${f.version} · ${f.predicts}`), scatter,
          ...(runId ? [this.provenanceLink('model_run', runId, `provenance-predictions-${key}`)] : [])], 'pp-run-meta'),
        named(span(describeFit(f) + (f.pairs.length > 12 ? `; first 12 of ${f.pairs.length} shown` : ''), 'pp-muted'), `prediction-fit-${key}`),
        t,
      ], 'pp-run'), `predictions-${key}`));
    }
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
      const link = linkTo(it.item, `attention-open-${slug(it.item)}`,
        () => void openAttentionItem(portfolioCode, programCode, it.item, flagged)
        .catch((e) => grok.shell.error(e?.message ?? String(e))), {title: 'Open the program data at this row'});
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
      const b = linkTo(on ? 'Unfollow' : 'Follow', `follow-series-${slug(code)}`,
        () => void this.toggleFollow(ctx, {kind: 'series', code}).catch(this.failed),
        {title: on ? 'Stop being notified of new data in this series' : 'Be notified when new results land in this series'});
      // what is new in the series: the Activity section narrowed to it
      const a = linkTo('Activity', `series-activity-${slug(code)}`,
        () => void this.setFeedSeries(ctx, code, true).catch(this.failed),
        {title: 'What happened in this series: its work-unit changes, new data and compound roles'});
      return [b, a];
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
        ...this.predictions.filter((p) => p.compound === x.compound.code).map((p) =>
          named(span(`${p.predicts} ≈ ${describePrediction(p)}`, 'pp-muted'), `compound-prediction-${slug(x.compound.code)}-${slug(p.predicts)}`)),
      ], 'pp-unit'), `compound-${slug(x.compound.code)}`)),
    ], 'pp-units'), 'program-compounds');
  }

  /** The model runs saved to the program, newest first: title, then the model (its card opens on
   * a click), the stage at run, unit operation, run type, purpose, author's date. */
  private savedRuns(items: SavedRunItem[]): HTMLElement {
    return named(divV(items.map(({run, model, stage, unitOperation, runType}) => {
      const modelLink = model ?
        linkTo(model.name, `model-link-${slug(run.model_name)}`,
          () => void showModelCard(model.code), {title: 'Open the model card'}) :
        named(span(run.model_name, 'pp-muted'), `model-link-${slug(run.model_name)}`);
      if (!model)
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
          this.provenanceLink('model_run', run.id, `provenance-run-${slug(run.title ?? run.run_id ?? run.id)}`),
          ...this.anchorLinks({kind: 'run', id: run.id, label: run.title ?? run.run_id ?? run.model_name},
            `run-${slug(run.title ?? run.run_id ?? run.id)}`),
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
    row.append(...this.anchorLinks({kind: 'criterion', id: c.id, label: c.name}, `criterion-${slug(c.code)}`));
    if (editable) {
      if (c.basis === 'recorded' && c.state !== 'met' && c.state !== 'waived')
        row.append(named(button('Met', () => void this.mark(programId, c.id, 'met').catch(this.failed)), `mark-met-${slug(c.code)}`));
      if (c.state !== 'waived' && c.state !== 'met')
        row.append(named(button('Waive', () => void this.mark(programId, c.id, 'waived').catch(this.failed)), `mark-waived-${slug(c.code)}`));
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
    // Compute2's Compare Runs over the candidates, also across models (F-29)
    const compare = named(button('Compare in Compute2…', () => {
      const program = this.programs.find((p) => p.id === programId)?.code ?? '';
      openCompareRunsDialog({program, criterion: {code: c.code, name: c.name}, runs,
        compareHere: runs.refused.length ? undefined : () => choose.click()});
    }), `compare-runs-${slug(c.code)}`);
    compare.disabled = runs.candidates.length < 2;
    compare.title = runs.candidates.length < 2 ? 'Comparing needs two candidate runs' :
      'The candidate runs in Compute2\'s Compare Runs, also runs of different models';
    return named(divV([
      named(divH([
        runs.selected ? span(`Run: ${describeSelection(runs.selected)}`) :
          span(`No run selected · ${runs.candidates.length} candidate${runs.candidates.length === 1 ? '' : 's'}`, 'pp-muted'),
        choose,
        compare,
      ], 'pp-criterion-run'), `selected-run-${slug(c.code)}`),
      ...(runs.superseded.length ? [named(divV(runs.superseded.map((e) =>
        span(`Superseded: ${describeSelection(e)}`, 'pp-gap'))), `run-history-${slug(c.code)}`)] : []),
    ]), `criterion-run-${slug(c.code)}`);
  }

  /** Waiting for approval (F-30): the program's pending requests, each with Approve and Reject for
   * an approver who did not propose it and Withdraw for the proposer, then the last few decided.
   * Shown in both modes; nothing when the portfolio asks for no approvals and nothing was asked. */
  private approvalSection(ctx: ProgramContext, approval: ReturnType<typeof approvalOf>): HTMLElement[] {
    const {program} = ctx;
    const rows = this.approvals.rows;
    if (!approval && !rows.length)
      return [];
    const me = DG.User.current().id;
    const pending = rows.filter((r) => r.status === 'pending');
    const closed = rows.filter((r) => r.status !== 'pending').slice(0, 3);
    const done = async (line: string) => {
      this.approvalLine.set(program.id, line);
      await this.refreshBoard();
      await this.openProgram(program.id, this.viewedStage);
    };
    const item = (r: ApprovalRow, name: string) => requestItem(r, {names: this.approvals.names, name, onDone: done,
      canReview: ctx.approver && (r.source === 'data' || r.proposed_by !== me) && !reviewsOf(r.reviews).some((x) => x.by === me),
      canWithdraw: r.source !== 'data' && r.proposed_by === me});
    return [
      named(span(`Waiting for approval (${pending.length})`, 'pp-heading'), 'approvals-heading'),
      named(span(describeApproval(approval), 'pp-muted'), 'approval-setting'),
      ...(this.approvalLine.has(program.id) ? [named(span(this.approvalLine.get(program.id)!, 'pp-muted'), 'approval-result')] : []),
      pending.length ? named(divV(pending.map((r, i) => item(r, `approval-${i}`)), 'pp-units'), 'pending-approvals') :
        named(span('Nothing is waiting for approval', 'pp-muted'), 'no-pending-approvals'),
      ...(closed.length ? [named(divV(closed.map((r, i) => item(r, `approval-closed-${i}`)), 'pp-units'), 'closed-approvals')] : []),
    ];
  }

  /** One card per document criterion of the current and the next gate: its state, the documents
   * over it (each opens; a draft has Approve for approvers), what is missing, the run selected for
   * it and the decisions that cite it; Link document adds one. Documents no card shows follow. */
  private documents(ctx: ProgramContext, d: {rows: DocumentRow[], cards: DocumentCard[], contributor: boolean,
    changes?: Map<string, EvidenceChange[]>}): HTMLElement[] {
    const {program} = ctx;
    const canLink = ctx.approver || d.contributor;
    const vocabulary = vocabularyOf(ctx);
    // a portfolio that asks for approval of documents: submit, and another approver approves (F-30)
    const submits = needsApproval(approvalOf(ctx.portfolio.approval), 'document');
    const failed = (e: any) => grok.shell.error(e?.message ?? String(e));
    const reopen = () => this.openProgram(program.id, this.viewedStage);
    const linkButton = (criterion: {id: string, code: string, name: string} | null, name: string) => {
      const b = named(button(criterion ? 'Link document…' : 'Link a document…', () => void openLinkDocumentDialog({program,
        criterion, vocabulary, approver: ctx.approver && !submits, onLinked: reopen}).catch(failed)), name);
      b.disabled = !canLink;
      b.title = canLink ? 'Link a file of the program folder or a path, or upload one, as a document' :
        'Only contributors and approvers link documents';
      return b;
    };
    // the evidence behind a criterion, generated into a new draft version of a document (F-40)
    const generateButton = (c: DocumentCard, name: string) => {
      const b = named(button('Generate evidence table', () => {
        b.disabled = true;
        void generateEvidenceTableFor(program, c.criterionId).then((r) => {
          grok.shell.info(summarizeGenerated(program, r));
          return reopen();
        }).catch((e) => {
          b.disabled = false;
          failed(e);
        });
      }), name);
      b.disabled = !canLink;
      b.title = canLink ? 'The run selected for it, the newest snapshot, the gate\'s data conditions and the decisions citing it, ' +
        'as a table in the program folder, linked as a new draft version' : 'Only contributors and approvers generate documents';
      return b;
    };
    const line = (row: DocumentRow, name: string) => {
      const open = linkTo(describeDocument({title: row.title, version: row.doc_version, status: row.status as DocumentStatus}),
        `${name}-open`, () => void openDocumentFile(row).catch(failed),
        {cls: row.status === 'superseded' ? 'pp-gap' : undefined, title: row.path});
      const items: HTMLElement[] = [open, this.provenanceLink('document', row.id, `provenance-${name}`),
        ...this.anchorLinks({kind: 'document', id: row.id, label: describeDocument({title: row.title, version: row.doc_version,
          status: row.status as DocumentStatus})}, name)];
      // a generated evidence table (F-40) opens as a table too
      if (generatedOf(row)) {
        const table = linkTo('Table', `table-${name}`,
          () => void evidenceTableOf(row).then((df) => grok.shell.addTableView(df)).catch(failed),
          {title: 'Open the evidence table this document was generated with'});
        items.push(table);
        // its evidence changed since it was written: stale, and a new version offered (F-41)
        const changes = d.changes?.get(row.id) ?? [];
        if (changes.length && row.status !== 'superseded') {
          items.push(named(badge('stale', {variant: 'warning'}), `stale-${name}`));
          const why = named(span(describeChanges(changes), 'pp-gap'), `${name}-changes`);
          items.push(why);
          if (canLink && row.status !== 'submitted') {
            const again = linkTo('Regenerate', `regenerate-${name}`, () => void generateEvidenceTableFor(program, generatedOf(row)!.criterion).then((r) => {
              grok.shell.info(summarizeGenerated(program, r));
              return reopen();
            }).catch(failed), {title: 'Generate the table again from the evidence now, as a new draft version'});
            items.push(again);
          }
        }
      }
      if ((row.status === 'draft' || row.status === 'rejected') && submits && canLink) {
        const submit = linkTo('Submit for approval', `submit-${name}`, () => void submitDocument(program, row.id).then((x) => {
          this.approvalLine.set(program.id, x.line);
          return reopen();
        }).catch(failed), {title: `${ctx.portfolio.code} asks for approval of documents: another approver approves it`});
        items.push(submit);
      } else if (row.status === 'draft' && ctx.approver && !submits) {
        const approve = linkTo('Approve', `approve-${name}`,
          () => void setDocumentStatus(program, row.id, 'approved').then(reopen).catch(failed),
          {title: 'Approve this version: meets the criterion unless waived, supersedes the older approved one'});
        items.push(approve);
      }
      return named(divH(items, 'pp-unit'), name);
    };
    const byId = new Map(d.rows.map((r) => [r.id, r]));
    const variant = (s: string): BadgeVariant => s === 'met' ? 'success' : s === 'waived' ? 'accent' :
      s === 'in_progress' ? 'warning' : s === 'not_met' ? 'error' : 'default';
    const cards = d.cards.map((c) => {
      const key = slug(c.code);
      return named(divV([
        divH([badge(c.state.replace('_', ' '), {variant: variant(c.state)}), span(`${c.name}${c.required ? '' : ' (optional)'}`),
          span(`${c.gate === 'current' ? 'current gate' : 'next gate'} · ${c.stage}`, 'pp-muted')], 'pp-criterion'),
        ...c.documents.map((x, i) => line(byId.get(x.id)!, `document-${key}-${i}`)),
        ...(c.missing ? [named(span(c.missing, 'pp-gap'), `document-missing-${key}`)] : []),
        ...(c.run ? [named(span(`Run: ${c.run}`, 'pp-gap'), `document-run-${key}`)] : []),
        ...(c.citedBy.length ? [named(span(`Cited by: ${c.citedBy.join(', ')}`, 'pp-gap'), `document-cited-${key}`)] : []),
        divH([linkButton({id: c.criterionId, code: c.code, name: c.name}, `link-document-${key}`),
          generateButton(c, `generate-evidence-${key}`)], 'pp-actions'),
      ], 'pp-document-card'), `document-card-${key}`);
    });
    const shown = new Set(d.cards.flatMap((c) => c.documents.map((x) => x.id)));
    const others = d.rows.filter((r) => !shown.has(r.id));
    return [
      named(span(`Documents (${d.rows.length})`, 'pp-heading'), 'documents-heading'),
      ...(cards.length ? [named(divV(cards, 'pp-units'), 'document-cards')] : []),
      ...(others.length ? [named(divV([span('Other documents', 'pp-muted'), ...others.map((r, i) => line(r, `other-document-${i}`))],
        'pp-units'), 'other-documents')] : []),
      divH([linkButton(null, 'link-other-document')], 'pp-actions'),
    ];
  }

  /** The program's snapshots, newest first, each opening read-only; Snapshot now takes one. */
  private snapshots(ctx: ProgramContext, rows: ProgramSnapshotRow[]): HTMLElement[] {
    const take = named(button('Snapshot now', async () => {
      try {
        const r = await takeSnapshot(ctx.program);
        grok.shell.info([`${ctx.program.code}: snapshot taken`, ...r.notes].join('\n'));
        await this.openProgram(ctx.program.id, this.viewedStage);
      } catch (e: any) {
        grok.shell.error(e?.message ?? String(e));
      }
    }), 'snapshot-button');
    take.title = 'Freeze the stage, the criteria, the selected runs and the program dataset now; written once';
    return [
      named(span(`Snapshots (${rows.length})`, 'pp-heading'), 'snapshots-heading'),
      divH([take], 'pp-actions'),
      rows.length ? named(divV(rows.map((s, i) => {
        const el = linkTo(`${s.label} · ${describeSnapshot(summaryOf(s))}`, `snapshot-${i}`,
          () => void openSnapshotDialog(s.id).catch((e) => grok.shell.error(e?.message ?? String(e))),
          {cls: 'pp-history-row', title: 'Open the snapshot, read-only'});
        return divH([el, this.provenanceLink('snapshot', s.id, `provenance-snapshot-${i}`)], 'pp-unit');
      })), 'program-snapshots') : named(span('No snapshots yet', 'pp-muted'), 'no-snapshots'),
    ];
  }

  /** What the next decision will rest on: the gate's readiness (one tick for its criteria), the run
   * selected for each criterion, and links added here. Everything is ticked; untick to leave it out. */
  private evidencePicker(programId: string, suggested: EvidenceItem[], stageName: string): HTMLElement {
    const box = (name: string, label: string, keys: string[], title?: string) =>
      check(name, label, keys.every((k) => !this.picked.off.has(k)), {tooltipText: title,
        onChanged: (on) => keys.forEach((k) => on ? this.picked.off.delete(k) : this.picked.off.add(k))}).root;
    const readiness = suggested.filter((e) => e.role === 'readiness');
    const rows: HTMLElement[] = [];
    if (readiness.length) {
      rows.push(box('evidence-readiness', `Readiness of ${stageName}: ${readiness.length} criteri${readiness.length === 1 ? 'on' : 'a'}`,
        readiness.map(evidenceKey), readiness.map((e) => `${e.name}: ${String(e.context?.state ?? '').replace(/_/g, ' ')}`).join('\n')));
    }
    for (const e of suggested.filter((x) => x.kind === 'model_run'))
      rows.push(box(`evidence-run-${slug(String(e.context?.criterion ?? e.target))}`, `Run: ${e.name} (for ${e.context?.criterion ?? '?'})`, [evidenceKey(e)]));
    for (const e of suggested.filter((x) => x.kind === 'file' && x.context?.document))
      rows.push(box(`evidence-document-${slug(String(e.context?.criterion ?? e.target))}`, `Document: ${e.name}` +
        ` (v${e.context?.version ?? '?'}, for ${e.context?.criterion ?? '?'})`, [evidenceKey(e)], e.target));
    this.picked.links.forEach((e, i) => rows.push(box(`evidence-link-${i}`, `${EVIDENCE_LABELS[e.kind]}: ${e.name}`, [evidenceKey(e)], e.target)));
    rows.push(check('evidence-snapshot', 'Snapshot now', !!this.picked.snapshot, {
      tooltipText: 'Freeze the stage, criteria, selected runs and the program dataset now, before the decision, and attach it',
      onChanged: (on) => this.picked.snapshot = on}).root);
    const link = text('evidence-link-input', undefined, {placeholder: 'Link or file path to attach'});
    link.root.classList.add('pp-grow');
    const add = named(button('Attach', () => {
      const item = linkEvidence(link.value.value);
      const problems = link.value.value.trim() ? validateEvidence(item) : ['Type a link or a file path'];
      if (problems.length) {
        grok.shell.warning(problems.join('; '));
        return;
      }
      if (!this.picked.links.some((x) => evidenceKey(x) === evidenceKey(item)))
        this.picked.links.push(item);
      this.picked.off.delete(evidenceKey(item));
      this.redraw?.();
    }), 'add-evidence-link-button');
    add.title = 'Attach a document, report or page to the next decision';
    return named(divV([
      span('Evidence for the decision', 'pp-muted'),
      ...rows,
      named(divH([link.root, add], 'pp-actions pp-role-form'), 'evidence-link'),
    ], 'pp-evidence-picker'), 'decision-evidence');
  }

  /** Opens the provenance of a portal object; a program node in it opens on this page. */
  provenance(kind: string, id: string): void {
    void openProvenanceDialog(kind, id, {openProgram: async (programId) => {
      const p = await portalDb.programs.get(programId);
      if (!p)
        throw new Error('This program is gone, or it is not shared with you');
      await this.showProgram((await portalDb.portfolios.get(p.portfolio_id)).code, programId);
    }})
      .catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  /** A small Provenance link beside a row. */
  private provenanceLink(kind: string, id: string, name: string): HTMLElement {
    return linkTo('Provenance', name, () => this.provenance(kind, id),
      {cls: 'pp-gap', title: 'What this rests on, and what rests on it', stop: true});
  }

  /** One piece of a decision's evidence: what it is, and a click to open it. A target that is gone
   * shows what was recorded. */
  private evidenceLink(ctx: ProgramContext, l: EvidenceLine, name: string): HTMLElement {
    const recorded = () => grok.shell.info([`${EVIDENCE_LABELS[l.row.kind]}: ${l.row.name}`, `target: ${l.row.target}`,
      ...Object.entries(l.row.context ?? {}).map(([k, v]) => `${k}: ${v}`)].join('\n'));
    const open = async () => {
      if (l.state === 'gone')
        return recorded();
      switch (l.row.kind) {
      case 'url':
        window.open(l.row.target, '_blank', 'noopener');
        return;
      case 'criterion': {
        const c = ctx.criteria.find((x) => x.id === l.row.target);
        return c ? this.openProgram(ctx.program.id, c.stageId === ctx.readiness?.stageId ? null : c.stageId) : recorded();
      }
      case 'model_run': {
        const call = await grok.dapi.functions.calls.find(l.row.target).catch(() => null);
        if (call)
          grok.shell.o = call;
        else
          recorded();
        return;
      }
      case 'snapshot':
        await openSnapshotDialog(l.row.target);
        return;
      case 'file':
        if ((l.row.context as any)?.document)
          await openDocumentFile({title: l.row.name, path: l.row.target, doc_version: (l.row.context as any).version,
            status: (l.row.context as any).status});
        else if (/\.(csv|tsv|txt)$/i.test(l.row.target))
          grok.shell.addTableView(await grok.dapi.files.readCsv(l.row.target));
        else
          recorded();
        return;
      default:
        return recorded();
      }
    };
    // a target gone since it was attached still opens what was recorded, in the plain text colour
    return linkTo(l.text, name, () => void open().catch((e) => grok.shell.error(e?.message ?? String(e))), {
      cls: l.state === 'gone' ? 'pp-gap pp-link-gone' : 'pp-gap',
      title: l.state === 'gone' ? 'Gone since it was attached; click for what was recorded' : l.row.target});
  }

  /** The program's own viewers, contributors and approvers, on top of the portfolio's groups. The
   * portfolio's approvers add and remove them; everyone sees them and their own access. */
  private programRoles(ctx: ProgramContext): HTMLElement[] {
    const {program} = ctx;
    const access = ctx.approver ? 'approver' : ctx.programRole ?? null;
    const out: HTMLElement[] = [
      named(span('Program roles', 'pp-heading'), 'program-roles-heading'),
      named(span(`You: ${access ? `${access}${ctx.programRole && ctx.programRole === access ? ' (program role)' : ' (portfolio)'}` :
        'no role beyond the portfolio\'s groups'}`, 'pp-muted'), 'my-access'),
    ];
    const change = async (groupId: string, role: ProgramRole | null) => {
      try {
        const r = await setProgramRole(program, groupId, role);
        this.roleLine.set(program.id, [r.line, ...r.notes].join('; '));
      } catch (e: any) {
        this.roleLine.set(program.id, e?.message ?? String(e));
      }
      await this.openProgram(program.id, this.viewedStage).catch(this.failed);
    };
    out.push(named(divV(this.roles.rows.length ? this.roles.rows.map((r) => named(divH([
      span(r.name),
      badge(PROGRAM_ROLE_LABELS[r.role as ProgramRole] ?? r.role, {variant: r.role === 'approver' ? 'success' : 'accent'}),
      ...(this.roles.manage ? [(() => {
        const remove = linkTo('Remove', `remove-role-${slug(r.name)}`, () => void change(r.group_id, null));
        return remove;
      })()] : []),
    ], 'pp-unit'), `program-role-${slug(r.name)}`)) :
      [span('None: the portfolio\'s groups decide', 'pp-muted')], 'pp-units'), 'program-roles'));
    if (this.roles.manage) {
      // u2 has a user picker and a one-source entity picker, not one over users and groups: a login or group name
      const who = text('role-principal-input', undefined, {placeholder: 'Group or user login'});
      who.root.classList.add('pp-grow');
      const role = choice('role-choice', undefined, PROGRAM_ROLES.map((r) => ({code: r, name: PROGRAM_ROLE_LABELS[r]})), 'viewer',
        {inline: true});
      const add = named(button('Add', async () => {
        try {
          const g = await resolvePrincipal(who.value.value);
          await change(g.id, choiceValue(role) as ProgramRole);
        } catch (e: any) {
          this.roleLine.set(program.id, e?.message ?? String(e));
          this.redraw?.();
        }
      }), 'add-role-button');
      add.title = 'Give this group or user the role on this program, and share its workspace folder to match';
      out.push(named(divH([who.root, labelled('role-input', 'as', role), add], 'pp-actions pp-role-form'), 'add-role'));
    }
    if (this.roleLine.has(program.id))
      out.push(named(span(this.roleLine.get(program.id)!, 'pp-muted'), 'role-result'));
    return out;
  }

  /** DMTA status (F-24): per series, compounds at each stage and the cycle time; the full status
   * per compound opens as a table. Discovery only; nothing when no compound is in the loop. */
  private dmtaSection(ctx: ProgramContext): HTMLElement[] {
    const {status, error} = this.dmta;
    const heading = named(span('DMTA', 'pp-heading'), 'dmta-heading');
    if (error)
      return [heading, named(span(`The DMTA status could not be read: ${error}`, 'pp-reason'), 'dmta-error')];
    if (!status || !status.rows.length)
      return [];
    heading.textContent = `DMTA (${status.rows.length})`;
    // after the click is done: the app view takes the focus back on a click inside it
    const open = named(button('Open status', () => setTimeout(() => grok.shell.addTableView(dmtaFrame(status, `${ctx.program.code} DMTA status`)))),
      'dmta-open-table');
    open.title = 'Each compound: designed, requested, made, tested, analyzed, with dates and sources';
    return [
      heading,
      named(divV(status.cycles.map((c) => named(span(describeCycle(c)), `dmta-series-${slug(c.series)}`)), 'pp-units'), 'dmta-cycles'),
      ...(status.notes.length ? [named(span(status.notes.join('; '), 'pp-muted'), 'dmta-notes')] : []),
      divH([open], 'pp-actions'),
    ];
  }

  /** Design sets and selection sets (F-23): what was proposed per series, ranked, and chosen for
   * synthesis; hit lists picked from campaigns. Shown where the portfolio has series, or sets exist. */
  private designSection(ctx: ProgramContext, units: {kind: WorkUnitKind, rows: Record<string, any>[]}[],
    compounds: {link: any, compound: any}[]): HTMLElement[] {
    const {program} = ctx;
    const d = this.designs;
    const series = units.find((u) => u.kind === 'series')?.rows.map((r) => String(r.code)) ?? [];
    const campaigns = units.find((u) => u.kind === 'campaign')?.rows.map((r) => String(r.code)) ?? [];
    if (!enabledWorkUnits(ctx.portfolio).includes('series') && !d.sets.length && !d.selections.length)
      return [];
    const reopen = () => this.openProgram(program.id, this.viewedStage);
    const newDesign = named(button('New design set…', () => openNewDesignSetDialog(program, series, reopen)), 'new-design-set-button');
    newDesign.disabled = d.refused.length > 0 || !series.length;
    newDesign.title = d.refused.join('; ') || (!series.length ? 'The program has no series' :
      'Proposed structures for a series, from a table open in the platform: scored, ranked, chosen for synthesis');
    const newSelection = named(button('New selection set…', () => openNewSelectionSetDialog(program, campaigns,
      compounds.map((c) => String(c.compound.code)).sort(), reopen)), 'new-selection-set-button');
    newSelection.disabled = d.refused.length > 0 || !compounds.length;
    newSelection.title = d.refused.join('; ') || (!compounds.length ? 'The program holds no compounds' :
      'A hit list: compounds picked from a campaign, with the rationale');
    const sets = d.sets.map((x) => {
      const open = linkTo(x.set.name, `design-set-${slug(x.set.code)}`,
        () => openDesignSetDialog(program, x, d.refused, reopen));
      return named(divH([open, span(`series ${x.set.series_code} · ${x.set.source}`, 'pp-muted'),
        span(describeDesignSet(summarizeDesignSet(x.proposals)))], 'pp-unit'), `design-set-row-${slug(x.set.code)}`);
    });
    const hits = d.selections.map((x) => {
      const list = selectionCompounds(x.compounds);
      const row = named(divH([span(x.name), span(`${x.campaign_code ? `from ${x.campaign_code} · ` : ''}${list.length} compound${list.length === 1 ? '' : 's'}`,
        'pp-muted'), span(x.rationale, 'pp-muted')], 'pp-unit'), `selection-set-${slug(x.code)}`);
      row.title = list.join(', ');
      return row;
    });
    return [
      named(span(`Design sets (${d.sets.length})`, 'pp-heading'), 'design-sets-heading'),
      divH([newDesign, newSelection], 'pp-actions'),
      named(divV(sets.length ? sets : [span('No design sets yet', 'pp-muted')], 'pp-units'), 'design-sets'),
      ...(hits.length ? [named(span(`Selection sets (${hits.length})`, 'pp-heading'), 'selection-sets-heading'),
        named(divV(hits, 'pp-units'), 'selection-sets')] : []),
    ];
  }

  /** Discuss and Ask SME beside a criterion, a saved run or a document (F-51). */
  private anchorLinks(anchor: Anchor, name: string): HTMLElement[] {
    const ctx = this.pageCtx;
    if (!ctx)
      return [];
    const reopen = () => this.openProgram(ctx.program.id, this.viewedStage);
    const discuss = linkTo('Discuss', `discuss-${name}`,
      () => openDiscussDialog({program: ctx.program, anchor, onPosted: reopen}),
      {cls: 'pp-gap', title: 'Start or continue a thread about this in the program\'s discussion', stop: true});
    const ask = linkTo('Ask SME', `ask-sme-${name}`,
      () => openAskSmeDialog({portfolio: ctx.portfolio, program: ctx.program,
      context: this.smeContext(ctx, anchor), onAsked: reopen}),
    {cls: 'pp-gap', title: 'Ask an expert about this; the context goes with the question', stop: true});
    return [discuss, ask];
  }

  /** What a question carries: where the program stands, the item it is about, and its address. */
  private smeContext(ctx: ProgramContext, anchor: Anchor | null): SmeContext {
    const doc = anchor?.kind === 'document' ? this.pageDocs.get(anchor.id) : undefined;
    return smeContextFrom({...ctx, anchor, anchorUrl: doc?.path ?? null,
      appUrl: `${window.location.origin}${programAddress(ctx.portfolio.code, ctx.program.code)}`});
  }

  /** Opens what a thread is about: a criterion on the gate panel, a run's or document's provenance. */
  private openAnchor(ctx: ProgramContext, a: Anchor): void {
    if (a.kind === 'criterion') {
      const c = ctx.criteria.find((x) => x.id === a.id);
      const el = c ? this.detail.querySelector(`[data-u2-name="criterion-${slug(c.code)}"]`) : null;
      if (el)
        el.scrollIntoView({block: 'center'});
      else
        grok.shell.info(`${a.label}: a criterion of another gate; click its stage on the track to see it`);
      return;
    }
    this.provenance(a.kind === 'run' ? 'model_run' : 'document', a.id);
  }

  /** Discussion (F-51): the program's threads, newest first, each with what it is about (a click
   * opens it), its posts and the decisions that cite it; a box to post, Ask SME, and per thread
   * Reply, Cite in the next decision (approvers) and Open in Chats (the platform's page). */
  private discussionSection(ctx: ProgramContext, transitions: StageTransitionRow[],
    evidence: Map<string, EvidenceLine[]>): HTMLElement[] {
    const {program} = ctx;
    const x = this.discussion;
    const heading = named(span('Discussion', 'pp-heading'), 'discussion-heading');
    if (x.error)
      return [heading, named(span(`The discussion could not be read: ${x.error}`, 'pp-reason'), 'discussion-error')];
    if (!x.d || x.d.reason)
      return [heading, named(span(x.d?.reason ?? 'No discussion', 'pp-muted'), 'discussion-unavailable')];
    const d = x.d;
    const reopen = () => this.openProgram(program.id, this.viewedStage);
    heading.textContent = `Discussion (${d.threads.reduce((n, t) => n + t.posts.length, 0)})`;
    const input = composer({name: 'discussion-input', sendLabel: 'Post',
      placeholder: `Post to ${program.code}'s discussion. @ to mention someone`,
      onSend: (markup) => postToProgram(program, markup).then(() => reopen()).catch((e) => {
        this.failed(e);
        throw e;
      })});
    const post = input.root.querySelector<HTMLElement>('.u2-msg-send');
    if (post) {
      post.dataset.u2Name = 'discussion-post-button';
      post.title = 'Posted to the platform chat on the program\'s folder: whoever may open the folder reads it';
    }
    const ask = named(button('Ask SME…', () => openAskSmeDialog({portfolio: ctx.portfolio, program,
      context: this.smeContext(ctx, null), onAsked: reopen})), 'discussion-ask-sme-button');
    const origin = window.location.origin;
    const shown = x.all ? d.threads : d.threads.slice(0, 5);
    const items = shown.map((t: Thread, i) => {
      const url = d.chatId ? threadUrl(origin, t.posts[0].chatId, t.key) : '';
      const a = t.anchor;
      const title = a ? linkTo(threadTitle(t), `thread-${i}-title`,
        () => this.openAnchor(ctx, a), {title: 'Open what this thread is about'}) :
        named(span(threadTitle(t)), `thread-${i}-title`);
      const decisions = url ? decisionsCiting(url, new Map([...evidence].map(([k, v]) => [k, v.map((l) => l.row)]))) : [];
      const led = decisions.map((id) => transitions.find((tr) => tr.id === id)).filter((tr): tr is StageTransitionRow => !!tr)
        .map((tr) => `${(tr.decided_on ?? tr.created_on)?.format?.('YYYY-MM-DD') ?? ''} · ${tr.decision === 'record' ?
          describeRecord(tr.outcome, ctx.stages.find((s) => s.id === tr.to_stage_id)?.name ?? '', null) : tr.decision}`);
      const posts = t.posts.slice(-3).map((p, j) => named(divV([
        divH([span(p.byName ?? '', 'pp-muted'), span(p.on.slice(0, 16).replace('T', ' '), 'pp-muted'),
          ...(p.question ? [badge('question', {variant: 'warning'})] : [])], 'pp-unit'),
        span(((b) => b.length > 400 ? `${b.slice(0, 399)}…` : b)(plainMentions(p.body)), 'pp-post-body'),
      ], 'pp-post'), `thread-${i}-post-${j}`));
      const reply = linkTo('Reply', `thread-${i}-reply`,
        () => openDiscussDialog({program, anchor: t.anchor, onPosted: reopen}));
      const actions: HTMLElement[] = [reply];
      if (ctx.approver && url) {
        const cite = linkTo('Cite in the next decision', `thread-${i}-cite`, () => {
          const item: EvidenceItem = {kind: 'url', target: url, name: `Discussion: ${threadTitle(t)}`, role: 'supporting'};
          if (!this.picked.links.some((l) => evidenceKey(l) === evidenceKey(item)))
            this.picked.links.push(item);
          this.picked.off.delete(evidenceKey(item));
          grok.shell.info(`${threadTitle(t)}: attached to the next decision`);
          this.redraw?.();
        }, {title: 'Attach this thread to the next gate decision\'s evidence; the thread then shows the decision'});
        actions.push(cite);
      }
      if (url) {
        const chats = linkTo('Open in Chats', `thread-${i}-chats`, () => window.open(url, '_blank', 'noopener'),
          {title: 'The platform\'s own page of this chat'});
        actions.push(chats);
      }
      return named(divV([
        divH([title, span(`${t.posts.length} post${t.posts.length === 1 ? '' : 's'}`, 'pp-muted')], 'pp-unit'),
        ...(led.length ? [named(span(`Led to: ${led.join('; ')}`, 'pp-gap'), `thread-${i}-decisions`)] : []),
        ...(t.posts.length > 3 ? [span(`${t.posts.length - 3} earlier`, 'pp-muted')] : []),
        ...posts,
        divH(actions, 'pp-actions'),
      ], 'pp-observation'), `thread-${i}`);
    });
    const more = d.threads.length > 5 ? [linkTo(x.all ? 'Fewer threads' : `All ${d.threads.length} threads`, 'discussion-all', () => {
      x.all = !x.all;
      this.redraw?.();
    })] : [];
    return [
      heading,
      named(divV([input.root, divH([ask], 'pp-actions')]), 'discussion-form'),
      named(divV(items.length ? items : [span('No posts yet', 'pp-muted')], 'pp-units'), 'discussion-threads'),
      ...more,
    ];
  }

  /** Observations (F-26): curated conclusions, current first; retracted ones stay, marked. Add
   * observation for contributors and approvers; Retract for its author or an approver. */
  private observationsSection(ctx: ProgramContext): HTMLElement[] {
    const {program} = ctx;
    const o = this.observations;
    const me = DG.User.current().id;
    const reopen = () => this.openProgram(program.id, this.viewedStage);
    const add = named(button('Add observation…', () => void openObservationDialog({portfolio: ctx.portfolio, program,
      onAdded: reopen}).catch(this.failed)), 'add-observation-button');
    add.disabled = o.refused.length > 0;
    add.title = o.refused.join('; ') || 'A curated conclusion on the program, a series, campaign or study, or a compound';
    const current = o.list.filter((x) => x.status === 'current').length;
    const items = o.list.map((x, i) => {
      const who = [o.names.get(x.by ?? '') ?? '', x.on?.slice(0, 10) ?? ''].filter((s) => s).join(', ');
      const head = divH([
        span(`${subjectLabel(x.subject, x.code)} · ${x.title}`, x.status === 'retracted' ? 'pp-muted' : ''),
        ...(x.status === 'retracted' ? [badge('retracted', {variant: 'default'})] : []),
        span(who, 'pp-muted'),
      ], 'pp-unit');
      if (x.status === 'current' && (x.by === me || o.approver)) {
        const retract = linkTo('Retract', `retract-observation-${i}`, () => openRetractDialog(x, reopen));
        head.append(retract);
      }
      return named(divV([
        head,
        span(x.finding, x.status === 'retracted' ? 'pp-muted' : ''),
        ...(x.evidence.length ? [span(`Rests on: ${x.evidence.map((e) => e.name).join('; ')}`, 'pp-muted')] : []),
        ...(x.status === 'retracted' ? [span(`Retracted${x.retractedBy ? ` by ${o.names.get(x.retractedBy) ?? ''}` : ''}` +
          `${x.retractedOn ? ` ${x.retractedOn.slice(0, 10)}` : ''}${x.retraction ? `: ${x.retraction}` : ''}`, 'pp-gap')] : []),
      ], 'pp-observation'), `observation-${i}`);
    });
    return [
      named(span(`Observations (${current})`, 'pp-heading'), 'observations-heading'),
      divH([add], 'pp-actions'),
      named(divV(items.length ? items : [span('No observations yet', 'pp-muted')], 'pp-units'), 'observations'),
    ];
  }

  /** Activity: the program's feed, newest first, filtered by kind and narrowed to one series. */
  private activity(ctx: ProgramContext, series: string[]): HTMLElement[] {
    const f = this.feed;
    const heading = named(span('Activity', 'pp-heading'), 'feed-heading');
    if (f.error)
      return [heading, named(span(`The activity could not be read: ${f.error}`, 'pp-reason'), 'feed-error')];
    const scoped = f.series ? filterFeed(f.entries, {series: f.series, seriesCompounds: f.compounds.get(f.series) ?? []}) : f.entries;
    const shown = filterFeed(scoped, {kinds: f.kind ? [f.kind as FeedKind] : null});
    heading.textContent = `Activity (${shown.length})`;
    const controls: HTMLElement[] = [feedSelect('feed-kind-input', 'Show', kindOptions(scoped), f.kind, (v) => {
      f.kind = v;
      f.page = 0;
      this.redraw?.();
    })];
    if (series.length) {
      controls.push(feedSelect('feed-series-input', 'Series', [{code: '', name: `Whole ${(ctx.portfolio.program_noun ?? 'program').toLowerCase()}`},
        ...series.map((c) => ({code: c, name: c}))], f.series, (v) => void this.setFeedSeries(ctx, v, false)));
    }
    return [heading, named(divH(controls, 'pp-actions'), 'feed-filters'),
      ...feedList({entries: shown, names: f.names, page: f.page, prefix: 'feed',
        onPage: (n) => {
          f.page = n;
          this.redraw?.();
        },
        onTrace: (e) => this.provenance(e.ref!.kind, e.ref!.id)})];
  }

  /** Narrows Activity to one series (its compounds from the program data, read once), or back to the program. */
  private async setFeedSeries(ctx: ProgramContext, series: string, scroll: boolean): Promise<void> {
    const f = this.feed;
    f.series = series;
    f.page = 0;
    if (series && !f.compounds.has(series))
      f.compounds.set(series, await seriesCompounds(ctx.portfolio, ctx.program, series).catch(() => []));
    if (this.mode === 'simple')
      this.moreOpen = true;
    this.redraw?.();
    if (scroll)
      this.detail.querySelector('[data-u2-name="feed-heading"]')?.scrollIntoView({block: 'start'});
  }

  /** Opens the program with [code] in the open portfolio; false when it has none. */
  async openProgramByCode(code: string): Promise<boolean> {
    const p = this.programs.find((x) => x.code === code);
    if (p)
      await this.openProgram(p.id);
    return !!p;
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

  /** Waiting for me: the approval requests the user may review, across portfolios. */
  waitingForMe(): void {
    void openWaitingDialog({openProgram: (code, programId) => this.showProgram(code, programId)})
      .catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  /** My feed: everything the user follows; a program in it opens here. */
  /** Searching the decisions of every portfolio the user can see. */
  decisions(): void {
    void openDecisionSearchDialog({portfolios: this.portfolios.map((p) => ({code: p.code, name: p.name})),
      current: this.portfolio?.code ?? null,
      openProgram: (code, programId) => this.showProgram(code, programId)}).catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  myFeed(): void {
    void openMyFeedDialog({
      openProgram: async (code, programId) => {
        const p = code ? null : await portalDb.programs.get(programId);
        if (!code && !p)
          throw new Error('This program is gone, or it is not shared with you');
        await this.showProgram(code || (await portalDb.portfolios.get(p!.portfolio_id)).code, programId);
      },
      trace: (kind, id) => this.provenance(kind, id),
    }).catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  /** The open portfolio as a portfolio file. */
  async exportPortfolio(): Promise<void> {
    if (this.portfolio)
      await openExportDialog(this.portfolio.code).catch((e) => grok.shell.error(e?.message ?? String(e)));
  }

  /** The open portfolio's workspace: the orphan report, parking and the parked folders. */
  async workspace(): Promise<void> {
    if (this.portfolio)
      await openWorkspaceDialog(this.portfolio).catch((e) => grok.shell.error(e?.message ?? String(e)));
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

/** The portal behind each open app view, by the view's root: the shell hands out a new wrapper
 * for a view each time, so a property set on the view object does not come back. */
const portals = new WeakMap<HTMLElement, PortalView>();

/** The Project Portal app view already open, and its portal; null when none is. */
export function openPortal(): {view: DG.ViewBase, portal: PortalView} | null {
  for (const view of Array.from(grok.shell.views)) {
    const portal = portals.get(view.root);
    if (portal)
      return {view, portal};
  }
  return null;
}

/** The portal on screen: the open app view's, else the one previewed from Browse. */
function shownPortal(): {view: DG.ViewBase, portal: PortalView} | null {
  const preview = grok.shell.preview;
  const previewed = preview ? portals.get(preview.root) : undefined;
  return previewed ? {view: preview!, portal: previewed} : openPortal();
}

/** Browse's Project Portal node: a node per portfolio the user can see. A click opens the
 * portfolio in the portal on screen, or in a new preview of the app. */
export async function portfolioTree(treeNode: DG.TreeViewGroup): Promise<void> {
  const portfolios = (await portalDb.portfolios.query({limit: ROW_LIMIT}).catch(() => [] as PortfolioRow[]))
    .sort((a, b) => a.name.localeCompare(b.name));
  treeNode.items.length = 0;
  for (const p of portfolios) {
    const node = treeNode.item(p.name);
    node.icon = portalIcon('portfolio');
    node.root.dataset.u2Name = `browse-portfolio-${slug(p.code)}`;
    node.root.title = `${p.code} · ${p.setting ?? ''}`;
    node.onSelected.subscribe(() => void openInPortal(p.code).catch((e) => grok.shell.error(e?.message ?? String(e))));
  }
}

async function openInPortal(code: string): Promise<void> {
  const shown = shownPortal();
  if (shown) {
    if (shown.view.root !== grok.shell.preview?.root)
      grok.shell.v = shown.view;
    await shown.portal.openPortfolio(code);
    return;
  }
  grok.shell.addPreview(await portalAppView(code));
}

export async function portalAppView(portfolioCode?: string, programCode?: string): Promise<DG.ViewBase> {
  const portal = new PortalView();
  const content = new Control(portal.root);
  const sync = named(button('Sync', () => void portal.sync()), 'sync-button');
  sync.title = 'Run this portfolio\'s registry feeds';
  const history = named(button('Sync history', () => void portal.syncHistory()), 'sync-history-button');
  history.title = 'This portfolio\'s schedule, and what each sync changed';
  const workspaceButton = named(button('Workspace', () => void portal.workspace()), 'workspace-button');
  workspaceButton.title = 'This portfolio\'s workspace: folders no program claims or whose program ended, and parking them';
  const exportButton = named(button('Export', () => void portal.exportPortfolio()), 'export-button');
  exportButton.title = 'This portfolio as a portfolio file, to diff or load on another instance';
  const table = named(button('Table', () => void openPortfolioTable(portal.portfolioCode)), 'table-button');
  table.title = 'This portfolio as a table, with filters and analytics';
  const all = named(button('All portfolios', () => void openPortfolioTable(null)), 'all-portfolios-button');
  all.title = 'Every portfolio you can see as one table';
  const feed = named(button('My feed', () => portal.myFeed()), 'my-feed-button');
  feed.title = 'What happened in the programs and series you follow, newest first';
  const waiting = named(button('Waiting for me', () => portal.waitingForMe()), 'waiting-button');
  const decisions = named(button('Decisions', () => portal.decisions()), 'decisions-button');
  decisions.title = 'Search gate decisions by program, gate, rationale or evidence, across the portfolios you can see';
  waiting.title = 'Proposed decisions and submitted documents and model cards you may approve';
  const view = appView({name: 'Project Portal', content, ribbon: [[sync, history, workspaceButton, exportButton], [table, all], [feed, waiting, decisions]],
    status: portal.status,
    path: portal.path});
  portal.boardButtons.set('sync', sync);
  portal.boardButtons.set('sync_history', history);
  portal.boardButtons.set('export', exportButton);
  portal.boardButtons.set('workspace', workspaceButton);
  portal.boardButtons.set('all_portfolios', all);
  portal.applyBoardButtons();
  await portal.load(portfolioCode);
  if (programCode && portal.portfolioCode === portfolioCode && !(await portal.openProgramByCode(programCode)))
    portal.status.value = `${portfolioCode} has no program ${programCode}, or it is not shared with you`;
  (view as any).portal = portal;
  portals.set(view.root, portal);
  return view;
}
