/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {textColumn} from './columns';
import {loadAllPortfolios, loadPortfolioSet, programReadiness as readinessOf, readSources, routingAudit, summarizePlan, summarizeRouting, syncPortfolio as syncPortfolioPlans, syncWorkUnits, syncCompounds, setCompoundRole as setRole, programCompound, portfolioByCode, programContext, stagesOf} from './portal';
import {portalDb} from './generated/db';
import {noNewFoldersReason} from './model/workspace-plan';
import {buildMissingFolders, buildPortfolioWorkspace, summarizeFolderGrowth, summarizeWorkspace, workspaceOrphans as orphansOf} from './workspace';
import {summarizeReadiness} from './model/readiness';
import {summarizeWorkUnitPlan} from './model/work-units';
import {CompoundRole, summarizeCompoundPlan} from './model/compounds';
import {portalAppView, portfolioTree} from './views/portal-view';
import {installStyles} from './views/styles';
import {linkTo} from './views/u2-controls';
import {registerProgramRecordHandler} from './program-record';
import {portfolioTableView, TABLE_SEGMENT} from './views/portfolio-table-view';
import {portfolioTables} from './portfolio-table';
import {criterionRuns, saveModelRun as saveRun, summarizeSavedRun} from './model-runs';
import {RunPurpose} from './model/model-runs';
import {loadCatalogues} from './catalogue';
import {addModelCardPanes, modelCardPanelWidget, modelCardWidget, workspaceModelCardWidget} from './views/model-card';
import {openSaveRunDialog} from './views/save-run-dialog';
import {trackModelRuns} from './save-run';
import {openProgramTable, programTable} from './program-data';
import {attentionFrame, attentionRulesOf, programAttentionOf} from './attention';
import {recordedSync, summarizeDue, syncDue} from './sync-runs';
import {syncHistoryFrame} from './views/sync-history';
import {span} from '@datagrok-libraries/u2';
import {autoAdvance, autoAdvanceOf, summarizeAutoAdvanceResult} from './auto-advance';
import {diffPortfolioRefs, exportPortfolioConfig, exportProfileBundle, importProfileBundle} from './profile';
import {summarizeDiff} from './model/profile';
import {missingAssayFrame, recommendByCode, recommendationFrame} from './recommend';
import {checkNewDataOf, followByCode} from './new-data';
import {namedRolesOf, programRolesFrame, resolvePrincipal, setProgramRole as setRoleOf} from './program-roles';
import {ProgramRole} from './model/program-roles';
import {addEvidenceTo, decisionEvidenceFrame, evidenceItemFor, summarizeAdded} from './evidence';
import {EvidenceKind} from './model/evidence';
import {describeSnapshotRow, snapshotsFrame, snapshotsOf, takeSnapshot} from './snapshots';
import {documentsFrame, documentsOf, linkDocumentTo, loadExampleDocuments as loadExamples} from './documents';
import {describeDocument} from './model/documents';
import {provenanceOf} from './provenance';
import {dataPackageFrame, dataPackageOf} from './data-package';
import {manifestCsv} from './model/data-package';
import {comparePlanFrame, criterionComparePlan} from './compare-runs';
import {addObservation as addObs, observationsFrame, observationsOf, retractObservation as retractObs} from './observations';
import {describeObservation, observationsAbout, ObservationSubject} from './model/observations';
import {chooseForSynthesis as choose, createDesignSet as newDesignSet, createSelectionSet as newSelectionSet, designSetsFrame,
  designSetsOf, loadExampleDesignSets as loadDesignExamples, selectionSetsFrame, selectionSetsOf} from './design-sets';
import {cyclesFrame, dmtaFrame, dmtaOf} from './dmta';
import {DesignSource, describeDesignSet, proposalColumns, proposalsFromRows, summarizeDesignSet} from './model/design-sets';
import {checkEvidenceChangesOf, generateEvidenceTableFor, summarizeGenerated} from './evidence-tables';
import {openProvenanceDialog} from './views/provenance-dialog';
import {predictionsFrame, scoreCompounds as scoreProgram} from './predictions';
import {feedFrame, feedUserNames, myFeedOf, programFeedOf, programOf} from './feed';
import {parseKinds} from './model/feed';
import {proposeDecision as propose, requestsFrame, requestsOf, reviewRequest as review, submitDocument as submitDoc,
  submitModelCard as submitCard, userNames, waitingForMe as waitingOf, withdrawRequest as withdraw} from './approvals';
import {reviewsOf} from './model/approvals';
import {secureWriteTables} from './write-grants';
import {moveRecordsOutOfRoot as moveRecords, moveToSecuredTables as moveSecured, ofMyPortfolios} from './secured-rows';
import {accessFrame, dataAccessChecks} from './data-access';
import {summarizeAccess} from './model/data-access';
import {linkProgramRows as linkRows} from './program-rows';
import {decisionsOf, recordGateDecision} from './decision-records';
import {DecisionOutcome, describeDecision, searchDecisions} from './model/decision-records';
import {findablePrograms, recentPrograms as savedRecents} from './find';
import {FindProgram, findPrograms as find, matchedBy, resolveRecents} from './model/find';
import {parkedFoldersOf, parkedFrame, parkOrphans as park, summarizeParkResult, unparkFolder as unpark} from './parking';
import {ApprovalRow} from './generated/db';
import {pinByCode} from './pinning';
import {pinnedFrame, pinnedTiles} from './pinned-programs';
import {MyProgramsWidget} from './views/my-programs-widget';
export * from './package.g';

import {newestFirstBy} from './model/time';
import {seedExamples as seedExampleHistory} from './examples';
import {askSme as askExpert, discussionFrame, discussionOf, postToProgram, resolveAnchor, smeContextOf} from './discussion';
import {AnchorKind, describeAnchor} from './model/discussion';
import {filterValue, ROW_LIMIT} from './rows';

export const _package = new DG.Package();

/** After a run of a catalogued model: a balloon that opens the save-to-program dialog for it. */
function offerSave(call: DG.FuncCall): void {
  installStyles();
  const link = linkTo('Save to program…', 'save-run-link',
    () => void openSaveRunDialog({call}).catch((e) => grok.shell.error(e?.message ?? String(e))));
  const box = document.createElement('div');
  box.append(span(`${call.func.friendlyName ?? call.func.name} ran. `), link);
  grok.shell.info(box, {timeout: 20});
}

//name: portalInit
//tags: init
//description: Once the package loads: offers to save each run of a catalogued model to a program, adds Save to program to a run's context menu, lets a favorite program open its page from the Workspace, and shows a catalogued model's card in its context panel
export async function portalInit(): Promise<void> {
  registerProgramRecordHandler();
  addModelCardPanes();
  trackModelRuns(offerSave);
  grok.events.onContextMenu.subscribe((e: any) => {
    const item = e?.args?.item;
    if (item instanceof DG.FuncCall && e.args.menu)
      e.args.menu.item('Save to program…', () => void openSaveRunDialog({call: item})
        .catch((err) => grok.shell.error(err?.message ?? String(err))));
  });
}

//name: portalAutostart
//description: Loads the package when the platform starts, so its init runs before a model function is selected: the Model card pane in the context panel and the offer to save a run need it
//meta.role: autostart
export function portalAutostart(): void {
  // the work is in portalInit, which the platform runs once the package loads
}

//name: saveRunToProgram
//description: Opens the save-to-program dialog for a run: portfolio, program, unit operation, run type and purpose, filled in from the model's catalogue
//input: funccall call {description: A call that has run}
export async function saveRunToProgram(call: DG.FuncCall): Promise<void> {
  installStyles();
  await openSaveRunDialog({call});
}

//name: Project Portal
//description: Several portfolios on one instance, each on its own stage-gate model and registry
//tags: app
//meta.icon: images/portal.svg
//input: string path {meta.url: true; optional: true}
//output: view result
export async function projectPortalApp(path?: string): Promise<DG.ViewBase> {
  installStyles();
  const segments = path?.split('/').filter((s) => s.length > 0) ?? [];
  // .../Portal/<code>/table: the portfolio table; .../Portal/table: every portfolio's
  if (segments[segments.length - 1] === TABLE_SEGMENT) {
    const code = segments[segments.length - 2];
    const known = code ? await portalDb.portfolios.getByKey({code}) : null;
    return portfolioTableView(known ? code : null);
  }
  // .../Portal/<code>/<program>: the program's page in its portfolio
  const [prev, last] = segments.slice(-2).map((x) => {
    try {
      return decodeURIComponent(x);
    } catch {
      return x;
    }
  });
  if (prev && last && !(await portalDb.portfolios.getByKey({code: last}).catch(() => null)) &&
    await portalDb.portfolios.getByKey({code: prev}).catch(() => null))
    return portalAppView(prev, last);
  return portalAppView(last ?? prev);
}

//name: projectPortalTreeBrowser
//description: The portfolios you can see, one node each under Project Portal in Browse; a click opens the portfolio
//input: dynamic treeNode
//meta.role: appTreeBrowser
//meta.app: Project Portal
export async function projectPortalTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
  await portfolioTree(treeNode);
}

//name: portfolioTable
//description: The programs of a portfolio (or, without one, of every portfolio you can see) as a table: stage, days in stage, readiness for the next gate, compounds by role, work units, workspace
//input: string portfolio {optional: true}
//output: dataframe result
export async function portfolioTable(portfolio?: string): Promise<DG.DataFrame> {
  return (await portfolioTables(portfolio || undefined)).programs;
}

//name: programData
//description: One program's governed data as a table named after it: the portfolio's program dataset, or with dataset one of the further datasets it declares. The query runs as you, on a connection shared with you
//input: string portfolio
//input: string program
//input: string dataset {optional: true; description: Dataset code; default: the program dataset}
//output: dataframe result
export async function programData(portfolio: string, program: string, dataset?: string): Promise<DG.DataFrame> {
  return (await programTable(portfolio, program, dataset || undefined)).df;
}

// A function registers under its TypeScript name, so programData was callable only as programDataTable
// before 0.3.0; the old name stays for scripts that use it.
//name: programDataTable
//description: The same as programData, under the name earlier builds registered it by
//input: string portfolio
//input: string program
//input: string dataset {optional: true; description: Dataset code; default: the program dataset}
//output: dataframe result
export async function programDataTable(portfolio: string, program: string, dataset?: string): Promise<DG.DataFrame> {
  return programData(portfolio, program, dataset);
}

//name: openProgramData
//description: Opens one program's governed data as a table view named after the program, with the viewers the portfolio declares or the defaults for its setting
//input: string portfolio
//input: string program
//input: string dataset {optional: true; description: Dataset code; default: the program dataset}
export async function openProgramData(portfolio: string, program: string, dataset?: string): Promise<void> {
  await openProgramTable(portfolio, program, dataset || undefined);
}

//name: programAttention
//description: The compounds (or lots) of one program that trip the portfolio's attention rules: one row per item and rule, with severity and message. The program dataset is read as you
//input: string portfolio
//input: string program
//output: dataframe result
export async function programAttention(portfolio: string, program: string): Promise<DG.DataFrame> {
  const r = await programAttentionOf(portfolio, program);
  return attentionFrame(r.attention, attentionRulesOf(r.portfolio), `${r.program.code} · Needs attention`);
}

//name: recommendModels
//description: Catalogued models that fit a program, best first: keyed on the card's appliesTo (task, stage, modality, work-unit kind; a key left out matches any), models predicting an assay the gate lacks first. Stage defaults to the program's next
//input: string portfolio
//input: string program
//input: string stage {optional: true; description: Stage code whose gate the missing assays come from}
//input: string series {optional: true; description: Series code: only its rows, and models for series}
//input: string task {optional: true; description: Only models declaring this task (potency, ADMET, physchem, ...)}
//output: dataframe result
export async function recommendModels(portfolio: string, program: string, stage?: string, series?: string,
  task?: string): Promise<DG.DataFrame> {
  return recommendationFrame(await recommendByCode(portfolio, program, {stage: stage || undefined,
    series: series || undefined, task: task || undefined}));
}

//name: missingAssays
//description: The measurements a program's gate still lacks: per open data condition and column, the rows (compounds, lots) without it, and the suggested models that predict it. Stage defaults to the program's next
//input: string portfolio
//input: string program
//input: string stage {optional: true; description: Stage code}
//input: string series {optional: true; description: Series code: only its rows}
//output: dataframe result
export async function missingAssays(portfolio: string, program: string, stage?: string, series?: string): Promise<DG.DataFrame> {
  return missingAssayFrame(await recommendByCode(portfolio, program, {stage: stage || undefined, series: series || undefined}));
}

//name: checkNewData
//description: New-data notification: reads the program dataset of each followed program of a portfolio once, compares it with each follow's last snapshot (rows, filled values, the watermark column when declared) and tells followers in the platform's notifications. Runs after every sync too
//input: string code {description: Portfolio code}
//output: string report
export async function checkNewData(code: string): Promise<string> {
  return (await checkNewDataOf(await portfolioByCode(code))).summary;
}

//name: followProgram
//description: Follow a program, or one of its series, to be notified when new results land in its program dataset. The data now is the baseline
//input: string portfolio
//input: string program
//input: string series {optional: true; description: Series code}
//output: string report
export async function followProgram(portfolio: string, program: string, series?: string): Promise<string> {
  return followByCode(portfolio, program, series || null);
}

//name: unfollowProgram
//description: Stop following a program, or one of its series
//input: string portfolio
//input: string program
//input: string series {optional: true; description: Series code}
//output: string report
export async function unfollowProgram(portfolio: string, program: string, series?: string): Promise<string> {
  return followByCode(portfolio, program, series || null, true);
}

//name: scoreCompounds
//description: Batch scoring: runs a catalogued predictive model once over a program's dataset, saves the run to the program, and stores one prediction per compound and column the card's appliesTo.predicts names. Idempotent per model version
//input: string portfolio
//input: string program
//input: string model {description: The model's catalogue code, e.g. ProjectPortal:SolubilityEstimate}
//output: string report
export async function scoreCompounds(portfolio: string, program: string, model: string): Promise<string> {
  return (await scoreProgram(portfolio, program, model)).summary;
}

//name: programPredictions
//description: A program's stored predictions, the latest per compound, column and model, with the measured value where the program dataset has one
//input: string portfolio
//input: string program
//input: string model {optional: true; description: Only this model}
//output: dataframe result
export async function programPredictions(portfolio: string, program: string, model?: string): Promise<DG.DataFrame> {
  return predictionsFrame(portfolio, program, model || null);
}

//name: programFeed
//description: A program's activity, newest first, from the portal's own records: gate decisions, evidence, runs saved and selected, documents, snapshots, compound role changes, registry and work-unit changes seen by the sync, new data. Kinds are comma-separated; a series narrows it to what is about that series
//input: string portfolio
//input: string program
//input: string series {optional: true; description: Series code}
//input: string kinds {optional: true; description: Comma-separated: decision, evidence, run, selection, prediction, document, snapshot, compound, work_unit, registry, new_data}
//input: int limit {optional: true; default: 100}
//output: dataframe result
export async function programFeed(portfolio: string, program: string, series?: string, kinds?: string,
  limit?: number): Promise<DG.DataFrame> {
  const {portfolio: pf, program: p} = await programOf(portfolio, program);
  const entries = (await programFeedOf(pf, p, {series: series || null, kinds: parseKinds(kinds)})).slice(0, limit || 100);
  return feedFrame(entries, await feedUserNames(entries), `${p.code} · Activity${series ? ` · ${series}` : ''}`);
}

//name: myFeed
//description: The activity of everything you follow, newest first: each followed program's feed, or only a followed series' part of it. Kinds are comma-separated, as for programFeed
//input: string kinds {optional: true}
//input: int limit {optional: true; default: 100}
//output: dataframe result
export async function myFeed(kinds?: string, limit?: number): Promise<DG.DataFrame> {
  const {entries} = await myFeedOf({kinds: parseKinds(kinds)});
  const shown = entries.slice(0, limit || 100);
  return feedFrame(shown, await feedUserNames(shown), 'My feed');
}

/** A program of a portfolio by code, under row security. */
async function programByCode(portfolio: string, program: string) {
  const pf = await portfolioByCode(portfolio);
  const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(pf.id)} and code = ${filterValue(program)}`, limit: ROW_LIMIT}))[0];
  if (!p)
    throw new Error(`${portfolio} has no program ${program}, or it is not shared with you`);
  return p;
}

//name: setProgramRole
//description: Gives a group, or a user by login, a role on one program (viewer, contributor, approver), on top of the portfolio's groups; none removes it. Only the portfolio's approvers may. The program's workspace folder is shared to match
//input: string portfolio
//input: string program
//input: string principal {description: Group name or user login}
//input: string role {choices: ["viewer", "contributor", "approver", "none"]}
//output: string report
export async function setProgramRole(portfolio: string, program: string, principal: string, role: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const group = await resolvePrincipal(principal);
  const r = await setRoleOf(p, group.id, role === 'none' || !role ? null : role as ProgramRole);
  return [r.line, ...r.notes.map((n) => `  ${n}`)].join('\n');
}

//name: programRoles
//description: The roles of one program: each group or user with its role (viewer, contributor, approver), on top of the portfolio's groups
//input: string portfolio
//input: string program
//output: dataframe result
export async function programRoles(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  return programRolesFrame(await namedRolesOf(p.id), `${p.code} · Program roles`);
}

//name: decisionEvidence
//description: A program's gate decisions with the evidence each rests on, one row per piece: kind, frozen name, target (a soft reference), and whether the target still exists under that name (live, renamed, gone, unchecked). With decision, only that gate-log entry
//input: string portfolio
//input: string program
//input: string decision {optional: true; description: Gate-log entry id (stage_transition); default: every decision}
//output: dataframe result
export async function decisionEvidence(portfolio: string, program: string, decision?: string): Promise<DG.DataFrame> {
  return decisionEvidenceFrame(await programByCode(portfolio, program), decision || null);
}

//name: addEvidence
//description: Adds evidence to a program's gate decision (default: its latest): a saved run by run_id or title, a criterion by code, a URL, a file path. Its name and context are frozen now. The program's approvers only
//input: string portfolio
//input: string program
//input: string kind {choices: ["model_run", "criterion", "url", "file", "snapshot"]}
//input: string target {description: run_id, criterion code, URL, file path or snapshot id}
//input: string name {optional: true; description: Name to freeze; default: the target's own}
//input: string decision {optional: true; description: Gate-log entry id; default: the latest decision}
//output: string report
export async function addEvidence(portfolio: string, program: string, kind: string, target: string, name?: string,
  decision?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const item = await evidenceItemFor(p, kind as EvidenceKind, target, name || undefined);
  return summarizeAdded(p, await addEvidenceTo(p, [item], decision || null));
}

//name: recordGateDecisionRecord
//description: Records a gate decision with its evidence where the registry owns the stage: outcome, rationale and the gate's suggested evidence (readiness, selected runs, documents). The program does not move: the registry sets the stage. The program's approvers only
//input: string portfolio
//input: string program
//input: string outcome {choices: ["go", "hold", "no_go"]}
//input: string rationale
//input: string gate {optional: true; description: Code of the stage the gate leads into; default: the next stage}
//output: string report
export async function recordGateDecisionRecord(portfolio: string, program: string, outcome: string, rationale: string,
  gate?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const gateStageId = gate ? (await stagesOf(p.stage_model_id)).find((s) => s.code === gate)?.id ?? 'unknown' : null;
  const r = await recordGateDecision(p.id, {outcome: outcome as DecisionOutcome, rationale, gateStageId});
  return [r.line, ...r.evidence.errors].join('\n');
}

//name: findDecisions
//description: Gate decisions across the portfolios the current user can see, newest first: words in the rationale, outcome or evidence, a program (code or name), a gate (stage code or name), a portfolio. Registry moves are left out unless asked for
//input: string text {optional: true}
//input: string portfolio {optional: true}
//input: string program {optional: true}
//input: string gate {optional: true}
//input: bool withSync = false
//output: dataframe decisions
export async function findDecisions(text?: string, portfolio?: string, program?: string, gate?: string,
  withSync: boolean = false): Promise<DG.DataFrame> {
  const found = searchDecisions(await decisionsOf(portfolio || null), {text, program, gate, withoutSync: !withSync});
  const df = DG.DataFrame.fromColumns([
    textColumn('portfolio', found.map((d) => d.portfolio)),
    textColumn('program', found.map((d) => d.programCode)),
    DG.Column.fromStrings('on', found.map((d) => d.on.slice(0, 10))),
    textColumn('decision', found.map((d) => d.decision)),
    textColumn('outcome', found.map((d) => d.outcome ?? '')),
    textColumn('gate', found.map((d) => d.to?.code ?? '')),
    textColumn('rationale', found.map((d) => d.rationale ?? '')),
    textColumn('evidence', found.map((d) => d.evidence.map((e) => e.name).join('; '))),
    textColumn('line', found.map(describeDecision)),
  ]);
  df.name = 'Decisions';
  return df;
}

//name: snapshotProgram
//description: Freezes a program now: stage, status, the state of each criterion, the runs selected, and the program dataset as a CSV file in the program's workspace folder, fingerprinted. Written once. Contributors and approvers
//input: string portfolio
//input: string program
//input: string label {optional: true; description: Default: program code and time}
//input: string note {optional: true}
//output: string report
export async function snapshotProgram(portfolio: string, program: string, label?: string, note?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const r = await takeSnapshot(p, {label: label || undefined, note: note || undefined});
  return [`${p.code}: snapshot "${r.snapshot.label}" taken: ${describeSnapshotRow(r.snapshot)}` +
    (r.snapshot.table_file ? `; table in ${r.snapshot.table_file}` : ''), ...r.notes.map((n) => `  ${n}`)].join('\n');
}

//name: programSnapshots
//description: A program's snapshots, newest first: when, label, stage, status, summary, the table file
//input: string portfolio
//input: string program
//output: dataframe result
export async function programSnapshots(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  return snapshotsFrame(await snapshotsOf(p.id), `${p.code} · Snapshots`);
}

//name: programDocuments
//description: A program's documents: the criterion each answers, kind, title, version, status, the file's path
//input: string portfolio
//input: string program
//output: dataframe result
export async function programDocuments(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  const ctx = await programContext(p.id);
  return documentsFrame(await documentsOf(p.id), ctx.criteria, `${p.code} · Documents`);
}

//name: linkDocument
//description: Links a file to a program as a document, or changes the document over that file. Approving meets the document criterion it answers, unless waived, and supersedes the older approved version. Contributors link drafts; approvers approve
//input: string portfolio
//input: string program
//input: string path {description: The file's full path: System:AppData/..., or a file in the program's folder}
//input: string criterion {optional: true; description: The document criterion it answers, by code}
//input: string kind {optional: true; description: From the portfolio's document kinds; default: the criterion's}
//input: string title {optional: true}
//input: string version {optional: true}
//input: string status {optional: true; choices: ["draft", "approved", "superseded"]; description: Default: draft, or the document's status now}
//output: string report
export async function linkDocument(portfolio: string, program: string, path: string, criterion?: string, kind?: string,
  title?: string, version?: string, status?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const r = await linkDocumentTo(p, {path, criterion: criterion || null, kind: kind || undefined, title: title || undefined,
    version: version || undefined, status: (status || undefined) as any});
  return [`${p.code}: ${describeDocument({title: r.document.title, version: r.document.doc_version, status: r.document.status as any})} ` +
    `${r.created ? 'linked' : 'changed'}` + (r.criterion ? `; ${r.criterion.code} ${r.criterion.state}` : '') +
    (r.superseded ? `; ${r.superseded} older version(s) superseded` : ''), ...r.notes.map((n) => `  ${n}`)].join('\n');
}

//name: generateEvidenceTable
//description: Generates the evidence table of a program's document criterion (the run selected for it with its inputs and outputs, the newest snapshot, the gate's data conditions, the decisions citing it) as a CSV and a Markdown file in the program's workspace folder, linked as a new draft version of a document over the criterion. Contributors and approvers
//input: string portfolio
//input: string program
//input: string criterion {description: The document criterion, by code}
//output: dataframe result
export async function generateEvidenceTable(portfolio: string, program: string, criterion: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  const r = await generateEvidenceTableFor(p, criterion);
  r.df.setTag('portal.report', summarizeGenerated(p, r));
  return r.df;
}

//name: checkEvidenceChanges
//description: Compares each generated evidence table of a portfolio's programs with the evidence now (another run selected, a newer snapshot, a data verdict changed, a decision added or gone), marks it stale or current, and tells a stale document's author and the program's followers once in the platform's notifications. The sync does this after every run
//input: string code
//output: string report
export async function checkEvidenceChanges(code: string): Promise<string> {
  return (await checkEvidenceChangesOf(await portfolioByCode(code))).summary;
}

//name: provenance
//description: What a portal object rests on and what rests on it, depth steps each way: data, runs, models, documents, criteria, decisions, the program. Returns {start, nodes, edges}; an edge runs from what a node rests on to the node
//input: string kind {choices: ["program", "decision", "evidence", "model_run", "snapshot", "document", "compound"]}
//input: string id {description: The record's id, or a run's run_id, or a compound's code}
//input: int depth {optional: true; default: 4}
//output: object result
export async function provenance(kind: string, id: string, depth?: number): Promise<object> {
  const p = await provenanceOf(kind, id, depth || 4);
  return {start: p.start, nodes: p.nodes, edges: p.edges};
}

//name: dataPackage
//description: A program's data package manifest: one line per criterion (to the next gate), document, evidence table, saved run, dataset and snapshot, with stage, criterion, status, version, source, by, on and a reference. Read from the records now
//input: string portfolio
//input: string program
//output: dataframe result
export async function dataPackage(portfolio: string, program: string): Promise<DG.DataFrame> {
  return dataPackageFrame(await dataPackageOf((await programByCode(portfolio, program)).id));
}

//name: dataPackageCsv
//description: A program's data package manifest as CSV text, as the program page's Download CSV saves it
//input: string portfolio
//input: string program
//output: string csv
export async function dataPackageCsv(portfolio: string, program: string): Promise<string> {
  return manifestCsv((await dataPackageOf((await programByCode(portfolio, program)).id)).lines);
}

//name: criterionRunsToCompare
//description: A run criterion's candidate runs, grouped by model, as Compute2's Compare Runs would take them (F-29): model, run, run id, saved on, selected. Compare Runs takes no runs as input, so this is the list to add there
//input: string portfolio
//input: string program
//input: string criterion {description: The criterion code, e.g. process_description}
//output: dataframe result
export async function criterionRunsToCompare(portfolio: string, program: string, criterion: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  const ctx = await programContext(p.id);
  const c = ctx.criteria.find((x) => x.code === criterion && x.artifact === 'model_run');
  if (!c)
    throw new Error(`${program} has no criterion ${criterion} that rests on a run`);
  return comparePlanFrame(criterionComparePlan(await criterionRuns(p.id, c.id)), `${program} · ${criterion} runs`);
}

//name: addObservation
//description: Records an observation (F-26): a curated conclusion on the program, one of its series, campaigns or studies, or one of its compounds, as you, with optional evidence (a saved run by run_id or title, a snapshot id, a URL or a file path). Contributors and approvers
//input: string portfolio
//input: string program
//input: string subject {choices: ["program", "series", "compound", "campaign", "study"]}
//input: string code {optional: true; description: The series, campaign, study or compound code (none for the program)}
//input: string title
//input: string finding
//input: string evidenceKind {optional: true; choices: ["model_run", "snapshot", "url", "file"]}
//input: string evidenceTarget {optional: true; description: run_id or title, snapshot id, URL or file path}
//output: string report
export async function addObservation(portfolio: string, program: string, subject: string, code: string | undefined,
  title: string, finding: string, evidenceKind?: string, evidenceTarget?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const evidence = evidenceKind && evidenceTarget ? [await evidenceItemFor(p, evidenceKind as EvidenceKind, evidenceTarget)] : [];
  const o = await addObs(p.id, {subject: subject as ObservationSubject, code: code || null, title, finding, evidence});
  return `${p.code}: observation recorded on ${describeObservation(o)} (${o.id})`;
}

//name: programDiscussion
//description: A program's discussion (F-51): the platform chat on its workspace folder, as rows newest first: thread (General, or the criterion, run or document it is about), on, by, question, text, id
//input: string portfolio
//input: string program
//output: dataframe result
export async function programDiscussion(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  const d = await discussionOf(p);
  if (d.reason)
    throw new Error(d.reason);
  return discussionFrame(d.threads, `${p.code} discussion`);
}

//name: postToDiscussion
//description: Posts to a program's discussion as you (F-51), in general or about a criterion (code), a saved run (run_id or title) or a document (title). Anyone who may open the program's folder may post
//input: string portfolio
//input: string program
//input: string text
//input: string aboutKind {optional: true; choices: ["criterion", "run", "document"]}
//input: string about {optional: true; description: The criterion's code, the run's run_id or title, the document's title}
//output: string report
export async function postToDiscussion(portfolio: string, program: string, text: string, aboutKind?: string,
  about?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const anchor = aboutKind && about ? await resolveAnchor(p.id, aboutKind as AnchorKind, about) : null;
  const {post} = await postToProgram(p, text, {anchor});
  return `${p.code}: posted to the discussion${anchor ? ` about ${describeAnchor(anchor)}` : ''} (${post.id})`;
}

//name: askSme
//description: Asks an expert about a program (F-51): the question, with the portfolio, program, stage, next gate and readiness and the item it is about, goes into the program's discussion; the people the portfolio names (askSme; its approvers by default), or the logins given, are told in the platform's notifications where the server allows it. Nothing is sent outside the platform
//input: string portfolio
//input: string program
//input: string question
//input: string aboutKind {optional: true; choices: ["criterion", "run", "document"]}
//input: string about {optional: true}
//input: string logins {optional: true; description: Comma-separated logins to ask instead of the portfolio's choice}
//output: string report
export async function askSme(portfolio: string, program: string, question: string, aboutKind?: string, about?: string,
  logins?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const anchor = aboutKind && about ? await resolveAnchor(p.id, aboutKind as AnchorKind, about) : null;
  const to = (logins ?? '').split(',').map((x) => x.trim()).filter((x) => x);
  const r = await askExpert(p.id, question, await smeContextOf(p.id, anchor), to.length ? {logins: to} : undefined);
  return `${p.code}: ${r.line}`;
}

//name: retractObservation
//description: Retracts an observation: kept, marked retracted, with who, when and why. Its author or an approver
//input: string id
//input: string reason
//output: string report
export async function retractObservation(id: string, reason: string): Promise<string> {
  const o = await retractObs(id, reason);
  return `observation retracted: ${describeObservation(o)}`;
}

//name: observations
//description: A program's observations, current first, newest first: about, title, finding, evidence, status, by, on, id. Narrowed to a subject when given
//input: string portfolio
//input: string program
//input: string subject {optional: true; choices: ["program", "series", "compound", "campaign", "study"]}
//input: string code {optional: true}
//output: dataframe result
export async function observations(portfolio: string, program: string, subject?: string, code?: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  const list = observationsAbout(await observationsOf(p.id), (subject || null) as ObservationSubject | null, code || null);
  return observationsFrame(list, await userNames(list.map((o) => o.by)), `${p.code} observations`);
}

//name: createDesignSet
//description: Proposes structures for a series of a program (F-23): one row each of [proposals] (a Molecule column or one named smiles/structure; label/id/name; a score column), ranked by score. Contributors and approvers
//input: string portfolio
//input: string program
//input: string series
//input: string name
//input: string source {choices: ["enumerated", "generative", "manual"]}
//input: dataframe proposals
//input: string scoreColumn {optional: true; description: Default: a column named score}
//input: string scoreName {optional: true; description: What the score is, e.g. predicted pIC50}
//input: bool lowerIsBetter {optional: true}
//output: string report
export async function createDesignSet(portfolio: string, program: string, series: string, name: string, source: string,
  proposals: DG.DataFrame, scoreColumn?: string, scoreName?: string, lowerIsBetter?: boolean): Promise<string> {
  const p = await programByCode(portfolio, program);
  const cols = proposalColumns(proposals.columns.toList().map((c) => ({name: c.name, semType: c.semType})), scoreColumn || null);
  const rows = Array.from({length: proposals.rowCount}, (_, i) =>
    Object.fromEntries(proposals.columns.names().map((c) => [c, proposals.get(c, i)])));
  const r = await newDesignSet(p, {name, series, source: source as DesignSource, scoreName: scoreName || null,
    higherIsBetter: !lowerIsBetter, proposals: proposalsFromRows(rows, cols, `${series}-D`)});
  return `${p.code}: design set ${r.set.code} created: ${describeDesignSet(summarizeDesignSet(r.proposals))}`;
}

//name: designSets
//description: A program's design proposals, by set (newest first) and rank: set, series, label, structure, score, rank, chosen, compound
//input: string portfolio
//input: string program
//output: dataframe result
export async function designSets(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  return designSetsFrame(await designSetsOf(p.id), `${p.code} design sets`);
}

//name: chooseForSynthesis
//description: Marks proposals of a design set chosen for synthesis, or not, with who and when. Contributors and approvers
//input: string portfolio
//input: string program
//input: string designSet {description: The design set code}
//input: string labels {description: Proposal labels, comma-separated}
//input: bool chosen {optional: true; description: Default true. False unmarks}
//output: string report
export async function chooseForSynthesis(portfolio: string, program: string, designSet: string, labels: string,
  chosen?: boolean): Promise<string> {
  const p = await programByCode(portfolio, program);
  const set = (await designSetsOf(p.id)).find((x) => x.set.code === designSet);
  if (!set)
    throw new Error(`${p.code} has no design set ${designSet}`);
  const n = await choose(p, set.set.id, labels.split(',').map((x) => x.trim()).filter((x) => x), chosen !== false);
  return `${p.code} ${designSet}: ${n} proposal(s) ${chosen !== false ? 'chosen for synthesis' : 'no longer chosen'}`;
}

//name: createSelectionSet
//description: A hit list (F-23): compounds of a program picked from a campaign, with the rationale. Contributors and approvers
//input: string portfolio
//input: string program
//input: string name
//input: string campaign {optional: true}
//input: string rationale
//input: string compounds {description: Compound codes, comma-separated}
//output: string report
export async function createSelectionSet(portfolio: string, program: string, name: string, campaign: string | undefined,
  rationale: string, compounds: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const r = await newSelectionSet(p, {name, campaign: campaign || null, rationale,
    compounds: compounds.split(',').map((x) => x.trim()).filter((x) => x)});
  return `${p.code}: selection set ${r.code} created`;
}

//name: selectionSets
//description: A program's selection sets (hit lists), newest first: name, campaign, compounds, rationale, on
//input: string portfolio
//input: string program
//output: dataframe result
export async function selectionSets(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  return selectionSetsFrame(await selectionSetsOf(p.id), `${p.code} selection sets`);
}

//name: loadExampleDesignSets
//description: Loads the synthetic design sets and selection set the package ships into EXS-KIN1. EX-* portfolios only; idempotent
//output: string report
export async function loadExampleDesignSets(): Promise<string> {
  return loadDesignExamples(_package);
}

//name: dmtaStatus
//description: Where each compound of a program stands in the DMTA loop (F-24): designed, requested, made, tested, analyzed, with dates, design label, series, cycle days designed to tested, and the sources. Derived from design sets, the registry, the assay results, observations and the portfolio's dmta sources
//input: string portfolio
//input: string program
//output: dataframe result
export async function dmtaStatus(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  return dmtaFrame(await dmtaOf(await portfolioByCode(portfolio), p), `${p.code} DMTA status`);
}

//name: dmtaCycles
//description: Per series of a program: compounds reaching each DMTA stage and the median cycle time designed to tested, in days
//input: string portfolio
//input: string program
//output: dataframe result
export async function dmtaCycles(portfolio: string, program: string): Promise<DG.DataFrame> {
  const p = await programByCode(portfolio, program);
  return cyclesFrame(await dmtaOf(await portfolioByCode(portfolio), p), `${p.code} DMTA cycles`);
}

//name: showProvenance
//description: Opens the provenance of a portal object in a dialog: two indented trees, each node clickable
//input: string kind {choices: ["program", "decision", "evidence", "model_run", "snapshot", "document", "compound"]}
//input: string id
export async function showProvenance(kind: string, id: string): Promise<void> {
  installStyles();
  await openProvenanceDialog(kind, id);
}

//name: proposeDecision
//description: Proposes a gate decision in a portfolio that asks for approval of it (approval.decisions). It takes effect once another approver approves it with reviewRequest. Approvers only
//input: string portfolio
//input: string program
//input: string decision {choices: ["advance", "recycle", "hold", "terminate"]}
//input: string rationale
//input: string stage {optional: true; description: The stage code an advance or recycle goes to; default: the next or previous stage}
//output: string report
export async function proposeDecision(portfolio: string, program: string, decision: string, rationale: string,
  stage?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const to = stage ? (await portalDb.stages.query({filter: `stage_model_id = ${filterValue(p.stage_model_id)} and code = ${filterValue(stage)}`, limit: ROW_LIMIT}))[0] : null;
  if (stage && !to)
    throw new Error(`${stage} is not a stage of ${program}'s model`);
  return (await propose(p.id, decision as any, rationale, {toStageId: to?.id ?? null})).line;
}

//name: submitDocument
//description: Submits a program's document (by its file path) for approval, in a portfolio that asks for approval of documents: draft or rejected becomes submitted. Contributors and approvers
//input: string portfolio
//input: string program
//input: string path {description: The document's file path, as linked}
//input: string rationale {optional: true}
//output: string report
export async function submitDocument(portfolio: string, program: string, path: string, rationale?: string): Promise<string> {
  const p = await programByCode(portfolio, program);
  const doc = (await documentsOf(p.id)).find((d) => d.path === path.trim());
  if (!doc)
    throw new Error(`${program} has no document at ${path}`);
  const r = await submitDoc(p, doc.id, rationale || undefined);
  return [r.line, ...(r.notified ? [`  approvers: ${r.notified}`] : []), ...(r.unplaced ? [`  visibility: ${r.unplaced}`] : [])].join('\n');
}

//name: submitModelCard
//description: Submits a model card for approval, in a portfolio whose catalogue holds the model and that asks for approval of model cards: draft or rejected becomes in review. Contributors and approvers of that portfolio
//input: string portfolio
//input: string model {description: The model's catalogue code (Package:Function)}
//input: string version {optional: true; description: Default: the model's current version}
//input: string rationale {optional: true}
//output: string report
export async function submitModelCard(portfolio: string, model: string, version?: string, rationale?: string): Promise<string> {
  const r = await submitCard(await portfolioByCode(portfolio), model, version || null, rationale || undefined);
  return [r.line, ...(r.notified ? [`  approvers: ${r.notified}`] : []), ...(r.unplaced ? [`  visibility: ${r.unplaced}`] : [])].join('\n');
}

//name: reviewRequest
//description: Approves or rejects an approval request with a rationale. On the last approval it takes effect: the decision is recorded (evidence, workspace), the document or model card approved. The proposer may not review their own request. Approvers only
//input: string request {description: The request's id}
//input: string verdict {choices: ["approve", "reject"]}
//input: string rationale
//output: string report
export async function reviewRequest(request: string, verdict: string, rationale: string): Promise<string> {
  const r = await review(request, verdict as any, rationale);
  return [r.line, ...(r.notified ? [`  proposer: ${r.notified}`] : [])].join('\n');
}

//name: withdrawRequest
//description: Withdraws a pending approval request; only its proposer may. A document or model card goes back to draft
//input: string request {description: The request's id}
//output: string report
export async function withdrawRequest(request: string): Promise<string> {
  return `${(await withdraw(request)).title}: withdrawn`;
}

//name: approvalRequests
//description: The approval requests of a portfolio, or of one of its programs: kind, what, status, progress, who proposed it and when, the reviews, the result
//input: string portfolio
//input: string program {optional: true}
//input: string status {optional: true; choices: ["pending", "approved", "rejected", "withdrawn"]}
//output: dataframe result
export async function approvalRequests(portfolio: string, program?: string, status?: string): Promise<DG.DataFrame> {
  const pf = await portfolioByCode(portfolio);
  const rows: ApprovalRow[] = await ofMyPortfolios(program ? await requestsOf((await programByCode(portfolio, program)).id) :
    (await portalDb.approvals.query({filter: `portfolio_id = ${filterValue(pf.id)}`, limit: ROW_LIMIT}))
      .sort(newestFirstBy((x) => x.proposed_on ?? x.created_on)));
  const shown = rows.filter((r) => !status || r.status === status);
  return requestsFrame(shown, await userNames(shown.flatMap((r) => [r.proposed_by, ...reviewsOf(r.reviews).map((x) => x.by)])),
    `${program || portfolio} · Approval requests`);
}

//name: waitingForMe
//description: The approval requests the current user may review now, across portfolios, oldest first
//output: dataframe result
export async function waitingForMe(): Promise<DG.DataFrame> {
  const rows = await waitingOf();
  return requestsFrame(rows, await userNames(rows.flatMap((r) => [r.proposed_by, ...reviewsOf(r.reviews).map((x) => x.by)])),
    'Waiting for me');
}

//name: loadExampleDocuments
//description: Links the example documents the package ships to the EX-BIOC programs they name. EX-* portfolios only; idempotent
//output: string report
export async function loadExampleDocuments(): Promise<string> {
  return loadExamples(_package);
}

//name: seedExamples
//description: Seeds the example history into the EX-* portfolios: gate decisions with rationale and evidence, criteria met, observations, a follow, the example documents and two example discussion threads, dated in the past. Admins only; EX-* only; idempotent. A dry run (the default) only says what it would write
//input: bool dryRun = true
//output: string report
export async function seedExamples(dryRun: boolean = true): Promise<string> {
  return (await seedExampleHistory(_package, dryRun)).report;
}

//name: loadPortfolios
//description: Loads every portfolio file the package ships (files/portfolios/*.json): groups, stage models, portfolios, feeds. Idempotent.
//output: string report
export async function loadPortfolios(): Promise<string> {
  return loadAllPortfolios(_package);
}

//name: loadPortfolioConfigs
//description: Loads portfolios handed over by another package: a JSON array of portfolio files, checked as a set. Idempotent.
//input: string configs
//output: string report
export async function loadPortfolioConfigs(configs: string): Promise<string> {
  return loadPortfolioSet(JSON.parse(configs));
}

//name: securePortalTables
//description: Sets the server's table grants: Edit on program_role and parked_folder to each portfolio's approvers, on approval to its contributors, on sync_log to its approvers; View on parked_folder to the approvers; and All users' View, Edit and Delete taken from the tables that let every user read other portfolios' records (approval_request, sync_run, follow, compound, parked_folder). Loading a portfolio does the same for that portfolio. Needs Share on the tables (admins). Idempotent
//input: bool dryRun = true
//output: string report
export async function securePortalTables(dryRun: boolean = true): Promise<string> {
  const lines = await secureWriteTables(undefined, dryRun);
  return lines.length ? lines.join('\n') : 'table permissions already in place';
}

//name: moveToSecuredTables
//description: Moves one portfolio's records off the tables every user could read: approval requests into approval and sync runs into sync_log (same ids), each linked into its program folder or the portfolio's root so only the portfolio's people read it; compound values onto the portfolio's program_compound links, and compounds whose every link has them emptied. Follows are not copied (a copy would belong to whoever copied it). Admins. Idempotent
//input: string code {description: Portfolio code}
//input: bool dryRun = true
//output: string report
export async function moveToSecuredTables(code: string, dryRun: boolean = true): Promise<string> {
  const lines = await moveSecured(code, dryRun);
  return lines.length ? lines.join('\n') : `${code}: nothing to move`;
}

//name: moveRecordsOutOfRoot
//description: Moves the sync runs and model-card requests linked into a portfolio's workspace root into _Sync history and _Approvals under it, so the root lists the program folders only. Each is linked into its folder first, then its link in the root removed; the records stay. Dry run by default. Applying needs one of the portfolio's approvers. Idempotent
//input: string portfolio
//input: bool dryRun = true
//output: string report
export async function moveRecordsOutOfRoot(portfolio: string, dryRun: boolean = true): Promise<string> {
  return (await moveRecords(portfolio, dryRun)).join('\n');
}

//name: dataAccessAudit
//description: Read-only: for each data connection a portfolio's datasets and feeds use, whether each of its groups (viewers, contributors, approvers) has View and Query on it, directly or through a parent group, with the fix for each gap. Every portfolio when no code is given. Grants nothing
//input: string code {optional: true; description: Portfolio code}
//output: string report
export async function dataAccessAudit(code?: string): Promise<string> {
  return summarizeAccess(await dataAccessChecks(code || undefined));
}

//name: dataAccessTable
//description: The data access audit as a table: one row per portfolio, connection and group, with view, query, status and fix
//input: string code {optional: true; description: Portfolio code}
//output: dataframe result
export async function dataAccessTable(code?: string): Promise<DG.DataFrame> {
  return accessFrame(await dataAccessChecks(code || undefined));
}

const findFrame = (rows: {p: FindProgram, matched: string}[], name: string): DG.DataFrame => {
  const df = DG.DataFrame.fromColumns([
    textColumn('portfolio', rows.map((r) => r.p.portfolio)),
    textColumn('program', rows.map((r) => r.p.code)),
    textColumn('name', rows.map((r) => r.p.name)),
    textColumn('matched', rows.map((r) => r.matched)),
  ]);
  df.name = name;
  return df;
};

//name: findPrograms
//description: The programs a text finds across every portfolio the current user can see: by code, name, target or a registry alias (a value in the program's details), best match first
//input: string query
//output: dataframe programs
export async function findPrograms(query: string): Promise<DG.DataFrame> {
  const hits = find(query, await findablePrograms());
  return findFrame(hits.map((h) => ({p: h.program, matched: matchedBy(h) ?? h.term.field})), `Find ${query}`);
}

//name: recentPrograms
//description: The programs the current user opened last in the portal that they can still see, newest first
//output: dataframe programs
export async function recentPrograms(): Promise<DG.DataFrame> {
  return findFrame(resolveRecents(savedRecents(), await findablePrograms()).map((p) => ({p, matched: ''})), 'Recent programs');
}

//name: linkProgramRows
//description: Makes a portfolio's program rows visible to whoever may view their workspace folders (the portfolio's groups and the program roles) by linking each row into its folder. View only. Creates no folder. Idempotent
//input: string portfolio
//input: bool dryRun = true
//output: string report
export async function linkProgramRows(portfolio: string, dryRun: boolean = true): Promise<string> {
  const lines = await linkRows(portfolio, dryRun);
  return lines.length ? lines.join('\n') : 'no programs';
}

//name: loadModelCatalogue
//description: Loads a model catalogue handed over without a portfolio file: unit operations, run types, models and their cards, as JSON. Checked before loading; idempotent
//input: string catalogue {description: JSON object with unitOperations, runTypes, models}
//input: string portfolio {optional: true; description: The portfolio the models belong to}
//output: string report
export async function loadModelCatalogue(catalogue: string, portfolio?: string): Promise<string> {
  return loadCatalogues([{portfolio: portfolio || null, catalogue: JSON.parse(catalogue)}]);
}

//name: modelCard
//description: The card of a catalogued model: the question it answers, intended use, inputs with their validity range, outputs, applicability, performance, lineage, approval
//input: string model {description: The model's function nqName, e.g. ProjectPortal:FedBatchTiter}
//output: widget result
export async function modelCard(model: string): Promise<DG.Widget> {
  installStyles();
  return modelCardWidget(model);
}

//name: workspaceModelCard
//friendlyName: Model card
//description: The model card beside a model's preview in the Workspace tab on Home: approval, portfolios, question, inputs and outputs. Nothing for a function the portal's catalogue does not list
//input: func model
//output: widget result
//meta.role: workspaceModelPanel
export async function workspaceModelCard(model: DG.Func): Promise<DG.Widget | null> {
  return workspaceModelCardWidget(model);
}

//name: modelCardPanel
//friendlyName: Model card
//description: The Model card of a function: its card, labelled with the user's portfolios whose catalogue lists it. Nothing for a function those catalogues do not list. The context panel shows it as the Model card pane of a selected model (Model Hub, Browse)
//input: func model
//output: widget result
//meta.role: panel
export async function modelCardPanel(model: DG.Func): Promise<DG.Widget | null> {
  return modelCardPanelWidget(model);
}

/** Programs first, so work units find the programs they hang under. A real sync is recorded in the
 * portfolio's sync history; a dry run is not. The report is the same either way. */
async function syncOne(code: string, dryRun = false): Promise<string> {
  if (!dryRun) {
    const r = await recordedSync(code, 'manual');
    if (r.error)
      throw new Error(r.error);
    return r.report;
  }
  const plans = await syncPortfolioPlans(code, dryRun);
  const programs = plans.map(summarizePlan);
  const units = (await syncWorkUnits(code, dryRun)).map(summarizeWorkUnitPlan);
  const compounds = await syncCompounds(code, dryRun);
  // the folders a real sync would make, for the programs that have none before it
  const fresh = plans.flatMap((p) => p.inserts).filter((v) => !noNewFoldersReason(v)).map((v) => v.code);
  const folders = summarizeFolderGrowth(await buildMissingFolders(code, true, fresh));
  // what data would advance, judged on the programs as they stand before the sync
  const advance = autoAdvanceOf(await portfolioByCode(code)) ? [summarizeAutoAdvanceResult(await autoAdvance(code, true))] : [];
  return [...programs, ...units, ...(compounds ? [summarizeCompoundPlan(compounds)] : []), ...(folders ? [folders] : []),
    ...advance].join('\n');
}

//name: syncPortfolio
//description: Runs one portfolio's registry feeds, then its work-unit and compound feeds, and applies (or with dryRun only plans) the result. A real sync then builds the workspace folders of programs that have none (not of ended or parked ones); a dry run lists them
//input: string code
//input: bool dryRun = false
//output: string report
export async function syncPortfolio(code: string, dryRun: boolean = false): Promise<string> {
  return syncOne(code, dryRun);
}

//name: syncAllPortfolios
//description: Syncs every loaded portfolio, discovery before development so lineage resolves in one pass
//output: string report
export async function syncAllPortfolios(): Promise<string> {
  const portfolios = await portalDb.portfolios.query({limit: ROW_LIMIT});
  const order = (s?: string) => s === 'discovery' ? 0 : s === 'development' ? 1 : 2;
  const out: string[] = [];
  for (const p of portfolios.sort((a, b) => order(a.setting) - order(b.setting)))
    out.push(await syncOne(p.code));
  return out.join('\n');
}

//name: syncDuePortfolios
//description: Syncs the portfolios whose schedule is due, and only those: a portfolio without a schedule is never touched. Each run is recorded in the portfolio's sync history; owners are notified of a run that changed something
//input: string codes {optional: true; description: Comma-separated portfolio codes to limit the check to}
//output: string report
export async function syncDuePortfolios(codes?: string): Promise<string> {
  const list = (codes ?? '').split(',').map((c) => c.trim()).filter((c) => c.length > 0);
  return summarizeDue(await syncDue(list.length ? list : undefined));
}

//name: autoAdvancePortfolio
//description: Stage advance driven by data, for a portfolio that declares autoAdvance: each active program whose next gate is met from its data alone moves one stage, never past autoAdvance.through, logged with source data and the evidence. With dryRun (the default) it only says what would move and why the rest stay
//input: string code
//input: bool dryRun = true
//output: string report
export async function autoAdvancePortfolio(code: string, dryRun: boolean = true): Promise<string> {
  return summarizeAutoAdvanceResult(await autoAdvance(code, dryRun));
}

//name: exportPortfolio
//description: A loaded portfolio as a portfolio file (JSON): stage model, criteria, feeds, groups by name, datasets, rules, schedule, its model catalogue. Loading it gives the same portfolio
//input: string code
//output: string json
export async function exportPortfolio(code: string): Promise<string> {
  return JSON.stringify(await exportPortfolioConfig(code), null, 2);
}

//name: diffPortfolios
//description: What changes from portfolio a to portfolio b, in readable lines: vocabulary (stages, criteria, work units, catalogue codes), workspace tree (workstreams, stage folders) and other fields. Each is a loaded portfolio's code or a portfolio file's JSON
//input: string a {description: Portfolio code (loaded) or portfolio JSON}
//input: string b {description: Portfolio code (loaded) or portfolio JSON}
//output: string report
export async function diffPortfolios(a: string, b: string): Promise<string> {
  const r = await diffPortfolioRefs(a, b);
  return summarizeDiff(r.title, r.diff);
}

//name: exportProfile
//description: A profile bundle for moving portfolios to another instance: the loaded portfolios named (default: all) as portfolio files in one JSON
//input: string codes {optional: true; description: Comma-separated portfolio codes}
//output: string json
export async function exportProfile(codes?: string): Promise<string> {
  const list = (codes ?? '').split(',').map((c) => c.trim()).filter((c) => c.length > 0);
  return JSON.stringify(await exportProfileBundle(list.length ? list : undefined), null, 2);
}

//name: importProfile
//description: Imports a profile bundle: checks its portfolios as a set (refused whole when they do not validate), diffs each against the loaded one, and unless dryRun (the default) loads them. Loading adds and updates; it never removes
//input: string bundle {description: The JSON exportProfile returns}
//input: bool dryRun = true
//output: string report
export async function importProfile(bundle: string, dryRun: boolean = true): Promise<string> {
  return (await importProfileBundle(bundle, dryRun)).report;
}

//name: syncHistory
//description: A portfolio's sync runs, newest first: when, by schedule or by hand, who, counts of created, updated and flagged rows, unresolved spellings, error
//input: string portfolio
//output: dataframe result
export async function syncHistory(portfolio: string): Promise<DG.DataFrame> {
  return syncHistoryFrame(await portfolioByCode(portfolio));
}

//name: programReadiness
//description: Readiness of one program for a stage (default: its next stage): every criterion met, waived or open, data conditions evaluated over the program dataset, and the gap list
//input: string portfolio
//input: string program
//input: string stage {optional: true}
//output: string report
export async function programReadiness(portfolio: string, program: string, stage?: string): Promise<string> {
  const r = await readinessOf(portfolio, program, stage || undefined);
  return summarizeReadiness(r.program.code, r.stage.name, r.readiness);
}

//name: setCompoundRole
//description: Changes a compound's role in a program (explored, hit, lead, backup, candidate, tool, discontinued) and logs who, when and why
//input: string portfolio
//input: string program
//input: string compound {description: Registration code}
//input: string role
//input: string rationale
//output: string report
export async function setCompoundRole(portfolio: string, program: string, compound: string, role: string,
  rationale: string): Promise<string> {
  const pf = await portfolioByCode(portfolio);
  const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(pf.id)} and code = ${filterValue(program)}`, limit: ROW_LIMIT}))[0];
  if (!p)
    throw new Error(`${portfolio} has no program ${program}, or it is not shared with you`);
  // the program's own link first: the compound table is closed to most users
  const linked = await programCompound(p.id, compound);
  const c = linked ? [linked.compound] : await portalDb.compounds.query({filter: `code = ${filterValue(compound)}`, limit: ROW_LIMIT}).catch(() => []);
  if (c.length !== 1)
    throw new Error(c.length ? `${compound} is registered in ${c.length} registries` : `no compound ${compound} is known to the portal`);
  await setRole(p.id, c[0].id, role as CompoundRole, rationale);
  return `${program}/${compound}: ${role}`;
}

//name: portfolioRoutingAudit
//description: For each source system, the keys no portfolio feed takes and the keys two feeds both take
//output: string report
export async function portfolioRoutingAudit(): Promise<string> {
  return summarizeRouting(await routingAudit(await readSources(_package)));
}

//name: routingAuditOf
//description: The routing audit for source systems handed over by another package: a JSON array of sources
//input: string sources
//output: string report
export async function routingAuditOf(sources: string): Promise<string> {
  return summarizeRouting(await routingAudit(JSON.parse(sources)));
}

//name: buildWorkspace
//description: Creates the portfolio's Spaces tree (root, a folder per program, its workstreams) and shares the root with the portfolio's groups. Additive and idempotent; dryRun previews what would be created
//input: string code
//input: bool dryRun = false
//input: string program {optional: true; description: Build only this program's folders (its code)}
//output: string report
export async function buildWorkspace(code: string, dryRun: boolean = false, program?: string): Promise<string> {
  let only: string | undefined;
  if (program) {
    const pf = await portfolioByCode(code);
    const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(pf.id)} and code = ${filterValue(program)}`, limit: ROW_LIMIT}))[0];
    if (!p)
      throw new Error(`${code} has no program ${program}, or it is not shared with you`);
    only = p.id;
  }
  return summarizeWorkspace(await buildPortfolioWorkspace(code, dryRun, only));
}

//name: workspaceOrphans
//description: Reports, and never touches, the folders a prune-or-park pass would look at: root folders no program of the portfolio claims, and folders of programs that have ended or left their feed
//input: string code
//output: string report
export async function workspaceOrphans(code: string): Promise<string> {
  return summarizeWorkspace(await orphansOf(code));
}

//name: parkOrphans
//description: Moves the folders workspaceOrphans lists (or those named) into _Parked/<date> under the portfolio's root and records each move; never deletes. Dry run by default. Applying needs one of the portfolio's approvers
//input: string code
//input: bool dryRun = true
//input: string folders {optional: true; description: Comma-separated paths or space ids of the orphans to park; default: all}
//output: string report
export async function parkOrphans(code: string, dryRun: boolean = true, folders?: string): Promise<string> {
  const only = (folders ?? '').split(',').map((f) => f.trim()).filter((f) => f.length > 0);
  return summarizeParkResult(await park(code, dryRun, only.length ? only : undefined));
}

//name: unparkFolder
//description: Moves a parked folder back to where it was, by its parked-folder record id or its space id. Approvers only
//input: string id
//output: string report
export async function unparkFolder(id: string): Promise<string> {
  const r = await unpark(id);
  return [r.line, ...r.notes.map((n) => `  ${n}`)].join('\n');
}

//name: parkedFolders
//description: A portfolio's parked folders: what, why, from where, to where, when, and whether it was restored
//input: string code
//input: bool all = false {description: Also the folders restored since}
//output: dataframe result
export async function parkedFolders(code: string, all: boolean = false): Promise<DG.DataFrame> {
  return parkedFrame(await parkedFoldersOf(await portfolioByCode(code), all), `${code} · Parked folders`);
}

//name: pinProgram
//description: Adds the program to your favorites (a platform favorite on its record), so it shows in your Workspace on Home and its tile in My programs. With a group, adds it for everyone in that group: the portfolio's approvers only
//input: string portfolio
//input: string program
//input: string group {optional: true; description: Group name; default: you}
//output: string report
export async function pinProgram(portfolio: string, program: string, group?: string): Promise<string> {
  return pinByCode(portfolio, program, true, group || null);
}

//name: unpinProgram
//description: Takes the program out of your favorites, or a group's (its record, and its folder when pinned before 0.3.0), and its tile off Home
//input: string portfolio
//input: string program
//input: string group {optional: true; description: Group name; default: you}
//output: string report
export async function unpinProgram(portfolio: string, program: string, group?: string): Promise<string> {
  return pinByCode(portfolio, program, false, group || null);
}

//name: pinnedPrograms
//description: The programs in your favorites (yours or your groups') that you can see, across portfolios: stage on its track, status, readiness for the next gate, attention, new data
//output: dataframe result
export async function pinnedPrograms(): Promise<DG.DataFrame> {
  return pinnedFrame(await pinnedTiles());
}

//name: myProgramsWidget
//friendlyName: My programs
//description: Home widget: a tile per pinned program, with its stage, readiness for the next gate, attention and new data; a click opens its page
//output: widget result
//meta.role: dashboard
//meta.order: 2
export async function myProgramsWidget(): Promise<DG.Widget | null> {
  // a user who sees no portfolio gets no widget
  if (!(await portalDb.portfolios.count({}).catch(() => 0)))
    return null;
  installStyles();
  return new MyProgramsWidget();
}

//name: saveModelRun
//description: Saves a model run to a program: annotates the call (portfolio, program, stage at run, unit operation, run type, purpose), saves it, links it into the program's workspace (the current stage's folder when there is one) and records it as a model run a gate criterion can select
//input: funccall call {description: A call that has run}
//input: string portfolio
//input: string program
//input: string unitOperation {optional: true; description: Unit operation code}
//input: string runType {optional: true; description: Run type code}
//input: string purpose {optional: true; choices: ["screening", "characterization", "scale-up", "troubleshooting", "verification"]}
//input: string title {optional: true}
//input: string folder {optional: true; description: Folder under the program's folder, '/'-separated, by default the current stage's}
//output: string report
export async function saveModelRun(call: DG.FuncCall, portfolio: string, program: string, unitOperation?: string,
  runType?: string, purpose?: string, title?: string, folder?: string): Promise<string> {
  return summarizeSavedRun(await saveRun(call, portfolio, program, {unitOperation: unitOperation || undefined,
    runType: runType || undefined, purpose: (purpose || undefined) as RunPurpose | undefined,
    title: title || undefined, folder: folder || undefined}));
}

//name: buildAllWorkspaces
//description: Builds the workspace of every loaded portfolio that declares a workspaceRoot, one after another
//output: string report
export async function buildAllWorkspaces(): Promise<string> {
  const portfolios = (await portalDb.portfolios.query({limit: ROW_LIMIT})).filter((p) => p.workspace_root);
  const out: string[] = [];
  for (const p of portfolios.sort((a, b) => a.code.localeCompare(b.code)))
    out.push(summarizeWorkspace(await buildPortfolioWorkspace(p.code)));
  return out.join('\n');
}
