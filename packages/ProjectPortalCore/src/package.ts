/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {loadAllPortfolios, loadPortfolioSet, programReadiness as readinessOf, readSources, routingAudit, summarizePlan, summarizeRouting, syncPortfolio as syncPortfolioPlans, syncWorkUnits, syncCompounds, setCompoundRole as setRole, portfolioByCode} from './portal';
import {portalDb} from './generated/db';
import {buildPortfolioWorkspace, summarizeWorkspace, workspaceOrphans as orphansOf} from './workspace';
import {summarizeReadiness} from './model/readiness';
import {summarizeWorkUnitPlan} from './model/work-units';
import {CompoundRole, summarizeCompoundPlan} from './model/compounds';
import {portalAppView} from './views/portal-view';
import {installStyles} from './views/styles';
import {portfolioTableView, TABLE_SEGMENT} from './views/portfolio-table-view';
import {portfolioTables} from './portfolio-table';
import {saveModelRun as saveRun, summarizeSavedRun} from './model-runs';
import {RunPurpose} from './model/model-runs';
import {loadCatalogues} from './catalogue';
import {modelCardWidget} from './views/model-card';
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
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/u2/src/dg/domain/styles.js';

export const _package = new DG.Package();

/** After a run of a catalogued model: a balloon that opens the save-to-program dialog for it. */
function offerSave(call: DG.FuncCall): void {
  installStyles();
  const link = span('Save to program…', 'pp-link');
  link.dataset.u2Name = 'save-run-link';
  link.onclick = () => void openSaveRunDialog({call}).catch((e) => grok.shell.error(e?.message ?? String(e)));
  const box = document.createElement('div');
  box.append(span(`${call.func.friendlyName ?? call.func.name} ran. `), link);
  grok.shell.info(box, {timeout: 20});
}

//name: portalInit
//tags: init
//description: Once the package loads: offers to save each run of a catalogued model to a program, and adds Save to program to a run's context menu
export async function portalInit(): Promise<void> {
  trackModelRuns(offerSave);
  grok.events.onContextMenu.subscribe((e: any) => {
    const item = e?.args?.item;
    if (item instanceof DG.FuncCall && e.args.menu)
      e.args.menu.item('Save to program…', () => void openSaveRunDialog({call: item}));
  });
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
  return portalAppView(segments.pop());
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

/** A program of a portfolio by code, under row security. */
async function programByCode(portfolio: string, program: string) {
  const pf = await portfolioByCode(portfolio);
  const p = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${program}"`}))[0];
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
//input: string model {description: The model's function nqName, e.g. ProjectPortalCore:FedBatchTiter}
//output: widget result
export async function modelCard(model: string): Promise<DG.Widget> {
  installStyles();
  return modelCardWidget(model);
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
  const programs = (await syncPortfolioPlans(code, dryRun)).map(summarizePlan);
  const units = (await syncWorkUnits(code, dryRun)).map(summarizeWorkUnitPlan);
  const compounds = await syncCompounds(code, dryRun);
  // what data would advance, judged on the programs as they stand before the sync
  const advance = autoAdvanceOf(await portfolioByCode(code)) ? [summarizeAutoAdvanceResult(await autoAdvance(code, true))] : [];
  return [...programs, ...units, ...(compounds ? [summarizeCompoundPlan(compounds)] : []), ...advance].join('\n');
}

//name: syncPortfolio
//description: Runs one portfolio's registry feeds, then its work-unit and compound feeds, and applies (or with dryRun only plans) the result
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
  const portfolios = await portalDb.portfolios.query({});
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
  const p = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${program}"`}))[0];
  if (!p)
    throw new Error(`${portfolio} has no program ${program}, or it is not shared with you`);
  const c = (await portalDb.compounds.query({filter: `code = "${compound}"`}));
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
    const p = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${program}"`}))[0];
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
  const portfolios = (await portalDb.portfolios.query({})).filter((p) => p.workspace_root);
  const out: string[] = [];
  for (const p of portfolios.sort((a, b) => a.code.localeCompare(b.code)))
    out.push(summarizeWorkspace(await buildPortfolioWorkspace(p.code)));
  return out.join('\n');
}
