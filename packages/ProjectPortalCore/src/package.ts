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
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/u2/src/dg/domain/styles.js';

export const _package = new DG.Package();

//name: Project Portal
//description: Several portfolios on one instance, each on its own stage-gate model and registry
//tags: app
//input: string path {meta.url: true; optional: true}
//output: view result
export async function projectPortalApp(path?: string): Promise<DG.ViewBase> {
  installStyles();
  const code = path?.split('/').filter((s) => s.length > 0).pop();
  return portalAppView(code);
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

/** Programs first, so work units find the programs they hang under. */
async function syncOne(code: string, dryRun = false): Promise<string> {
  const programs = (await syncPortfolioPlans(code, dryRun)).map(summarizePlan);
  const units = (await syncWorkUnits(code, dryRun)).map(summarizeWorkUnitPlan);
  const compounds = await syncCompounds(code, dryRun);
  return [...programs, ...units, ...(compounds ? [summarizeCompoundPlan(compounds)] : [])].join('\n');
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
//output: string report
export async function buildWorkspace(code: string, dryRun: boolean = false): Promise<string> {
  return summarizeWorkspace(await buildPortfolioWorkspace(code, dryRun));
}

//name: workspaceOrphans
//description: Reports, and never touches, the folders a prune-or-park pass would look at: root folders no program of the portfolio claims, and folders of programs that have ended or left their feed
//input: string code
//output: string report
export async function workspaceOrphans(code: string): Promise<string> {
  return summarizeWorkspace(await orphansOf(code));
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
