/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {loadAllPortfolios, loadPortfolioSet, readSources, routingAudit, summarizePlan, summarizeRouting, syncPortfolio as syncPortfolioPlans} from './portal';
import {portalDb} from './generated/db';
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

//name: syncPortfolio
//description: Runs one portfolio's registry feeds and applies (or with dryRun only plans) the result
//input: string code
//input: bool dryRun = false
//output: string report
export async function syncPortfolio(code: string, dryRun: boolean = false): Promise<string> {
  return (await syncPortfolioPlans(code, dryRun)).map(summarizePlan).join('\n');
}

//name: syncAllPortfolios
//description: Syncs every loaded portfolio, discovery before development so lineage resolves in one pass
//output: string report
export async function syncAllPortfolios(): Promise<string> {
  const portfolios = await portalDb.portfolios.query({});
  const order = (s?: string) => s === 'discovery' ? 0 : s === 'development' ? 1 : 2;
  const out: string[] = [];
  for (const p of portfolios.sort((a, b) => order(a.setting) - order(b.setting)))
    out.push((await syncPortfolioPlans(p.code)).map(summarizePlan).join('\n'));
  return out.join('\n');
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
