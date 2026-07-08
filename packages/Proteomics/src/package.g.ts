import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
//meta.role: init
export async function initProteomics() : Promise<void> {
  await PackageFunctions.initProteomics();
}

//tags: autostart
//meta.autostartImmediate: true
export async function buildProteomicsTopMenu() : Promise<void> {
  await PackageFunctions.buildProteomicsTopMenu();
}

//top-menu: Proteomics | Import | Spectronaut Candidates...
export async function importSpectronautCandidates() : Promise<void> {
  await PackageFunctions.importSpectronautCandidates();
}

//top-menu: Proteomics | Import | Spectronaut Report...
export async function importSpectronaut() : Promise<void> {
  await PackageFunctions.importSpectronaut();
}

//top-menu: Proteomics | Import | MaxQuant...
export async function importMaxQuant() : Promise<void> {
  await PackageFunctions.importMaxQuant();
}

//top-menu: Proteomics | Import | FragPipe...
export async function importFragPipe() : Promise<void> {
  await PackageFunctions.importFragPipe();
}

//top-menu: Proteomics | Import | Generic Matrix...
export async function importGenericMatrix() : Promise<void> {
  await PackageFunctions.importGenericMatrix();
}

//name: annotateExperiment
export async function annotateExperiment() : Promise<void> {
  await PackageFunctions.annotateExperiment();
}

//name: setLog2Scale
export async function setLog2Scale() : Promise<void> {
  await PackageFunctions.setLog2Scale();
}

//name: normalizeProteomics
export async function normalizeProteomics() : Promise<void> {
  await PackageFunctions.normalizeProteomics();
}

//name: imputeMissingValues
export async function imputeMissingValues() : Promise<void> {
  await PackageFunctions.imputeMissingValues();
}

//name: differentialExpression
export async function differentialExpression() : Promise<void> {
  await PackageFunctions.differentialExpression();
}

//name: computeSpcStatus
export async function computeSpcStatus() : Promise<void> {
  await PackageFunctions.computeSpcStatus();
}

//name: showVolcanoPlot
export async function showVolcanoPlot() : Promise<void> {
  await PackageFunctions.showVolcanoPlot();
}

//name: showHeatmap
export async function showHeatmap() : Promise<void> {
  await PackageFunctions.showHeatmap();
}

//name: showPcaPlot
export async function showPcaPlot() : Promise<void> {
  await PackageFunctions.showPcaPlot();
}

//name: showGroupMeanCorrelation
export async function showGroupMeanCorrelation() : Promise<void> {
  await PackageFunctions.showGroupMeanCorrelation();
}

//name: showQcDashboard
export async function showQcDashboard() : Promise<void> {
  await PackageFunctions.showQcDashboard();
}

//name: showSpcDashboard
export async function showSpcDashboard() : Promise<void> {
  await PackageFunctions.showSpcDashboard();
}

//name: showAllVisualizations
export async function showAllVisualizations() : Promise<void> {
  await PackageFunctions.showAllVisualizations();
}

//name: enrichmentAnalysis
export async function enrichmentAnalysis() : Promise<void> {
  await PackageFunctions.enrichmentAnalysis();
}

//name: exportEnrichmentInputs
export async function exportEnrichmentInputs() : Promise<void> {
  await PackageFunctions.exportEnrichmentInputs();
}

//name: rankAbundance
export async function rankAbundance() : Promise<void> {
  await PackageFunctions.rankAbundance();
}

//name: enrichmentCharts
export async function enrichmentCharts() : Promise<void> {
  await PackageFunctions.enrichmentCharts();
}

//name: Proteomics Demo
//meta.demoPath: Proteomics | Differential Expression
//meta.isDemoDashboard: true
export async function proteomicsDemo() : Promise<void> {
  await PackageFunctions.proteomicsDemo();
}

//name: Proteomics Enrichment Demo
//meta.demoPath: Proteomics | Enrichment Analysis
//meta.isDemoDashboard: true
export async function proteomicsEnrichmentDemo() : Promise<void> {
  await PackageFunctions.proteomicsEnrichmentDemo();
}

//name: Proteomics | UniProt
//description: UniProt protein details
//input: string proteinId { semType: Proteomics-ProteinId }
//output: widget result
//meta.role: widgets,panel
export function uniprotPanelWidget(proteinId: string) : any {
  return PackageFunctions.uniprotPanelWidget(proteinId);
}

//name: shareAnalysisForReview
export async function shareAnalysisForReview() : Promise<void> {
  await PackageFunctions.shareAnalysisForReview();
}

//name: Proteomics | Shared Analysis
//description: Audit context for a shared analysis snapshot
//input: string proteinId { semType: Proteomics-ProteinId }
//output: widget result
//meta.role: widgets,panel
export function publishedAnalysisPanelWidget(proteinId: string) : any {
  return PackageFunctions.publishedAnalysisPanelWidget(proteinId);
}

//tags: autostart
//meta.autostartImmediate: true
export async function recoverPublishedProjectsOnStartup() : Promise<void> {
  await PackageFunctions.recoverPublishedProjectsOnStartup();
}
