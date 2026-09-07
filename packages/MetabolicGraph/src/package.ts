/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {parsePath, loadStateProxy, analysisExists, readAnalysis, handleReactionDataUpload} from './utils';
import map from './maps/E_coli_Core_metabolism_map.json';
import model from './maps/E_coli_core_cobra.json';
import type {MapData, CobraModelData} from '../escher_src/src/ts/types';
import {MetabolicGraphView, createEscherBuilder} from './view';
import * as aiTools from './ai-tools';

export const _package = new DG.Package();

declare global {
  // eslint-disable-next-line no-unused-vars
  interface Window {
    escher: any;
    builder: any
  }
}

//name: MetabolicGraph
//meta.role: app
//meta.icon: files/icons/metabolic.png
//description: Metabolic graph application
//input: string path {meta.url: true; optional: true}
//input: string filter {optional: true}
//output: view v
//meta.browsePath: Misc
export function metabolicGraphApp(path?: string, filter?: string): DG.ViewBase {
  const view = new MetabolicGraphView();
  setTimeout(() => {
    const b = createEscherBuilder(view, map as unknown as MapData, model as unknown as CobraModelData);
    handleReactionDataUpload(view, b);
    setTimeout(async () => {
      if (!b.map)
        return;
      const pathParsed = parsePath(path);
      if (pathParsed) {
        if (!await analysisExists(pathParsed)) {
          grok.shell.error(`File ${pathParsed}.json not found`);
          return;
        }
        loadStateProxy(b, await readAnalysis(pathParsed), pathParsed);
        const curHref = window.location.href;
        if (curHref && !curHref.endsWith(`/${pathParsed}`)) {
          // @ts-ignore
          if (history.replaceState) {
            const title = document.title;
            const obj = {Title: title, Url: `${curHref}/${pathParsed}`};
            history.replaceState(obj, obj.Title, obj.Url);
          }
        }
      }
    }, 500);
  }, 500);
  return view;
}

//name: EscherFileViewer
//meta.role: fileViewer
//meta.fileViewer: json
//meta.fileViewerCheck: Metabolicgraph:escherFileViewerCheck
//input: file file
//output: view v
export async function escherFileViewer(file: DG.FileInfo) {
  const view = new MetabolicGraphView();
  const mapJSON = await file.readAsString();
  const mapData = JSON.parse(mapJSON) as MapData;
  setTimeout(() => createEscherBuilder(view, mapData, null), 500);
  return view;
}

//name: escherFileViewerCheck
//input: string content
//output: bool result
export async function escherFileViewerCheck(content: string) {
  return (content?.length ?? 1e12) < 1e7 && !!content?.startsWith('[') && !!content?.includes('"https://escher.github.io/escher/jsonschema/1-0-0#"');
}

// 'metabolicViewFunction'-tagged: the AI assistant acts on the open Metabolic Graph view through these
// (collected via MetabolicGraphView.getFunctions). Bodies live in ai-tools.ts.

//name: getMetabolicState
//description: Current state of the metabolic graph view - model and map summary, FBA objectives, loaded flux data, selected metabolites, active time-course and analysis name. Call this first to understand what is loaded
//tags: metabolicViewFunction
//input: view view
//output: dynamic result
export function getMetabolicState(view: DG.ViewBase) {
  return aiTools.getMetabolicState(view);
}

//name: findMetabolicEntities
//description: Search the model reactions, metabolites and genes by id or name substring. Returns brief entries with bounds, objective and flux for reactions
//tags: metabolicViewFunction
//input: view view
//input: string query { description: Substring matched against ids and names. Pass an empty string to list everything }
//input: string kind { optional: true; choices: ["reaction", "metabolite", "gene"] }
//input: int limit { optional: true; description: Max matches returned, default 15 }
//output: dynamic result
export function findMetabolicEntities(view: DG.ViewBase, query: string, kind?: string, limit?: number) {
  return aiTools.findMetabolicEntities(view, query, kind, limit);
}

//name: getMetabolicEntityDetails
//description: Full detail of one reaction (stoichiometry, bounds, gene rule, flux), metabolite (formula, compartment, its reactions) or gene by id
//tags: metabolicViewFunction
//input: view view
//input: string entityId
//output: dynamic result
export function getMetabolicEntityDetails(view: DG.ViewBase, entityId: string) {
  return aiTools.getMetabolicEntityDetails(view, entityId);
}

//name: setReactionBounds
//description: Set flux bounds of model reactions. Pass a map of reaction id to a two-number list of lower and upper bound. Set both to 0 to knock a reaction out. Affects subsequent FBA and sampling
//tags: metabolicViewFunction
//input: view view
//input: map bounds
//output: dynamic result
export function setReactionBounds(view: DG.ViewBase, bounds: aiTools.BoundsMap) {
  return aiTools.setReactionBounds(view, bounds);
}

//name: runFba
//description: Run flux balance analysis and color the map by the resulting fluxes. Optionally pass objectives as a map of reaction id to maximize or minimize, replacing the current objectives. Returns the objective value and fluxes
//tags: metabolicViewFunction
//input: view view
//input: map objectives { optional: true }
//output: dynamic result
export async function runFba(view: DG.ViewBase, objectives?: aiTools.ObjectivesMap) {
  return aiTools.runFba(view, objectives);
}

//name: sampleReactionFluxes
//description: Sample the feasible flux space of the model. Colors the map by average flux and puts a flux histogram into every reaction tooltip. Returns the largest average fluxes
//tags: metabolicViewFunction
//input: view view
//input: int samples { optional: true; description: Number of samples, default 10000 }
//input: int bins { optional: true; description: Histogram bins in reaction tooltips, default 20 }
//input: int thinning { optional: true; description: Sampler thinning interval, default 10 }
//input: bool addDataFrame { optional: true; description: Also open the raw samples as a table view }
//output: dynamic result
export async function sampleReactionFluxes(view: DG.ViewBase, samples?: number, bins?: number,
  thinning?: number, addDataFrame?: boolean) {
  return aiTools.sampleReactionFluxes(view, samples, bins, thinning, addDataFrame);
}

//name: runTimeCourse
//description: Time-course simulation - interpolates reaction bounds from start to end over the steps, samples fluxes at every step, then shows an animation slider on the map and a table of average flux per step. Pass startBounds and endBounds as maps of reaction id to a two-number list of lower and upper bound, both listing the same reactions
//tags: metabolicViewFunction
//input: view view
//input: map startBounds
//input: map endBounds
//input: int steps { optional: true; description: Interpolation steps including both endpoints, default 10 }
//input: int samples { optional: true; description: Samples per step, default 10000 }
//output: dynamic result
export async function runTimeCourse(view: DG.ViewBase, startBounds: aiTools.BoundsMap, endBounds: aiTools.BoundsMap,
  steps?: number, samples?: number) {
  return aiTools.runTimeCourse(view, startBounds, endBounds, steps, samples);
}

//name: setReactionFluxes
//description: Color the map with custom per-reaction values. Pass values as a map of reaction id to number. Omit values to clear the current reaction data
//tags: metabolicViewFunction
//input: view view
//input: map values { optional: true }
//input: string source { optional: true; description: Label shown as the data source }
//output: dynamic result
export function setReactionFluxes(view: DG.ViewBase, values?: aiTools.FluxMap, source?: string) {
  return aiTools.setReactionFluxes(view, values, source);
}

//name: showMetabolicEntity
//description: Zoom the map to a reaction or metabolite and highlight it so the user sees it
//tags: metabolicViewFunction
//input: view view
//input: string entityId
//output: dynamic result
export function showMetabolicEntity(view: DG.ViewBase, entityId: string) {
  return aiTools.showMetabolicEntity(view, entityId);
}

//name: highlightShortestPath
//description: Find and highlight the k-th shortest path between two metabolites on the map. Direction matters
//tags: metabolicViewFunction
//input: view view
//input: string fromMetabolite
//input: string toMetabolite
//input: int k { optional: true; description: 1 for the shortest path, 2 for the second shortest and so on }
//output: dynamic result
export function highlightShortestPath(view: DG.ViewBase, fromMetabolite: string, toMetabolite: string, k?: number) {
  return aiTools.highlightShortestPath(view, fromMetabolite, toMetabolite, k);
}

//name: listMetabolicAnalyses
//description: List the saved metabolic analyses that loadMetabolicAnalysis can open
//tags: metabolicViewFunction
//input: view view
//output: dynamic result
export async function listMetabolicAnalyses(view: DG.ViewBase) {
  return aiTools.listMetabolicAnalyses(view);
}

//name: saveMetabolicAnalysis
//description: Save the current analysis (map, model, selection and flux data) under a name. Refuses to overwrite an existing analysis unless overwrite is true
//tags: metabolicViewFunction
//input: view view
//input: string name
//input: string description { optional: true }
//input: bool overwrite { optional: true }
//output: dynamic result
export async function saveMetabolicAnalysis(view: DG.ViewBase, name: string, description?: string, overwrite?: boolean) {
  return aiTools.saveMetabolicAnalysis(view, name, description, overwrite);
}

//name: loadMetabolicAnalysis
//description: Load a saved analysis by name, replacing the current map, model and flux data. ALWAYS confirm with the user first when unsaved work could be lost
//tags: metabolicViewFunction
//input: view view
//input: string name
//output: dynamic result
export async function loadMetabolicAnalysis(view: DG.ViewBase, name: string) {
  return aiTools.loadMetabolicAnalysis(view, name);
}

// //name: glpkFBA
// //description: Run FBA using GLPK
// export async function glpkFBA() {
//   console.time('FBA');
//   //const res1 = modelFromJsonData(icho as unknown as CobraModelData).sampleExtremePoints();
//   const res = await WorkerCobraSolver.get_extreme_points(model as unknown as CobraModelData);
//   console.timeEnd('FBA');
//   console.log(res);
// }

// //name: sampleReactionsWasm
// //description: Run FBA using GLPK
// export async function samplerWasm() {
//   console.time('sampler');
//   //const res1 = modelFromJsonData(icho as unknown as CobraModelData).sampleExtremePoints();
//   const cobraModel = model as unknown as CobraModelData;
//   const res = await WorkerCobraSolver.runSampling(cobraModel, 10000);
//   console.timeEnd('sampler');

//   const reactions = cobraModel.reactions;
//   const columns = reactions.map((r, j) => {
//     const col = DG.Column.float(r.id, 10000);
//     col.init((i) => res[i * reactions.length + j]);
//     return col;
//   });
//   const table = DG.DataFrame.fromColumns(columns);
//   table.name = 'Sampler results';
//   grok.shell.addTableView(table);
//   // console.log(res);
// }
