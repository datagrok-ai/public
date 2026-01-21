/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {parsePath, loadStateProxy, sampleReactions, saveStateDialog, loadAnalisisDialog, handleReactionDataUpload, runFBADialog} from './utils';
import map from './maps/E_coli_Core_metabolism_map.json';
import model from './maps/E_coli_core_cobra.json';
import type {MapData, CobraModelData, SettingsType} from '../escher_src/src/ts/types';
import type {BuilderType, BuilderConstructor} from '../escher_src/src/Builder';
import { modelFromJsonData } from './FBA/cobraSolver';
import { WorkerCobraSolver } from './cobra';
import { sampleReactionsWasm } from './cobra/sampler-wrapper';

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
  const view = DG.View.create('d4-escher-container');
  view.name = 'Metabolic Graph';
  setTimeout(() => {
    //@ts-ignore
    const mapData = map as MapData;
    const modelData = model as unknown as CobraModelData;
    const Builder = window.escher.Builder as BuilderConstructor;
    const b = new Builder(mapData, modelData, null, window.escher.libs.d3_select('.d4-escher-container'),
      {scroll_behavior: 'zoom', fill_screen: false,
        never_ask_before_quit: true,
        samplingFunction: (mp: CobraModelData) => sampleReactions(mp, b),
        saveAction: () => saveStateDialog(b, undefined),
        loadAction: () => loadAnalisisDialog(b),
        runFBA: async () => {runFBADialog(b)},
        pathFindingDisabled: true,
      });
    handleReactionDataUpload(view, b);
    setTimeout(async () => {
      if (!b.map)
        return;
      const pathParsed = parsePath(path);
      if (pathParsed) {
        if (!_package.files.exists(`campaigns/${pathParsed}.json`)) {
          grok.shell.error(`File ${pathParsed}.json not found`);
          return;
        }
        const data = await _package.files.readAsText(`campaigns/${pathParsed}.json`);
        loadStateProxy(b, data, pathParsed);
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
  const view = DG.View.create('d4-escher-container');
  view.name = 'Metabolic Graph';
  const mapJSON = await file.readAsString();
  const mapData = JSON.parse(mapJSON) as MapData;
  const modelData = null;
  setTimeout(() => {
    const Builder = window.escher.Builder as BuilderConstructor;
    const b = new Builder(mapData, modelData, null, window.escher.libs.d3_select('.d4-escher-container'),
      {scroll_behavior: 'zoom', fill_screen: false,
        never_ask_before_quit: true,
        samplingFunction: (mp: CobraModelData) => sampleReactions(mp, b),
        saveAction: () => saveStateDialog(b, undefined),
        loadAction: () => loadAnalisisDialog(b),
        runFBA: async () => {runFBADialog(b)},
        pathFindingDisabled: true,
      });
  }, 500);
  return view;
}

//name: escherFileViewerCheck
//input: string content
//output: bool result
export async function escherFileViewerCheck(content: string) {
  return (content?.length ?? 1e12) < 1e7 && !!content?.startsWith('[') && !!content?.includes('"https://escher.github.io/escher/jsonschema/1-0-0#"');
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