/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {parsePath, loadStateProxy, sampleReactions, saveStateDialog, loadAnalisisDialog} from './utils';
import map from './maps/S5_iJO1366.Glycolysis_PPP_AA_Nucleotides.json';
import model from './maps/iJO1366.json';
import type {MapData, CobraModelData, SettingsType} from '../escher_src/src/ts/types';
import type {BuilderType, BuilderConstructor} from '../escher_src/src/Builder';

export const _package = new DG.Package();

declare global {
  // eslint-disable-next-line no-unused-vars
  interface Window {
    escher: any;
    builder: any
  }
}

//name: MetabolicGraph
//tags: app
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
      {scroll_behavior: 'zoom', fill_screen: false, saveAction: () => saveStateDialog(b, undefined),
        loadAction: () => loadAnalisisDialog(b), pathFindingDisabled: true,
        never_ask_before_quit: true});
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
    }, 200);
  }, 500);
  return view;
}

//name: EscherFileViewer
//tags: fileViewer
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
      {scroll_behavior: 'zoom', fill_screen: false, saveAction: () => saveStateDialog(b, undefined),
         pathFindingDisabled: true,
        loadAction: () => loadAnalisisDialog(b),
        never_ask_before_quit: true});
  }, 500);
  return view;
}

//name: escherFileViewerCheck
//input: string content
//output: bool result
export async function escherFileViewerCheck(content: string) {
  return (content?.length ?? 1e12) < 1e7 && !!content?.startsWith('[') && !!content?.includes('"https://escher.github.io/escher/jsonschema/1-0-0#"');
}
