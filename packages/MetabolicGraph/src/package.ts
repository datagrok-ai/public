/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import map from './escher/docs/_static/example_data/S5_iJO1366.Glycolysis_PPP_AA_Nucleotides.json';
import model from './escher/docs/_static/example_data/iJO1366.json';

export const _package = new DG.Package();

//name: MetabolicGraph
//tags: app
//description: Metabolic graph application
//input: string path {meta.url: true; optional: true}
//input: string filter {optional: true}
//output: view v
//meta.browsePath: Misc
export function metabolicGraphApp(path?: string, filter?: string): DG.ViewBase {
  const view = DG.View.create('d4-escher-container');
  view.name = 'Metabolic Graph App';
  setTimeout(() => {
    //@ts-ignore
    window.escher.Builder(map, model, null, window.escher.libs.d3_select('.d4-escher-container'),
      {scroll_behavior: 'zoom', fill_screen: false, never_ask_before_quit: true});
  }, 500);
  return view;
}
