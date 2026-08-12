/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

//name: Scaffold
//description: A customizable DMTA workflow skeleton with embedded AIDD — built on Datagrok.
//meta.role: app
//input: string path {meta.url: true; optional: true}
//output: view result
export async function scaffoldApp(path?: string): Promise<DG.ViewBase> {
  const view = DG.View.create();
  view.name = 'Scaffold';

  const vision = await _package.files.readAsText('vision.md');
  const content = ui.markdown(vision);
  content.style.maxWidth = '900px';
  content.style.margin = '0 auto';
  content.style.padding = '24px';

  const scroller = ui.div([content]);
  scroller.style.overflow = 'auto';
  scroller.style.height = '100%';

  view.root.appendChild(scroller);
  return view;
}
