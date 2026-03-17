/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {knowledgeFiles} from './knowledge';
import {renderMarkdown, RenderOptions} from './md-renderer';
import {createPipelineDiagram} from './pipeline-diagram';
import {glossary} from './glossary';

export const _package = new DG.Package();

function renderDiagram(src: string): HTMLElement | null {
  if (src === 'popPKPD_pipeline_personas_en.svg')
    return createPipelineDiagram(glossary);
  return null;
}

function openMdView(fileName: string): void {
  const text = knowledgeFiles[fileName];
  if (!text) {
    grok.shell.error(`File not found: ${fileName}`);
    return;
  }
  const title = fileName.replace(/_/g, ' ').replace('.md', '');
  const view = grok.shell.newView(title);
  const options: RenderOptions = {onMdLink: openMdView, renderDiagram, glossary};
  const content = renderMarkdown(text, options);
  view.append(ui.panel([content]));
}

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: PopPK
//description: Interactive presentation of population pharmacokinetics modeling
//tags: app
export function popPk() {
  openMdView('PopPK_Studio_Overview.md');
}
