import * as ui from 'datagrok-api/ui';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import {_package} from './package';

export function getElementWithLatexContent(latexText: string): HTMLIFrameElement {
  const generator = new HtmlGenerator({hyphenate: false});
  const htmlGenerator = parse(latexText, {generator: generator});
  const doc = htmlGenerator.htmlDocument(`${_package.webRoot}css/latex/`) as Document;

  const iframe = ui.iframe({width: '100%', height: '100%'});
  const rightMargin = doc.body.children[1];

  if ((rightMargin !== null) && (rightMargin !== null))
    (rightMargin as HTMLElement).style.display = 'none';

  const contentRegion = doc.body.children[0];

  if ((contentRegion !== null) && (contentRegion !== null)) {
    (contentRegion as HTMLElement).style.transform = 'scale(1.5)';
    (contentRegion as HTMLElement).style.transformOrigin = 'top';
  }

  iframe.srcdoc = doc.documentElement.outerHTML;
  iframe.style.border = 'transparent';

  return iframe;
}
