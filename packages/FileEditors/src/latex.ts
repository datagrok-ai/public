import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import '../css/latex/css/katex.css';
import '../css/latex/css/article.css';

// import katex from './katex/katex.min.js';
// import './katex/katex.min.css';

export async function texFilePreview(file: DG.FileInfo): Promise<DG.View> {
  try {
    const latexText = await file.readAsString();
    const generator = new HtmlGenerator({hyphenate: false});
    const doc = parse(latexText, {generator: generator}).htmlDocument();
    const div = ui.div([]);
    const v = DG.View.create();
    v.append(div);
    v.name = file.fileName;
    div.innerHTML = doc.documentElement.outerHTML;

    return v;
  } catch (err) {
    const v = DG.View.create();
    v.append(ui.div([ui.h2('The file is corrupted and cannot be opened!')]));
    v.name = file.fileName;

    if (err instanceof Error)
      grok.shell.error(err.message);

    return v;
  }
  // try {
  //   const latexText = await file.readAsString();
  //   const div = ui.div([]);

  //   katex.render(latexText, div, {
  //     throwOnError: false,
  //   });

  //   const v = DG.View.create();
  //   v.append(div);
  //   v.name = file.fileName;

  //   return v;
  // } catch (err) {
  //   const v = DG.View.create();
  //   v.append(ui.div([ui.h2('The file is corrupted and cannot be opened!')]));
  //   v.name = file.fileName;

  //   if (err instanceof Error)
  //     grok.shell.error(err.message);

  //   return v;
  // }
} // texFilePreview
