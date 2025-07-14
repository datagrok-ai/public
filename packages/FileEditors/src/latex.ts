// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import {_package} from './package';

//import '../css/latex/css/katex.css';
//import '../css/latex/css/article.css';


// import katex from './katex/katex.min.js';
// import './katex/katex.min.css';

export async function addTexToDiv(latexText: string, div: HTMLDivElement): Promise<void> {
  try {
    const generator = new HtmlGenerator({hyphenate: false});
    const htmlGenerator = parse(latexText, {generator: generator});

    const doc = htmlGenerator.htmlDocument(`${_package.webRoot}css/latex/`);

    // const linksArr = doc.head.querySelectorAll('link');
    // linksArr[0].href = `${_package.webRoot}css/latex/css/katex.css`;
    // linksArr[1].href = `${_package.webRoot}css/latex/css/article.css`;

    //console.log(doc.head.querySelectorAll('link'));
    // for (const el of doc.head.querySelectorAll('link'))
    //   el.href = el.href.replace(window.location.origin + window.location.pathname, `${_package.webRoot}css/latex/css/`);
      //console.log(el.href.toString());


    div.style.width = '100%';
    div.style.height = '100%';

    //v.name = file.fileName;

    // const host = div;//document.createElement('div');
    // const shadow = host.attachShadow({mode: 'open'});
    // shadow.innerHTML = doc.documentElement.outerHTML; //`<p>This is a shadow DOM component</p>`;
    // document.body.appendChild(host);

    //div.outerHTML = doc.documentElement.outerHTML;

    //doc.documentElement.outerHTML;
    // const html = `
    // <html>
    //   <head>
    //   ${_package.webRoot}/latex.ts
    //   </head>
    //   <body>
    //   </body>
    // </html>`;
    const iframe = document.createElement('iframe');
    //this.iframe = document.createElement('iframe');
    //const link = doc.createElement('link');
    //link.href = `${_package.webRoot}/latex.ts`;
    //ui.iframe({src: html, width: '100%', height: '100%'});
    div.append(iframe);

    iframe.srcdoc = doc.documentElement.outerHTML;

    //iframe.contentDocument.head.append(link);

    //const iframe = document.getElementById('my-iframe');
    // iframe.onload = function() {
    //   const doc = iframe.contentDocument || iframe.contentWindow.document;

    //   // Create and append a <link> element
    //   const link = doc.createElement('link');
    //   link.rel = 'stylesheet';
    //   link.href = 'path/to/your-styles.css'; // Your CSS path here

    //   doc.head.appendChild(link); // Insert into <head>
    // };

    //return v;
  } catch (err) {
    // const v = DG.View.create();
    // v.append(ui.div([ui.h2('The file is corrupted and cannot be opened!')]));
    // v.name = file.fileName;

    // if (err instanceof Error)
    //   grok.shell.error(err.message);

    // return v;
  }
} // texFilePreview
