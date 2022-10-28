import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class HelmWebEditor {
  host: any;
  editor: any;
  w = 200;
  h = 100;
  constructor() {
    this.host = ui.div([], {style: {width: `${this.w}px`, height: `${this.h}px`}});
    this.editor = new JSDraw2.Editor(this.host, {width: this.w, height: this.h, viewonly: true});
  }

  resizeEditor(w: number, h: number) {
    this.editor.setSize(w, h);
  }

  createWebEditor(value: string) {
    const editorView = ui.div();
    org.helm.webeditor.MolViewer.molscale = 0.8;
    const webEditor = new scil.helm.App(editorView, {
      showabout: false,
      mexfontsize: '90%',
      mexrnapinontab: true,
      topmargin: 20,
      mexmonomerstab: true,
      sequenceviewonly: false,
      mexfavoritefirst: true,
      mexfilter: true
    });
    const sizes = webEditor.calculateSizes();
    webEditor.canvas.resize(sizes.rightwidth - 100, sizes.topheight - 210);
    let s = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
    //@ts-ignore
    scil.apply(webEditor.sequence.style, s);
    //@ts-ignore
    scil.apply(webEditor.notation.style, s);
    s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + webEditor.toolbarheight) + 'px'};
    //@ts-ignore
    scil.apply(webEditor.properties.parent.style, s);
    webEditor.structureview.resize(sizes.rightwidth, sizes.bottomheight + webEditor.toolbarheight);
    webEditor.mex.resize(sizes.topheight - 80);
    setTimeout(function() {
      webEditor.canvas.helm.setSequence(value, 'HELM');
    }, 200);
    return {editorDiv: editorView, webEditor: webEditor};
  }
}

