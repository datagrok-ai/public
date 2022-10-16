import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class HelmWebEditor {
  editorView: any;
  webEditor: any;
  host: any;
  editor: any;
  w = 200;
  h = 100;
  constructor() {
    this.editorView = ui.div();
    org.helm.webeditor.MolViewer.molscale = 0.8;
    this.webEditor = new scil.helm.App(this.editorView, {
      showabout: false,
      mexfontsize: '90%',
      mexrnapinontab: true,
      topmargin: 20,
      mexmonomerstab: true,
      sequenceviewonly: false,
      mexfavoritefirst: true,
      mexfilter: true
    });
    const sizes = this.webEditor.calculateSizes();
    this.webEditor.canvas.resize(sizes.rightwidth - 100, sizes.topheight - 210);
    let s = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
    //@ts-ignore
    scil.apply(this.webEditor.sequence.style, s);
    //@ts-ignore
    scil.apply(this.webEditor.notation.style, s);
    s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + this.webEditor.toolbarheight) + 'px'};
    //@ts-ignore
    scil.apply(this.webEditor.properties.parent.style, s);
    this.webEditor.structureview.resize(sizes.rightwidth, sizes.bottomheight + this.webEditor.toolbarheight);
    this.webEditor.mex.resize(sizes.topheight - 80);
    this.host = ui.div([], {style: {width: `${this.w}px`, height: `${this.h}px`}});
    this.editor = new JSDraw2.Editor(this.host, {width: this.w, height: this.h, viewonly: true});
  }
}

