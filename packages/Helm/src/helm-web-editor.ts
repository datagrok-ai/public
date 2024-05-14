import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as JSDraw2 from 'JSDraw2';

import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';

export class HelmWebEditor implements IHelmWebEditor {
  editor: JSDraw2.Editor;
  host: HTMLDivElement;

  w = 200;
  h = 100;

  constructor() {
    this.host = ui.div([], {style: {width: `${this.w}px`, height: `${this.h}px`}});
    this.editor = new JSDraw2.Editor(this.host, {width: this.w, height: this.h, viewonly: true});
  }

  resizeEditor(w: number, h: number) {
    this.editor.setSize(w, h);
  }
}

