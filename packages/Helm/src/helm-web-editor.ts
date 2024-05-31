import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HelmType} from '@datagrok/js-draw-lite/src/types/org';
import {IEditor} from '@datagrok/js-draw-lite/src/types/jsdraw2';

import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';

import {JSDraw2Module} from './types';

declare const JSDraw2: JSDraw2Module;

export class HelmWebEditor implements IHelmWebEditor {
  editor: IEditor<HelmType>;
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

