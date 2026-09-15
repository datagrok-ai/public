import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Base} from './base';
import {callThings} from './utils.js';
import {shared} from './common';
import {something, other as renamed} from '@datagrok-libraries/utils/src/test';
import '@datagrok-libraries/utils/src/nowhere';
import {helper} from '@datagrok/plain/src/package';
import {Observable} from 'rxjs';
import './styles.css';
import {gone} from './missing';

// DG.Column in a comment does not count
export class DemoViewer extends DG.JsViewer implements DG.IDisposable {
  dispose(): void {
  }

  async render(): Promise<void> {
    const df = DG.DataFrame.fromCsv('a');
    ui.div([ui.input.string('x')]);
    grok.shell.info(df.name);
    grok.shell.info(DG.SEMTYPE.MOLECULE);
    grok.log.info(String(DG.LogLevel.Info));
    grok.chem.similarity('C', 'CC');
    const first = DG.Nowhere.x;
    const second = grok.nowhere.y;
    DG.U2.Control.forElement(this); DG.Utils.Thing.go();
    const common = await import('./common');
    console.log(first, second, common, shared, something, renamed, helper, callThings, Observable, gone);
  }
}

export class Child extends Base {
}

export class Lost extends Missing {
}
