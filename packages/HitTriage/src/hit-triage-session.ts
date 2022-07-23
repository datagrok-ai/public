/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class HitTriageSession {
  sourceQuery?: DG.FuncCall;
  sourceDataFrame?: DG.DataFrame;

  static demo(): HitTriageSession {
    const session = new HitTriageSession();
    session.sourceDataFrame = grok.data.demo.molecules(20000);
    return session;
  }
}