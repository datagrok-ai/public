/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class HitTriageSessionSummary {

}


export class HitTriageSession {
  project: string = 'New project';

  sourceQuery?: DG.FuncCall;
  sourceDataFrame?: DG.DataFrame;
  sourceMoleculeColumn: string = '';
  sourceType: string = 'file';
  sourceDescription: string = '';

  enrichedDataFrame?: DG.DataFrame;
  enrichmentSteps: string[] = [];
  enrichmentDescriptions: string[] = [];

  filterDescriptions: string[] = [];

  getSummary() {
    return {
      'From': this.sourceType,
      'Path': this.sourceDescription,
      'Total Molecules': this.sourceDataFrame!.rowCount,
      'Filter': this.filterDescriptions,
      'Result Molecules': this.sourceDataFrame!.filter.trueCount
    }
  }

  static demo(): HitTriageSession {
    const session = new HitTriageSession();
    session.project = 'Demo project';
    session.sourceDataFrame = grok.data.demo.molecules(20000);
    session.sourceMoleculeColumn = 'smiles';
    session.sourceDataFrame.meta.detectSemanticTypes().then((_) => {});
    session.sourceType = 'file';
    session.sourceDescription = 'AppData:/HitTriage/campaigns/bfg9000.csv';
    return session;
  }
}