/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from './hit-triage-app';
import {PepTriageTemplate, HitTriageTemplate} from './types';
import {PepTriageInfoView} from './pep-triage-views/info-view';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {_package} from '../package';

export class PepTriageApp extends HitTriageApp {
  private _seqColName?: string;

  constructor(c: DG.FuncCall) {
    super(c, 'PepTriage', (app) => new PepTriageInfoView(app as PepTriageApp));
  }

  get seqColName(): string | undefined {return this._seqColName;}

  protected override async prepareDataFrame(template: HitTriageTemplate): Promise<void> {
    const ptTemplate = template as PepTriageTemplate;
    if (!this.dataFrame)
      return;
    await this.dataFrame.meta.detectSemanticTypes();

    // 1. Find sequence column by configured name or by semtype
    const seqColName = ptTemplate.sequenceColumnName;
    let seqCol = seqColName ? this.dataFrame.col(seqColName) : null;
    if (!seqCol)
      seqCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    if (!seqCol) {
      throw new Error(`Sequence column '${seqColName ?? 'Macromolecule'}' not found in the dataset. ` +
        'Make sure the data contains a column with sequences.');
    }
    this._seqColName = seqCol.name;

    // 2. Initialize SeqHandler for proper notation detection
    const pg = DG.TaskBarProgressIndicator.create('Initializing sequences...');
    try {
      const seqHelper = await getSeqHelper();
      const tempSH = seqHelper.getSeqHandler(seqCol);
      await tempSH.refinerPromise;
      // re-init after refinement
      seqHelper.getSeqHandler(seqCol);

      // 3. Check if molecule column already exists
      const configuredMolColName = ptTemplate.moleculeColumnName;
      let molCol = configuredMolColName ? this.dataFrame.col(configuredMolColName) : null;
      if (!molCol)
        molCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);

      if (!molCol) {
        // 4. Convert sequences to molecules via Bio:toAtomicLevel
        pg.update(50, 'Converting sequences to molecules...');
        await grok.functions.call('Bio:toAtomicLevel', {table: this.dataFrame, seqCol: seqCol, nonlinear: true});
        await this.dataFrame.meta.detectSemanticTypes();
        molCol = this.dataFrame.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE).find((a) => a.name.toLowerCase().includes(seqCol.name.toLowerCase())) ?? null;
        if (!molCol) {
          grok.shell.warning('Failed to convert sequences to molecules. Calculations may not work.');
          return;
        }
      }

      // 5. Reorder columns: sequence first, molecule second, rest follows
      this.reorderColumns(seqCol.name, molCol.name);
    } catch (e) {
      _package.logger.error(e);
      grok.shell.error('Error preparing sequence data: ' + (e instanceof Error ? e.message : e));
    } finally {
      pg.close();
    }
  }

  private reorderColumns(seqColName: string, molColName: string) {
    if (!this.dataFrame)
      return;
    const seqCol = this.dataFrame.col(seqColName);
    const molCol = this.dataFrame.col(molColName);
    if (seqCol) {
      this.dataFrame.columns.remove(seqColName);
      this.dataFrame.columns.insert(seqCol, 0);
    }
    if (molCol) {
      this.dataFrame.columns.remove(molColName);
      this.dataFrame.columns.insert(molCol, 1);
    }
  }
}
