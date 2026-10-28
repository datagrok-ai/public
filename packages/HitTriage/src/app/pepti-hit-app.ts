/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {PeptiHitTemplate} from './types';
import {PeptiHitHelmColName} from './consts';
import {HitDesignApp} from './hit-design-app';
import {PeptiHitInfoView} from './pepti-hits-views/info-view';
import {createNewHitDesignDataFrame} from './accordeons/new-hit-design-campaign-accordeon';
import type {CampaignSessionOptions} from './shared/campaign-session';

/** Hit Design for peptides: a design is a HELM sequence, converted to the molecule the compute and the V-iD are of. */
export class PeptiHitApp extends HitDesignApp<PeptiHitTemplate> {
  _helmColName: string = PeptiHitHelmColName;
  constructor(c: DG.FuncCall) {
    super(c, 'PeptiHit', (app) => new PeptiHitInfoView(app as PeptiHitApp));
  }

  protected override newCampaignFrame(template: PeptiHitTemplate): DG.DataFrame {
    return createNewHitDesignDataFrame(template, true);
  }

  /** The Helm column, tagged as HELM for its renderer and editor. */
  protected override prepareTable(df: DG.DataFrame): void {
    this._helmColName = df.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)?.name ?? PeptiHitHelmColName;
    const helmCol = df.col(this._helmColName);
    if (!helmCol)
      throw new Error('No helm column found');
    helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    helmCol.setTag('units', 'helm');
    helmCol.setTag('.alphabetIsMultichar', 'true');
    helmCol.setTag('cell.renderer', 'helm');
  }

  /** Hit Design's, with the HELM editor on the Helm column, the V-iD registered after the conversion, no designer,
   * and the template's layout for any campaign without its own. */
  protected override sessionOptions(): CampaignSessionOptions {
    const layout = this.campaign?.layout ?? this.template?.layoutViewState;
    return {entryColumn: this.helmColName, editor: 'Helm:editMoleculeCell', duplicateLabel: 'Duplicate Sequence',
      registerFirst: false, designer: false, layout};
  }

  get helmColName(): string {
    return this._helmColName ??
        this.dataFrame!.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)?.name;
  }

  /** The sequence converted to the molecule (Bio's toAtomicLevel) in the row, then Hit Design's calculation of it. */
  async performSingleCellCalculations(cellIndex: number, cellValue?: string) {
    try {
      if (!cellValue)
        throw new Error('No cell value provided');
      const df = this.dataFrame!;
      const run = this.nextRun(cellIndex);
      const col = DG.Column.fromStrings(this.helmColName, [cellValue]);
      const table = DG.DataFrame.fromColumns([col]);
      col.semType = DG.SEMTYPE.MACROMOLECULE;
      col.setTag('units', 'helm');
      col.setTag('.alphabetIsMultichar', 'true');
      await grok.functions.call('Bio:toAtomicLevel', {table: table, seqCol: col, nonlinear: true});
      if (!this.isLatestRun(cellIndex, run) || this.dataFrame !== df)
        return;
      const molCol = table.columns.names().find((c) => c.toLowerCase() !== this.helmColName.toLowerCase());
      const mol = molCol ? table.col(molCol)?.get(0) : null;
      if (!mol)
        throw new Error('Failed to convert sequence to atomic level');
      df.set(this.molColName, cellIndex, mol);
      await super.performSingleCellCalculations(cellIndex, mol);
    } catch (e) {
      console.error(e);
    }
  }
}
