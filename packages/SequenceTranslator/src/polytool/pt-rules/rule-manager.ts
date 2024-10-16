import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {dfFromRules, getRules, Rules} from './pt-rules';
import {RulesForm} from './pt-rules-form';


export class RulesManager {
  // private adjustColWidths() {
  //   setTimeout(() => {
  //     if (this.tv?.grid) {
  //       this.tv!.grid.col(MONOMER_DF_COLUMN_NAMES.NAME)!.width = 100;
  //       this.tv!.grid.col(MONOMER_DF_COLUMN_NAMES.SYMBOL)!.width = 70;
  //     }
  //   }, 200);
  // }

  rules: Rules;
  ruleDataFrame: DG.DataFrame;
  rulesForm: RulesForm;

  private static instance: RulesManager;

  protected constructor(rules: Rules) {
    this.rules = rules;
    this.ruleDataFrame = dfFromRules(this.rules);

    this.rulesForm = new RulesForm();

    // this._newMonomerForm = new MonomerForm(monomerLibManamger, () => this.activeMonomerLib, async (scrollToRowSymbol?: string) => {
    //   const df = await this.getMonomersDf(this.libInput.value!);
    //   if (this.tv?.dataFrame) {
    //     this.tv.dataFrame = df;
    //     this.adjustColWidths();
    //     if (scrollToRowSymbol != undefined) {
    //       setTimeout(() => {
    //         const col = df.col(MONOMER_DF_COLUMN_NAMES.SYMBOL)!;
    //         const scrollToRow = col.toList().indexOf(scrollToRowSymbol);
    //         if (scrollToRow === -1) return;
    //         this.tv?.grid.scrollToCell(df.columns.byIndex(0), scrollToRow);
    //         df.currentRow = df.rows.get(scrollToRow);
    //       }, 500);
    //     }
    //   }
    //}, () => this.tv?.dataFrame);
  }

  getViewRoot(): DG.ViewBase {
    const tv = DG.TableView.create(this.ruleDataFrame, false);
    const form = this.rulesForm.getForm();
    tv.dockManager.dock(ui.divV([]), DG.DOCK_TYPE.LEFT, null, '', 0.4);
    return tv;
  }

  public static async getInstance(name: string): Promise<RulesManager> {
    if (!this.instance) {
      const rules = await getRules([name]);
      this.instance = new RulesManager(rules);
    }
    return this.instance;
  }
}

