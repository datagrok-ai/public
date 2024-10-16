import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getUnusedColName} from '@datagrok-libraries/bio/src/monomer-works/utils';

import {defaultErrorHandler} from '../utils/err-info';
import {doPolyToolUnrule} from './pt-unrule';
import {getRules, RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './pt-rules/pt-rules';
import {PT_ERROR_DATAFRAME, PT_UI_DIALOG_UNRULE, PT_UI_RULES_USED} from './const';
import {_package} from '../package';

type PolyToolUnruleSerialized = {
  rules: string[];
};

export async function getPolyToolUnruleDialog(srcCol?: DG.Column<string>): Promise<DG.Dialog> {
  const subs: Unsubscribable[] = [];
  const destroy = () => {
    for (const sub of subs) sub.unsubscribe();
  };
  try {
    let srcColVal: DG.Column<string> | undefined = srcCol;
    if (!srcColVal) {
      const srcColList = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
      if (!srcColList)
        throw new Error(PT_ERROR_DATAFRAME);
      srcColVal = srcColList[0];
    }
    const srcColInput = ui.input.column('Column', {
      table: srcColVal.dataFrame, value: srcColVal,
      filter: (col: DG.Column) => {
        if (col.semType !== DG.SEMTYPE.MACROMOLECULE) return false;
        const sh = _package.seqHelper.getSeqHandler(col);
        return sh.notation === NOTATION.HELM;
      }
    });
    let ruleFileList: string[];
    const ruleInputs = new RuleInputs(RULES_PATH, RULES_STORAGE_NAME, '.json', {
      onValueChanged: (value: string[]) => { ruleFileList = value;}
    });
    const rulesHeader = ui.inlineText([PT_UI_RULES_USED]);
    const rulesForm = await ruleInputs.getForm();

    const div = ui.divV([
      srcColInput,
      rulesHeader,
      rulesForm
    ]);

    const exec = async (): Promise<void> => {
      try {
        const ruleFileList = await ruleInputs.getActive();
        await polyToolUnrule(srcColInput.value!, ruleFileList);
      } catch (err: any) {
        defaultErrorHandler(err);
      }
    };

    const dialog = ui.dialog(PT_UI_DIALOG_UNRULE)
      .add(div)
      .onOK(() => { exec(); });
    subs.push(dialog.onClose.subscribe(() => {
      destroy();
    }));
    dialog.history(
      /* getInput */ (): PolyToolUnruleSerialized => {
        return {
          rules: ruleFileList,
        };
      },
      /* applyInput */ (x: PolyToolUnruleSerialized): void => {
        ruleInputs.setActive(ruleFileList);
      });
    return dialog;
  } catch (err: any) {
    destroy(); // on failing to build a dialog
    throw err;
  }
}

export async function polyToolUnrule(
  srcCol: DG.Column<string>, ruleFiles: string[]
): Promise<DG.Column> {
  const pi = DG.TaskBarProgressIndicator.create('PolyTool unrule...');
  try {
    const table = srcCol.dataFrame;
    const rules = await getRules(ruleFiles);
    const resHelmList = doPolyToolUnrule(srcCol.toList(), rules);
    const resHelmColName = `harmonized(srcCol.name)`;
    const resHelmCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, resHelmColName, resHelmList,);
    resHelmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    resHelmCol.meta.units = NOTATION.CUSTOM;
    if (table) {
      resHelmCol.name = getUnusedColName(table, resHelmColName);
      table.columns.add(resHelmCol, true);
      await grok.data.detectSemanticTypes(table);
    }
    return resHelmCol;
  } finally {
    pi.close();
  }
}
