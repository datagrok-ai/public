import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

import {Chain} from './conversion/pt-chain';
import {getPolyToolUnruleDialog} from './pt-unrule-dialog';
import {Rules} from './pt-rules';

import {_package} from '../package';
import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

export async function polyToolUnruleUI(): Promise<void> {
  let dialog: DG.Dialog;
  try {
    dialog = await getPolyToolUnruleDialog();
    dialog.show();
  } catch (err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.warning('To run PolyTool Unrule, open a dataframe with Helm');
    _package.logger.error(errMsg, undefined, errStack);
  }
}

/** Returns list of harmonized sequences. Covered with tests. */
export function doPolyToolUnrule(helms: string[], rules: Rules, helmHelper: IHelmHelper): string[] {
  const resHrzSeqList = new Array<string>(helms.length);
  for (let i = 0; i < helms.length; ++i) {
    if (!helms[i])
      resHrzSeqList[i] = '';
    else {
      const chain = Chain.parseHelm(helms[i], helmHelper);
      resHrzSeqList[i] = chain.getNotation();
    }
  }
  return resHrzSeqList;
}
