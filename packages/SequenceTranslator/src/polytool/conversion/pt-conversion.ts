
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {Rules, RuleReaction} from '../pt-rules';
import {InvalidReactionError, MonomerNotFoundError} from '../types';
import {Chain} from './pt-chain';

import {_package} from '../../package';

/** The main PolyTool convert engine. Returns list of Helms. Covered with tests. */
export function doPolyToolConvert(sequences: string[], rules: Rules, helmHelper: IHelmHelper): string[] {
  const helms = new Array<string>(sequences.length);
  for (let i = 0; i < sequences.length; i++) {
    try {
      if (sequences[i] == null) { helms[i] = ''; } else {
        const chain = Chain.fromNotation(sequences[i], rules, helmHelper);
        helms[i] = chain.getHelm();
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
      helms[i] = '';
    }
  }
  return helms;
}
