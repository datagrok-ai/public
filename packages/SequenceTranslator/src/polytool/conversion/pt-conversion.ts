
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {Rules} from './pt-rules';
import {Chain} from './pt-chain';
import {_package} from '../../package';

/** The main PolyTool convert engine. Returns list of Helms. Covered with tests. */
export function doPolyToolConvert(sequences: string[], rules: Rules, helmHelper: IHelmHelper): [string[], boolean[]] {
  const helms = new Array<string>(sequences.length);
  const isLinear = new Array<boolean>(sequences.length);
  for (let i = 0; i < sequences.length; i++) {
    try {
      if (sequences[i] == null) { helms[i] = ''; } else {
        const chain = Chain.fromSeparator(sequences[i], helmHelper);
        chain.applyRules(rules);
        isLinear[i] = chain.monomersUnderRules.length > 1 || chain.linkagesUnderRules.length > 0 ? false : true;
        helms[i] = chain.getHelm();
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
      helms[i] = '';
    }
  }
  return [helms, isLinear];
}
