import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {HelmService} from './utils/helm-service';

type HelmWindowType = Window & {
  $helmService?: HelmServiceBase,
}
declare const window: HelmWindowType;

export function _getHelmService(): HelmServiceBase {
  let res = window.$helmService;
  if (!res) res = window.$helmService = new HelmService();
  return res;
}
