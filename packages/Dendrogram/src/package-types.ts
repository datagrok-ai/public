import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {LoggerWrapper} from '@datagrok-libraries/bio/src/utils/logger';

export class DendrogramPackage extends DG.Package {
  constructor(opts: { debug: boolean } = {debug: false}) {
    super();
    // @ts-ignore
    super._logger = new LoggerWrapper(super.logger, opts.debug, 'Dendrogram: ');
  }
}
