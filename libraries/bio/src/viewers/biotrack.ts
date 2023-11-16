import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {IViewer} from './viewer';

export const BiotrackPropsDefault = new class {

}();

export type BiotrackProps = typeof BiotrackPropsDefault;

export interface IBiotrackViewer extends IViewer {

}

declare module 'datagrok-api/dg' {
  interface DataFramePlotHelper {
    // eslint-disable-next-line max-len
    fromType(viewerType: 'Biotrack', options: Partial<BiotrackProps>): Promise<DG.Viewer<BiotrackProps> & IBiotrackViewer>;
  }
}
