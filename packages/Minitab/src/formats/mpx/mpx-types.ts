import {MwxWorksheet} from '../mwx/mwx-types';

export interface MpxProject {
  creator: string;
  comments: string;
  history: string[];
  worksheets: MwxWorksheet[];
}
