import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface INewickHelper {
  newickToDf(newick: string, name: string): DG.DataFrame;
}