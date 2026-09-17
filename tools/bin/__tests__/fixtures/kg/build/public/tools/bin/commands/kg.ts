import * as DG from 'datagrok-api/dg';
import {Ids, fileId} from '../utils/ids';

export function kg(argv: string[]): string {
  return new Ids().parse(fileId(argv[0] ?? DG.SEMTYPE));
}
