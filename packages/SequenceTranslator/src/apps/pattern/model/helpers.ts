import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap} from '../../common/data-loading-utils/json-loader';

export function isOverhang(modification: string): boolean {
  const overhangSuffix = '(o)';
  return modification.endsWith(overhangSuffix);
}
