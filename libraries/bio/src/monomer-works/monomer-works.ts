import {HelmType, PolymerType} from '../helm/types';
import {HelmTypes, PolymerTypes} from '../helm/consts';

export function helmTypeToPolymerType(helmType: HelmType): PolymerType {
  let polymerType: PolymerType | undefined = undefined;
  switch (helmType) {
  case HelmTypes.BASE:
  case HelmTypes.SUGAR: // r - ribose, d - deoxyribose
  case HelmTypes.LINKER: // p - phosphate
  case HelmTypes.NUCLEOTIDE:
    // @ts-ignore
  case 'nucleotide':
    polymerType = PolymerTypes.RNA;
    break;
  case HelmTypes.AA:
    polymerType = PolymerTypes.PEPTIDE;
    break;
  case HelmTypes.CHEM:
    polymerType = PolymerTypes.CHEM;
    break;
  case HelmTypes.BLOB:
    polymerType = PolymerTypes.BLOB;
    break;
  default:
    polymerType = PolymerTypes.PEPTIDE;
    console.warn(`Unexpected HelmType '${helmType}'`);
  }
  return polymerType;
}
