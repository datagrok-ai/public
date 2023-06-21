import * as DG from 'datagrok-api/dg';
import {COL_NAMES, GENERATED_COL_NAMES, SEQUENCE_TYPES} from '../autostart/constants';
import * as grok from 'datagrok-api/grok';
import {removeEmptyRows} from '../utils/helpers';
import {parseStrandsFromDuplexCell, parseStrandsFromTriplexOrDimerCell} from './parse';
import {isValidSequence} from '../sdf-tab/sequence-codes-tools';
import {batchMolWeight, molecularWeight, saltMass, saltMolWeigth} from '../autostart/calculations';
import {weightsObj} from '../hardcode-to-be-eliminated/map';

export class SdfColumnsExistsError extends Error {
  constructor(message: string) {
    super();
  }
}

export function sdfAddColumns(
  df: DG.DataFrame, saltNamesList: string[], saltsMolWeightList: number[], onError: (rowI: number, err: any) => void
): DG.DataFrame {
  const sequenceCol = df.getCol(COL_NAMES.SEQUENCE);
  const saltCol = df.getCol(COL_NAMES.SALT);
  const equivalentsCol = df.getCol(COL_NAMES.EQUIVALENTS);
  const typeCol = df.getCol(COL_NAMES.TYPE);
  const chemistryNameCol = df.getCol(COL_NAMES.CHEMISTRY_NAME);

  if (GENERATED_COL_NAMES.some((colName) => df.columns.contains(colName)))
    throw new SdfColumnsExistsError('Columns already exist');

  df = removeEmptyRows(df, sequenceCol);

  df.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => {
    let res: string = '';
    try {
      res = ([SEQUENCE_TYPES.DUPLEX, SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) ?
        chemistryNameCol.get(i) :
        sequenceCol.get(i);
    } catch (err) {
      onError(i, err);
    }
    return res;
  });

  df.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => {
    let res: string = '';
    try {
      if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(typeCol.get(i))) {
        res = sequenceCol.get(i);
      } else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
        const obj = parseStrandsFromDuplexCell(sequenceCol.get(i));
        res = `${chemistryNameCol.get(i)}; duplex of SS: ${obj.SS} and AS: ${obj.AS}`;
      } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) {
        const obj = parseStrandsFromTriplexOrDimerCell(sequenceCol.get(i));
        res = `${chemistryNameCol.get(i)}; duplex of SS: ${obj.SS} and AS1: ${obj.AS1} and AS2: ${obj.AS2}`;
      }
    } catch (err) {
      onError(i, err);
    }
    return res;
  });

  df.columns.addNewFloat(COL_NAMES.COMPOUND_MOL_WEIGHT).init((i: number) => {
    let res: number = Number.NaN;
    try {
      if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(typeCol.get(i))) {
        res = (isValidSequence(sequenceCol.get(i), null).indexOfFirstNotValidChar == -1) ?
          molecularWeight(sequenceCol.get(i), weightsObj) :
          DG.FLOAT_NULL;
      } else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
        const obj = parseStrandsFromDuplexCell(sequenceCol.get(i));
        res = (Object.values(obj).every((seq) => isValidSequence(seq, null).indexOfFirstNotValidChar == -1)) ?
          molecularWeight(obj.SS, weightsObj) + molecularWeight(obj.AS, weightsObj) :
          DG.FLOAT_NULL;
      } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) {
        const obj = parseStrandsFromTriplexOrDimerCell(sequenceCol.get(i));
        res = (Object.values(obj).every((seq) => isValidSequence(seq, null).indexOfFirstNotValidChar == -1)) ?
          molecularWeight(obj.SS, weightsObj) + molecularWeight(obj.AS1, weightsObj) +
          molecularWeight(obj.AS2, weightsObj) :
          DG.FLOAT_NULL;
      }
    } catch (err) {
      onError(i, err);
    }
    return res;
  });

  df.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) => {
    let res: number = Number.NaN;
    try {
      res = saltMass(saltNamesList, saltsMolWeightList, equivalentsCol, i, saltCol);
    } catch (err) {
      onError(i, err);
    }
    return res;
  });

  df.columns.addNewFloat(COL_NAMES.SALT_MOL_WEIGHT).init((i: number) => {
    let res: number = Number.NaN;
    try {
      res = saltMolWeigth(saltNamesList, saltCol, saltsMolWeightList, i);
    } catch (err) {
      onError(i, err);
    }
    return res;
  });

  const compoundMolWeightCol = df.getCol(COL_NAMES.COMPOUND_MOL_WEIGHT);
  const saltMassCol = df.getCol(COL_NAMES.SALT_MASS);
  df.columns.addNewFloat(COL_NAMES.BATCH_MOL_WEIGHT).init((i: number) => {
    let res: number = Number.NaN;
    try {
      res = batchMolWeight(compoundMolWeightCol, saltMassCol, i);
    } catch (err) {
      onError(i, err);
    }
    return res;
  });

  return df;
}
