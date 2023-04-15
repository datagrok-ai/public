import * as DG from 'datagrok-api/dg';
import {COL_NAMES, GENERATED_COL_NAMES, SEQUENCE_TYPES} from '../model/registration/const';
import * as grok from 'datagrok-api/grok';
import {RegistrationSequenceParser} from '../model/registration/registration-sequence-parser';
import {isValidSequence} from '../model/code-converter/conversion-validation-tools';
import {batchMolWeight, molecularWeight, saltMass, saltMolWeigth} from '../model/registration/calculations';
import {weightsObj} from '../hardcode-to-be-eliminated/map';

export class SdfColumnsExistsError extends Error {
  constructor(message: string) {
    super(message);
  }
}

// todo: rename, refactor
export class RegistrationColumnsHandler {
  constructor(
    private df: DG.DataFrame,
    private onError: (rowI: number, err: any) => void
  ) { }

  addColumns(
    saltNamesList: string[], saltsMolWeightList: number[]
  ): DG.DataFrame {
    const sequenceCol = this.df.getCol(COL_NAMES.SEQUENCE);
    const saltCol = this.df.getCol(COL_NAMES.SALT);
    const equivalentsCol = this.df.getCol(COL_NAMES.EQUIVALENTS);
    const typeCol = this.df.getCol(COL_NAMES.TYPE);
    const chemistryNameCol = this.df.getCol(COL_NAMES.CHEMISTRY_NAME);

    this.validate();
    this.removeEmptyRows(sequenceCol);

    this.df.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => {
      let res: string = '';
      try {
        res = ([SEQUENCE_TYPES.DUPLEX, SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) ?
          chemistryNameCol.get(i) :
          sequenceCol.get(i);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

    this.df.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => {
      let res: string = '';
      const parser = new RegistrationSequenceParser();
      try {
        if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(typeCol.get(i))) {
          res = sequenceCol.get(i);
        } else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
          const obj = parser.getDuplexStrands(sequenceCol.get(i));
          res = `${chemistryNameCol.get(i)}; duplex of SS: ${obj.ss} and AS: ${obj.as}`;
        } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) {
          const obj = parser.getDimerStrands(sequenceCol.get(i));
          res = `${chemistryNameCol.get(i)}; duplex of SS: ${obj.ss} and AS1: ${obj.as1} and AS2: ${obj.as2}`;
        }
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

    this.df.columns.addNewFloat(COL_NAMES.COMPOUND_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      const parser = new RegistrationSequenceParser();
      try {
        if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(typeCol.get(i))) {
          res = (isValidSequence(sequenceCol.get(i), null).indexOfFirstInvalidChar == -1) ?
            molecularWeight(sequenceCol.get(i), weightsObj) :
            DG.FLOAT_NULL;
        } else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
          const obj = parser.getDuplexStrands(sequenceCol.get(i));
          const strands = Object.values(obj);
          res = strands.every((seq) => {
            // console.log('sequence:', seq);
            const invalidChar = isValidSequence(seq, null).indexOfFirstInvalidChar;
            const validity = invalidChar === -1;
            // console.log('validity:', validity);
            // console.log(`invalid char ${invalidChar}:`, seq.substr(invalidChar-5, 10));
            return validity;
          }) ?
            molecularWeight(obj.ss, weightsObj) + molecularWeight(obj.as, weightsObj) :
            DG.FLOAT_NULL;
        } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) {
          const obj = parser.getDimerStrands(sequenceCol.get(i));
          res = (Object.values(obj).every((seq) => isValidSequence(seq, null).indexOfFirstInvalidChar == -1)) ?
            molecularWeight(obj.ss, weightsObj) + molecularWeight(obj.as1, weightsObj) +
            molecularWeight(obj.as2, weightsObj) :
            DG.FLOAT_NULL;
        }
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

    this.df.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = saltMass(saltNamesList, saltsMolWeightList, equivalentsCol, i, saltCol);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

    this.df.columns.addNewFloat(COL_NAMES.SALT_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = saltMolWeigth(saltNamesList, saltCol, saltsMolWeightList, i);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

    const compoundMolWeightCol = this.df.getCol(COL_NAMES.COMPOUND_MOL_WEIGHT);
    const saltMassCol = this.df.getCol(COL_NAMES.SALT_MASS);
    this.df.columns.addNewFloat(COL_NAMES.BATCH_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = batchMolWeight(compoundMolWeightCol, saltMassCol, i);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

    return this.df;
  }

  private validate() {
    if (GENERATED_COL_NAMES.some((colName) => this.df.columns.contains(colName)))
      throw new SdfColumnsExistsError('Columns already exist');
  }

  private removeEmptyRows(colToCheck: DG.Column): void {
    for (let i = 0; i < this.df.rowCount; ++i) {
      if (colToCheck.getString(i) === '')
        this.df.rows.removeAt(i, 1, false);
    }
  }
}
