import * as DG from 'datagrok-api/dg';
import {COL_NAMES, GENERATED_COL_NAMES, SEQUENCE_TYPES} from './const';
import * as grok from 'datagrok-api/grok';
import {RegistrationSequenceParser} from './sequence-parser';
import {getBatchMolWeight, getMolWeight, getSaltMass, getSaltMolWeigth} from './calculations';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';
import {FormatDetector} from '../parsing-validation/format-detector';
import {SequenceValidator} from '../parsing-validation/sequence-validator';

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
  ) {
    this.sequenceCol = this.df.getCol(COL_NAMES.SEQUENCE);
    this.saltCol = this.df.getCol(COL_NAMES.SALT);
    this.equivalentsCol = this.df.getCol(COL_NAMES.EQUIVALENTS);
    this.typeCol = this.df.getCol(COL_NAMES.TYPE);
    this.chemistryNameCol = this.df.getCol(COL_NAMES.CHEMISTRY_NAME);

    this.validate();
    this.removeEmptyRows(this.sequenceCol);
  }

  private sequenceCol: DG.Column;
  private saltCol: DG.Column;
  private equivalentsCol: DG.Column;
  private typeCol: DG.Column;
  private chemistryNameCol: DG.Column;

  addColumns(
    saltNamesList: string[], saltsMolWeightList: number[]
  ): DG.DataFrame {
    this.addCompoundNameCol();
    this.addCompoundCommentsCol();
    this.addCompoundMolWeightCol();
    this.addSaltMassCol(saltNamesList, saltsMolWeightList);
    this.addSaltMolWeightCol(saltNamesList, saltsMolWeightList);
    this.addBatchMolWeightCol();

    return this.df;
  }

  private addCompoundNameCol(): void {
    this.df.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => {
      let cellValue: string = '';
      try {
        const seqTypeSpecified = [SEQUENCE_TYPES.DUPLEX, SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(this.typeCol.get(i));
        cellValue = (seqTypeSpecified) ?  this.chemistryNameCol.get(i) : this.sequenceCol.get(i);
      } catch (err) {
        this.onError(i, err);
      }
      return cellValue;
    });
  }

  private addCompoundCommentsCol(): void {
    this.df.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => {
      let res: string = '';
      const parser = new RegistrationSequenceParser();
      try {
        const seqType = this.typeCol.get(i);
        if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(seqType)) {
          res = this.sequenceCol.get(i);
        } else if (seqType == SEQUENCE_TYPES.DUPLEX) {
          const obj = parser.getDuplexStrands(this.sequenceCol.get(i));
          res = `${this.chemistryNameCol.get(i)}; duplex of SS: ${obj.ss} and AS: ${obj.as}`;
        } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(seqType)) {
          const obj = parser.getDimerStrands(this.sequenceCol.get(i));
          res = `${this.chemistryNameCol.get(i)}; duplex of SS: ${obj.ss} and AS1: ${obj.as1} and AS2: ${obj.as2}`;
        }
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });
  }

  private addCompoundMolWeightCol(): void {
    this.df.columns.addNewFloat(COL_NAMES.COMPOUND_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      const seqType = this.typeCol.get(i);
      const parser = new RegistrationSequenceParser();
      const codesToWeightsMap = MonomerLibWrapper.getInstance().getCodesToWeightsMap();
      try {
        if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(seqType)) {
          const seq = this.sequenceCol.get(i);
          res = (isValid(seq)) ? getMolWeight(seq, codesToWeightsMap) : DG.FLOAT_NULL;
        } else if (seqType == SEQUENCE_TYPES.DUPLEX) {
          const seq = this.sequenceCol.get(i);
          const obj = parser.getDuplexStrands(seq);
          const strands = Object.values(obj);
          res = strands.every((seq) => isValid(seq)) ?
            getMolWeight(obj.ss, codesToWeightsMap) + getMolWeight(obj.as, codesToWeightsMap) : DG.FLOAT_NULL;
        } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(seqType)) {
          const seq = this.sequenceCol.get(i);
          const obj = parser.getDimerStrands(seq);
          res = (Object.values(obj).every((seq) => isValid(seq))) ?
            getMolWeight(obj.ss, codesToWeightsMap) + getMolWeight(obj.as1, codesToWeightsMap) +
            getMolWeight(obj.as2, codesToWeightsMap) :
            DG.FLOAT_NULL;
        }
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

  }

  private addSaltMassCol(saltNamesList: string[], saltsMolWeightList: number[]): void {
    this.df.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = getSaltMass(saltNamesList, saltsMolWeightList, this.equivalentsCol, i, this.saltCol);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });

  }

  private addSaltMolWeightCol(saltNamesList: string[], saltsMolWeightList: number[]): void {
    this.df.columns.addNewFloat(COL_NAMES.SALT_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = getSaltMolWeigth(saltNamesList, this.saltCol, saltsMolWeightList, i);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });
  }

  private addBatchMolWeightCol(): void {
    const compoundMolWeightCol = this.df.getCol(COL_NAMES.COMPOUND_MOL_WEIGHT);
    const saltMassCol = this.df.getCol(COL_NAMES.SALT_MASS);
    this.df.columns.addNewFloat(COL_NAMES.BATCH_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = getBatchMolWeight(compoundMolWeightCol, saltMassCol, i);
      } catch (err) {
        this.onError(i, err);
      }
      return res;
    });
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

function isValid(sequence: string): boolean {
  const validator = new SequenceValidator(sequence);
  const format = (new FormatDetector(sequence).getFormat());
  return (format === null) ? false : validator.isValidSequence(format!);
}
