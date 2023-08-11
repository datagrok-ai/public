/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {GapSymbols, UnitsHandler} from './units-handler';
import {ISeqSplitted, SplitterFunc} from './macromolecule/types';
import {NOTATION} from './macromolecule/consts';
import {getSplitterForColumn, splitterAsHelm} from './macromolecule/utils';

export type ConvertFunc = (src: string) => string;

/** Class for handling conversion of notation systems in Macromolecule columns */
export class NotationConverter extends UnitsHandler {
  private _splitter: SplitterFunc | null = null;

  protected get splitter(): SplitterFunc {
    if (this._splitter === null)
      this._splitter = getSplitterForColumn(this.column);
    return this._splitter;
  }

  public toFasta(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.FASTA; }

  public toSeparator(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.SEPARATOR; }

  public toHelm(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.HELM; }

  /**
   *  Convert HELM string to FASTA/SEPARATOR
   *
   * @param {string} helmPolymer    A string to be converted
   * @param {string} tgtNotation    Target notation: FASTA or SEPARATOR
   * @param {string} tgtSeparator   Optional target separator (for HELM ->
   * @param {string | null} tgtGapSymbol   Optional target gap symbol
   * SEPARATOR)
   * @return {string} Converted string
   */
  public convertHelmToFastaSeparator(
    helmPolymer: string, tgtNotation: string, tgtSeparator?: string, tgtGapSymbol?: string
  ): string {
    if (!tgtGapSymbol) {
      tgtGapSymbol = (this.toFasta(tgtNotation as NOTATION)) ?
        GapSymbols[NOTATION.FASTA] :
        GapSymbols[NOTATION.SEPARATOR];
    }

    if (!tgtSeparator)
      tgtSeparator = (this.toFasta(tgtNotation as NOTATION)) ? '' : this.separator;

    const helmWrappersRe = /(R\(|D\(|\)|P)/g;
    const isNucleotide = helmPolymer.startsWith('DNA') || helmPolymer.startsWith('RNA');
    // items can be monomers or helms
    const helmItemsArray = this.splitter(helmPolymer);
    const tgtMonomersArray: string[] = [];
    for (let i = 0; i < helmItemsArray.length; i++) {
      let item = helmItemsArray[i];
      if (isNucleotide)
        item = item.replace(helmWrappersRe, '');
      if (item === GapSymbols[NOTATION.HELM])
        tgtMonomersArray.push(tgtGapSymbol!);
      else if (this.toFasta(tgtNotation as NOTATION) && item.length > 1) {
        // the case of a multi-character monomer converted to FASTA
        const monomer = '[' + item + ']';
        tgtMonomersArray.push(monomer);
      } else
        tgtMonomersArray.push(item);
    }
    return tgtMonomersArray.join(tgtSeparator);
  }

  /** Dispatcher method for notation conversion
   *
   * @param {NOTATION} tgtNotation   Notation we want to convert to
   * @param {string | null} tgtSeparator   Possible separator
   * @return {DG.Column}                Converted column
   */
  public convert(tgtNotation: NOTATION, tgtSeparator?: string): DG.Column {
    const convert: ConvertFunc = this.getConverter(tgtNotation, tgtSeparator);
    const newColumn = this.getNewColumn(tgtNotation, tgtSeparator);
    // assign the values to the newly created empty column
    newColumn.init((rowI: number) => { return convert(this.column.get(rowI)); });
    // newColumn.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
    return newColumn;
  }

  public constructor(col: DG.Column) {
    super(col);
  }

  public getConverter(tgtUnits: NOTATION, tgtSeparator: string | null = null): ConvertFunc {
    if (tgtUnits === NOTATION.SEPARATOR && !tgtSeparator)
      throw new Error(`Target separator is not specified for target units '${NOTATION.SEPARATOR}'.`);

    const srcUh = this;
    if (tgtUnits === NOTATION.FASTA)
      return function(src: string) { return convertToFasta(srcUh, src); };
    if (tgtUnits === NOTATION.HELM)
      return function(src: string) { return convertToHelm(srcUh, src); };
    else if (tgtUnits === NOTATION.SEPARATOR)
      return function(src: string) { return convertToSeparator(srcUh, src, tgtSeparator!); };
    else
      throw new Error();
  }
}

const helmWrappersRe = /[RD]\((\w)\)P?/g;


function convertToFasta(srcUh: UnitsHandler, src: string): string {
  const srcMList: ISeqSplitted = srcUh.isHelm() ? splitterAsHelmNucl(srcUh, src) : srcUh.getSplitter()(src);
  const tgtMList: string[] = new Array<string>(srcMList.length);
  for (const [srcM, mI] of wu.enumerate(srcMList)) {
    let m = srcM;
    if (srcUh.isHelm())
      m = srcM.replace(helmWrappersRe, '$1');

    if (srcUh.isGap(m))
      m = GapSymbols[NOTATION.FASTA];
    else if (m.length > 1)
      m = '[' + srcMList[mI] + ']';

    tgtMList[mI] = m;
  }
  return tgtMList.join('');
}

function convertToSeparator(srcUh: UnitsHandler, src: string, tgtSeparator: string): string {
  const srcMList: ISeqSplitted = srcUh.isHelm() ? splitterAsHelmNucl(srcUh, src) : srcUh.getSplitter()(src);
  const tgtMList: string[] = new Array<string>(srcMList.length);
  const isPT = src.startsWith('PEPTIDE');
  for (const [srcM, mI] of wu.enumerate(srcMList)) {
    let m: string | null = srcM;
    if (srcUh.isGap(m))
      m = GapSymbols[NOTATION.SEPARATOR];
    tgtMList[mI] = m;
  }
  return tgtMList.filter((m) => m !== null).join(tgtSeparator);
}

function convertToHelm(srcUh: UnitsHandler, src: string): string {
  const isDna = src.startsWith('DNA');
  const isRna = src.startsWith('RNA');
  const [prefix, leftWrapper, rightWrapper, postfix] = srcUh.getHelmWrappers();
  const srcS = srcUh.getSplitter()(src);

  const tgtMList: string[] = wu(srcS).map((srcM: string) => {
    let m: string = srcM;
    if (srcUh.isGap(m))
      m = GapSymbols[NOTATION.HELM];
    else if (isDna || isRna)
      m = m.replace(helmWrappersRe, '$1');
    else
      m = srcM.length == 1 ? `${leftWrapper}${srcM}${rightWrapper}` : `${leftWrapper}[${srcM}]${rightWrapper}`;
    return m;
  }).toArray();
  return `${prefix}${tgtMList.join('.')}${postfix}`;
}

/** Splits Helm sequence adjusting nucleotides to single char symbols. (!) Removes lone phosphorus. */
function splitterAsHelmNucl(srcUh: UnitsHandler, src: string): string[] {
  const srcMList: ISeqSplitted = srcUh.getSplitter()(src);
  const tgtMList: (string | null)[] = new Array<string>(srcMList.length);
  const isDna = src.startsWith('DNA');
  const isRna = src.startsWith('RNA');
  for (const [srcM, mI] of wu.enumerate(srcMList)) {
    let m: string | null = srcM;
    if (isDna || isRna) {
      m = m.replace(helmWrappersRe, '$1');
      m = m === 'P' ? null : m;
    }
    tgtMList[mI] = m;
  }
  return tgtMList.filter((m) => m !== null) as string[];
}
