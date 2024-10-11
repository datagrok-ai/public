import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

import {ISeqSplitted, SeqColStats} from './types';
import {NOTATION} from './consts';
import {HelmType} from '../../helm/types';
import {CellRendererBackBase} from '../cell-renderer-back-base';

export const SeqTemps = new class {
  /** Column's temp slot name for a SeqHandler object */
  seqHandler = `seq-handler`;
  notationProvider = `seq-handler.notation-provider`;
}();
export type ConvertFunc = (src: string) => string;
export type JoinerFunc = (src: ISeqSplitted) => string;

export interface ISeqHandler {
  get column(): DG.Column<string>;

  get alphabet(): string;
  get notation(): NOTATION;
  get separator(): string | undefined;

  get aligned(): string;
  get units(): string;
  get defaultBiotype(): HelmType;
  get maxLength(): number;
  get length(): number;
  get defaultGapOriginal(): string;
  get stats(): SeqColStats;
  get joiner(): JoinerFunc;

  get posList(): string[];

  isFasta(): boolean;
  isMsa(): boolean;
  isHelm(): boolean;
  isSeparator(): boolean;


  getSplitted(rowIdx: number, limit?: number): ISeqSplitted;
  getHelm(rowIdx: number, options?: any): Promise<DG.SemanticValue<string>>;

  getAlphabetSize(): number;
  getAlphabetIsMultichar(): boolean;
  getNewColumnFromList(name: string, seqList: string[]): DG.Column<string>;
  convert(tgtNotation: NOTATION, tgtSeparator?: string): DG.Column<string>;
  convertHelmToFastaSeparator(srcSeq: string, tgtNotation: string,
    tgtSeparator?: string, tgtGapOriginal?: string): string;
  getRegion(startIdx: number | null, endIdx: number | null, name: string): DG.Column<string>;
  getJoiner(opts?: { notation: NOTATION, separator?: string }): JoinerFunc;

  getDistanceFunctionName(): MmDistanceFunctionsNames;
  getConverter(tgtUnits: NOTATION, tgtSeparator?: string): ConvertFunc;
  getRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>): CellRendererBackBase<string>;
}
