import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';

import {HelmAtom} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';

import {
  HelmType, IWebEditorMonomer, MonomerSetType, MonomerType, PolymerType
} from '../helm/types';
import {
  HELM_REQUIRED_FIELD as REQ,
  HELM_RGROUP_FIELDS as RGP, HELM_OPTIONAL_FIELDS as OPT, HELM_POLYMER_TYPE,
} from '../utils/const';

export type RGroup = {
  [RGP.CAP_GROUP_SMILES]: string,
  [RGP.ALTERNATE_ID]: string,
  [RGP.CAP_GROUP_NAME]: string,
  [RGP.LABEL]: string,
}

/** https://github.com/PistoiaHELM/HELMMonomerSets/blob/master/HELMmonomerSchema.json */
export type Monomer = {
  [REQ.SYMBOL]: string,
  [REQ.NAME]: string,
  [REQ.MOLFILE]: string,
  [REQ.AUTHOR]: string,
  [REQ.ID]: number,
  [REQ.RGROUPS]: RGroup[],
  [REQ.SMILES]: string,
  [REQ.POLYMER_TYPE]: PolymerType,
  [REQ.MONOMER_TYPE]: MonomerType,
  [REQ.CREATE_DATE]: string | null,
  [OPT.NATURAL_ANALOG]?: string,
  [OPT.META]?: { [property: string]: any },

  lib?: IMonomerLib,
  wem?: IWebEditorMonomer,
};

export interface IMonomerLinkData {
  source: string;
  symbol: string;
  notes: string;
}

export interface IMonomerLink {
  get lib(): IMonomerLib;
  get symbol(): string;
}

export interface IMonomerSetPlaceholder {
  get symbol(): string;
  get polymerType(): PolymerType;
  get monomerType(): MonomerType;
  get monomerLinks(): IMonomerLinkData[];

  get monomers(): Monomer[];
}

export interface IMonomerSet {
  get source(): string | undefined;
  get error(): string | undefined;

  get description(): string;
  get placeholders(): IMonomerSetPlaceholder[];
}

export type MonomerLibSummaryType = { [polymerType: string]: number };

export type MonomerLibData = { [polymerType: string]: { [symbol: string]: Monomer } };

export interface IMonomerLibBase {
  get onChanged(): Observable<any>;

  /* Gets library monomer for sequence monomer */
  addMissingMonomer(polymerType: PolymerType, monomerSymbol: string): Monomer;

  getMonomer(polymerType: PolymerType | null, monomerSymbol: string): Monomer | null;

  /** HELMWebEditor expects null for HelmTypes.LINKER and R-Group count != 2 */
  getWebEditorMonomer(a: HelmAtom | HelmType, symbol?: string): IWebEditorMonomer | null;

  getRS(smiles: string): { [r: string]: string };
}

export interface IMonomerLib extends IMonomerLibBase {
  get source(): string | undefined;
  get error(): string | undefined;

  getMonomerMolsByPolymerType(polymerType: PolymerType): { [monomerSymbol: string]: string } | null;
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: PolymerType, element?: string): string[];
  getMonomerSymbolsByType(polymerType: PolymerType): string[];
  getPolymerTypes(): PolymerType[];
  update(lib: IMonomerLib): void;
  toJSON(): Monomer[];

  /** Summary string with lib monomer count by type
   * @deprecated Keep for backward compatibility */
  getSummary(): string;

  /** Summary with lib monomer count by type */
  getSummaryObj(): MonomerLibSummaryType;

  /** Gets dataframe with columns 'polymerType', 'count'. */
  getSummaryDf(): DG.DataFrame;

  getTooltip(biotype: HelmType, monomerSymbol: string): HTMLElement;

  // For monomer palettes
  getMonomerSet(biotype: HelmType): MonomerSetType | null;

  override(overrideData: MonomerLibData): IMonomerLibBase;
}
