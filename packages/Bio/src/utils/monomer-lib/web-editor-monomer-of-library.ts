import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HelmType, IMonomerColors, IWebEditorMonomer, MonomerType, PolymerType, WebEditorRGroups} from '@datagrok-libraries/bio/src/helm/types';
import {IMonomerLibBase, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {HELM_OPTIONAL_FIELDS as OPT, HELM_REQUIRED_FIELD as REQ, HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';

import {BrokenWebEditorMonomer, MissingWebEditorMonomer} from './web-editor-monomer-dummy';
import {naturalMonomerColors} from './monomer-colors';

export class LibraryWebEditorMonomer implements IWebEditorMonomer {
  public get rs(): number { return Object.keys(this.at).length; }

  public get issmiles(): boolean { return !!this.smiles; }

  public linecolor?: string;
  public backgroundcolor?: string;
  public textcolor?: string;

  /* eslint-disable max-params */
  protected constructor(
    public readonly id: string,
    public readonly m: string,
    public readonly n: string,
    public readonly na: string | undefined,
    public readonly type: PolymerType,
    public readonly mt: MonomerType,
    public readonly at: WebEditorRGroups,
    public readonly smiles?: string,
  ) /* eslint-enable max-params */ {}

  static fromMonomer(biotype: HelmType, monomer: Monomer, monomerLib: IMonomerLibBase): IWebEditorMonomer {
    let at: WebEditorRGroups = {};
    const symbol = monomer[REQ.SYMBOL];
    const smiles = monomer[REQ.SMILES];
    if (monomer.rgroups.length > 0) {
      monomer.rgroups.forEach((it) => {
        at[it[RGP.LABEL]] = it[RGP.CAP_GROUP_NAME];
      });
    } else if (smiles) {
      // Generate R-Groups from SMILES
      at = monomerLib.getRS(smiles);
    } else if (!monomer.lib) {
      // missing
      return new MissingWebEditorMonomer(biotype, symbol);
    } else {
      // broken
      return new BrokenWebEditorMonomer(biotype, symbol);
    }

    const res = new LibraryWebEditorMonomer(
      monomer[REQ.SYMBOL],
      monomer[REQ.MOLFILE],
      monomer[REQ.NAME],
      monomer[OPT.NATURAL_ANALOG],
      monomer[REQ.POLYMER_TYPE],
      monomer[REQ.MONOMER_TYPE],
      at);

    const colors = getMonomerColors(biotype, monomer);
    if (colors) {
      res.textcolor = colors?.textcolor;
      res.linecolor = colors?.linecolor;
      res.backgroundcolor = colors?.backgroundcolor;
    }

    return res;
  }
}

function getMonomerColors(biotype: HelmType, monomer: Monomer): IMonomerColors | null {
  const currentMonomerSchema = 'default';
  let monomerSchema: string = currentMonomerSchema;

  let res: any;
  if (monomer.meta && monomer.meta.colors) {
    const monomerColors: { [colorSchemaName: string]: any } = monomer.meta.colors;
    if (!(currentMonomerSchema in monomerColors)) monomerSchema = 'default';
    let res = monomerColors[monomerSchema];
  }

  if (!res) {
    const biotypeColors: { [symbol: string]: string } | undefined = naturalMonomerColors[biotype];
    const nColor = biotypeColors?.[monomer.symbol];
    if (nColor)
      res = {textColor: "#000000", lineColor: "#000000", backgroundColor: nColor};
  }

  const na = monomer[OPT.NATURAL_ANALOG];
  if (!res && na) {
    const biotypeColors: { [symbol: string]: string } | undefined = naturalMonomerColors[biotype];
    const naColor = biotypeColors?.[na];
    res = {textColor: "#000000", lineColor: "#000000", backgroundColor: naColor ?? "#FFFFFF",};
  }

  return !res ? null : {
    textcolor: res.text ?? res.textColor,
    linecolor: res.line ?? res.lineColor,
    backgroundcolor: res.background ?? res.backgroundColor
  } as IMonomerColors;
}
