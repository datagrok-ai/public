import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  HelmType, HelmMol,
  JSDraw2ModuleType, OrgType
} from '@datagrok-libraries/bio/src/helm/types';


import {Chain} from './pt-conversion';
import {getAvailableMonomers} from './utils';
import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes, PolyToolPlaceholders} from './types';

// For example keep monomers presented in HELMCoreLibrary.json only (not [NH2])
export const PT_HELM_EXAMPLE = 'PEPTIDE1{R.[Aca].T.G.H.F.G.A.A.Y.P.E.[meI]}$$$$';

/** Initialized by getHelmHelper via init Helm package */
declare const JSDraw2: JSDraw2ModuleType;
declare const org: OrgType;

function polyToolEnumeratorCore(m: HelmMol, position: number, monomerList: string[]): HelmMol[] {
  const resMolList: HelmMol[] = new Array<HelmMol>(monomerList.length);
  for (let i = 0; i < monomerList.length; i++) {
    const newSymbol = monomerList[i];
    const resM = resMolList[i] = m.clone() as HelmMol;
    const oldSymbol = resM.atoms[position].elem;
    resM.atoms[position].elem = newSymbol;

    const idOldSymbol = oldSymbol?.length > 1 ? `[${oldSymbol}]` : oldSymbol;
    const idNewSymbol = newSymbol?.length > 1 ? `[${newSymbol}]` : newSymbol;
    resM.name = `${m.name}-${idOldSymbol}${position + 1}${idNewSymbol}`;
  }
  return resMolList;
}

/**
 * @param {string} helm Molecule string Helm format
 * @param  placeholders Placeholders by zero-based position key
 * @returns {string[]} List of enumerated molecules in Helm format
 */
function getPtEnumeratorSingle(m: HelmMol, placeholders: PolyToolPlaceholders): HelmMol[] {
  const coreResList: HelmMol[][] = Object.entries(placeholders)
    .map(([p, monomerList]: [string, string[]]) => polyToolEnumeratorCore(m, parseInt(p), monomerList));
  const resMolList = coreResList.reduce((acc, posList) => acc.concat(posList), []);
  return resMolList;
}

function getPtEnumeratorMatrix(m: HelmMol, placeholders: PolyToolPlaceholders): HelmMol[] {
  let resMolList = [m];
  for (const [p, monomerList] of Object.entries(placeholders)) {
    const pos: number = parseInt(p);
    const posResMolList: HelmMol[][] = resMolList.map((m: HelmMol) => polyToolEnumeratorCore(m, pos, monomerList));
    resMolList = posResMolList.reduce((acc, l) => acc.concat(l), []);
  }
  return resMolList;
}

export function getPtEnumeratorHelm(helm: string, id: string, params: PolyToolEnumeratorParams): [string, string][] {
  const molHandler = new JSDraw2.MolHandler<HelmType>();
  const plugin = new org.helm.webeditor.Plugin(molHandler);
  org.helm.webeditor.IO.parseHelm(plugin, helm, new JSDraw2.Point(0, 0), undefined);
  const m = molHandler.m;
  m.name = id;

  let resMolList: HelmMol[];
  switch (params.type) {
  case PolyToolEnumeratorTypes.Single: {
    resMolList = getPtEnumeratorSingle(molHandler.m, params.placeholders);
    break;
  }
  case PolyToolEnumeratorTypes.Matrix: {
    resMolList = getPtEnumeratorMatrix(molHandler.m, params.placeholders);
    break;
  }
  }

  if (params.keepOriginal)
    resMolList = [m, ...resMolList];

  const resList = resMolList.map<[string, string]>((m: HelmMol) => { return [org.helm.webeditor.IO.getHelm(m)!, m.name!]; });
  return resList;
}
