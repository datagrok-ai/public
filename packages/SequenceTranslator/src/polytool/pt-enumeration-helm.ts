import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  HelmType, HelmMol,
  JSDraw2ModuleType, OrgType,
  IHelmEditorOptions
} from '@datagrok-libraries/bio/src/helm/types';


import {Chain} from './pt-conversion';
import {getAvailableMonomers} from './utils';
import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes, PolyToolPlaceholders, PolyToolPlaceholdersBreadth} from './types';

// For example keep monomers presented in HELMCoreLibrary.json only (not [NH2])
export const PT_HELM_EXAMPLE = 'PEPTIDE1{R.[Aca].T.G.H.F.G.A.A.Y.P.E.[meI]}$$$$';

/** Initialized by getHelmHelper via init Helm package */
declare const JSDraw2: JSDraw2ModuleType;
declare const org: OrgType;

function polyToolEnumeratorCore(m: HelmMol, start: number, end: number, monomerList: string[]): HelmMol[] {
  const resMolList: HelmMol[] = new Array<HelmMol>(monomerList.length * (end - start + 1));
  for (let monI: number = 0; monI < monomerList.length; ++monI) {
    const posCount = end - start + 1;
    for (let posI: number = 0; posI < posCount; ++posI) {
      const pos = start + posI;
      const newSymbol = monomerList[monI];
      const resM = resMolList[monI * posCount + posI] = m.clone() as HelmMol;
      const oldSymbol = resM.atoms[pos].elem;
      resM.atoms[pos].elem = newSymbol;

      const idOldSymbol = oldSymbol?.length > 1 ? `[${oldSymbol}]` : oldSymbol;
      const idNewSymbol = newSymbol?.length > 1 ? `[${newSymbol}]` : newSymbol;
      resM.name = `${m.name}-${idOldSymbol}${pos + 1}${idNewSymbol}`;
    }
  }
  return resMolList;
}

/**
 * @param {string} helm Molecule string Helm format
 * @param  placeholders Placeholders by zero-based position key
 * @returns {string[]} List of enumerated molecules in Helm format
 */
function getPtEnumeratorSingle(m: HelmMol, placeholders: PolyToolPlaceholders): HelmMol[] {
  const coreResList: HelmMol[][] = placeholders
    .map((ph) => polyToolEnumeratorCore(m, ph.position, ph.position, ph.monomers));
  const resMolList = coreResList.reduce((acc, posList) => acc.concat(posList), []);
  return resMolList;
}

function getPtEnumeratorMatrix(m: HelmMol, placeholders: PolyToolPlaceholders): HelmMol[] {
  let resMolList = [m];
  for (const ph of placeholders) {
    const phResMolList: HelmMol[][] = resMolList.map((m: HelmMol) => polyToolEnumeratorCore(m, ph.position, ph.position, ph.monomers));
    resMolList = phResMolList.reduce((acc, l) => acc.concat(l), []);
  }
  return resMolList;
}

function getPtEnumeratorBreadth(m: HelmMol, placeholdersBreadth: PolyToolPlaceholdersBreadth): HelmMol[] {
  let resMolList = [m];
  for (const phb of placeholdersBreadth) {
    const phResMolList: HelmMol[][] = resMolList.map((m: HelmMol) => polyToolEnumeratorCore(m, phb.start, phb.end, phb.monomers));
    resMolList = phResMolList.reduce((acc, l) => acc.concat(l), []);
  }
  return resMolList;
}

/** Returns list of Helm with id. Covered with tests. */
export function doPolyToolEnumerateHelm(
  helm: string, id: string, params: PolyToolEnumeratorParams
): [ /* helm */ string, /* id */ string][] {
  const molHandler = new JSDraw2.MolHandler<HelmType, IHelmEditorOptions>();
  const plugin = new org.helm.webeditor.Plugin(molHandler);
  org.helm.webeditor.IO.parseHelm(plugin, helm, new JSDraw2.Point(0, 0), undefined);
  const m = molHandler.m;
  m.name = id;

  let resMolList: HelmMol[] = [];
  if (params.placeholders) {
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
  }

  let resBreadthMolList: HelmMol[] = [];
  if (params.placeholdersBreadth) {
    resBreadthMolList = getPtEnumeratorBreadth(molHandler.m, params.placeholdersBreadth);
  }
  resMolList = resMolList.concat(resBreadthMolList);

  if (params.keepOriginal)
    resMolList = [m, ...resMolList];

  const resList = resMolList.map<[string, string]>((m: HelmMol) => { return [org.helm.webeditor.IO.getHelm(m)!, m.name!]; });
  return resList;
}
