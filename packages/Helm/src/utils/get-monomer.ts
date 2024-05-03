import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as org from 'org';
import scil from 'scil';

import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';

import {_package} from '../package';


export type GetMonomerResType = org.helm.WebEditorMonomer | null;

export type GetMonomerFunc = (a: org.helm.IAtom | org.helm.HelmType, name: string) => GetMonomerResType;
type GetMonomerOverridingFunc = (
  a: org.helm.IAtom | org.helm.HelmType, name: string, monomerLib: IMonomerLib,
  originalGetMonomer: (a: org.helm.IAtom | org.helm.HelmType, name: string) => GetMonomerResType) => GetMonomerResType;

export function getMonomerOverrideAndLogAlert(
  monomerLib: IMonomerLib, getMonomerOverriding: GetMonomerOverridingFunc, trigger: () => void, logger: ILogger,
): void {
  const logPrefix = 'Helm: overrideGetMonomerAndAlert()';
  const monomers = org.helm.webeditor.Monomers;
  const getMonomerOriginal = monomers.getMonomer.bind(monomers);
  const alertOriginal = scil.Utils.alert;
  try {
    org.helm.webeditor.Monomers.getMonomer = (
      a: org.helm.IAtom | org.helm.HelmType, name: string
    ): GetMonomerResType => {
      return getMonomerOverriding(a, name, monomerLib, getMonomerOriginal);
    };
    // Preventing alert message box for missing monomers with compressed Scilligence.JSDraw2.Lite.js
    scil.Utils.alert = (s: string): void => {
      logger.warning(`${logPrefix}, scil.Utils.alert() s = 's'.`);
    };
    trigger();
  } finally {
    monomers.getMonomer = getMonomerOriginal;
    scil.Utils.alert = alertOriginal;
  }
}
