import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as org from 'org';
import * as scil from 'scil';
import * as JSDraw2 from 'JSDraw2';
import Atom = JSDraw2.Atom;

import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {IMonomerLib, HelmType} from '@datagrok-libraries/bio/src/types';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {HelmService} from './utils/helm-service';
import {GetMonomerFunc, GetMonomerResType, getWebEditorMonomer} from './utils/get-monomer';


type HelmWindowType = Window & {
  $helmService?: HelmServiceBase,
}
declare const window: HelmWindowType;

export function _getHelmService(): HelmServiceBase {
  let res = window.$helmService;
  if (!res) res = window.$helmService = new HelmService();
  return res;
}

export function initHelmPatchDojo(): void {
  // patch window.dojox.gfx.svg.Text.prototype.getTextWidth hangs
  /** get the text width in pixels */
  // @ts-ignore
  window.dojox.gfx.svg.Text.prototype.getTextWidth = function() {
    const rawNode = this.rawNode;
    const oldParent = rawNode.parentNode;
    const _measurementNode = rawNode.cloneNode(true);
    _measurementNode.style.visibility = 'hidden';

    // solution to the "orphan issue" in FF
    let _width = 0;
    const _text = _measurementNode.firstChild.nodeValue;
    oldParent.appendChild(_measurementNode);

    // solution to the "orphan issue" in Opera
    // (nodeValue == "" hangs firefox)
    if (_text != '') {
      let watchdogCounter = 100;
      while (!_width && --watchdogCounter > 0) { // <-- hangs
        //Yang: work around svgweb bug 417 -- http://code.google.com/p/svgweb/issues/detail?id=417
        if (_measurementNode.getBBox)
          _width = parseInt(_measurementNode.getBBox().width);
        else
          _width = 68;
      }
    }
    oldParent.removeChild(_measurementNode);
    return _width;
  };
}

export const helmJsonReplacer = (key: string, value: any): any => {
  switch (key) {
  case '_parent': {
    return `${value.toString()}`;
  }
  default:
    return value;
  }
};

export class HelmPackage extends DG.Package {
  public getMonomerOriginal: GetMonomerFunc;
  public alertOriginal: (s: string) => void;

  /** Patches Pistoia Monomers.getMonomer method to utilize DG Bio monomer Lib */
  public initHelmPatchPistoia(monomerLib: IMonomerLib, logger: ILogger): void {
    const logPrefix: string = 'Helm: initHelmPatchPistoia()';
    const monomers = org.helm.webeditor.Monomers;

    this.getMonomerOriginal = monomers.getMonomer.bind(monomers);
    this.logger.debug(`${logPrefix}, this.getMonomerOriginal stored`);

    this.alertOriginal = scil.Utils.alert;

    org.helm.webeditor.Monomers.getMonomer = (
      a: Atom<HelmType> | HelmType, name: string
    ): GetMonomerResType => {
      const logPrefixInt = `${logPrefix}, org.helm.webeditor.Monomers.getMonomer()`;
      try {
        // logger.debug(`${logPrefixInt}, a: ${JSON.stringify(a, helmJsonReplacer)}, name: '${name}'`);

        // Creates monomers in lib
        const dgWem = getWebEditorMonomer(monomerLib, a, name);

        // // Returns null for gap
        // const oWem = this.getMonomerOriginal(a, name);
        // if (!oWem)
        //   logger.warning(`${logPrefixInt}, getMonomerOriginal( a: ${a}, name: ${name}) returns null`);
        return dgWem; //dgWem;
      } catch (err) {
        const [errMsg, errStack] = errInfo(err);
        logger.error(`${logPrefixInt}, Error: ${errMsg}`, undefined, errStack);
        throw err;
      }
    };
    scil.Utils.alert = (s: string): void => {
      logger.warning(`${logPrefix}, scil.Utils.alert() s = 's'.`);
    };

    // @ts-ignore, intercept with proxy to observe access and usage
    org.helm.webeditor.Monomers = new class {
      constructor(base: org.helm.IMonomers) {
        return new Proxy(base, {
          get(target: any, p: string | symbol, _receiver: any): any {
            return target[p];
          },
          set(target: any, p: string | symbol, newValue: any, _receiver: any): boolean {
            if (p != 'sugars' && p != 'linkers' && p != 'bases' && p != 'aas') {
              const k = 11;
            }
            target[p] = newValue;
            return true;
          },
          apply(target: any, thisArg: any, argArray: any[]): any {
            return target[thisArg](...argArray);
          },
        });
      }
    }(org.helm.webeditor.Monomers) as org.helm.IMonomers;

    // @ts-ignore, intercept with proxy to observe access and usage
    org.helm.webeditor = new class {
      constructor(base: org.helm.IOrgHelmWebEditor) {
        return new Proxy(base, {
          get(target: any, p: string | symbol, _receiver: any): any {
            return target[p];
          },
          set(target: any, p: string | symbol, newValue: any, _receiver: any): boolean {
            target[p] = newValue;
            return true;
          },
          apply(target: any, thisArg: any, argArray: any[]): any {
            return target[thisArg](...argArray);
          },
        });
      }
    }(org.helm.webeditor) as org.helm.IOrgHelmWebEditor;
  }
}
