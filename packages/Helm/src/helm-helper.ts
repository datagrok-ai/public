import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor, IWebEditorApp} from '@datagrok-libraries/bio/src/helm/types';
import {HelmWebEditor} from './helm-web-editor';

type HelmHelperWindowType = Window & {
  $helmHelper?: HelmHelper,
}
declare const window: HelmHelperWindowType;

export class HelmHelper implements IHelmHelper {
  createHelmWebEditor(): IHelmWebEditor {
    return new HelmWebEditor();
  }

  createWebEditorApp(host: HTMLDivElement, helm: string): IWebEditorApp {
    // @ts-ignore
    org.helm.webeditor.MolViewer.molscale = 0.8;
    // @ts-ignore
    const webEditorApp = new scil.helm.App(host, {
      showabout: false,
      mexfontsize: '90%',
      mexrnapinontab: true,
      topmargin: 20,
      mexmonomerstab: true,
      sequenceviewonly: false,
      mexfavoritefirst: true,
      mexfilter: true
    }) as IWebEditorApp;
    const sizes = webEditorApp.calculateSizes();
    webEditorApp.canvas.resize(sizes.rightwidth - 100, sizes.topheight - 210);
    let s = {width: sizes.rightwidth - 100 + 'px', height: sizes.bottomheight + 'px'};
    //@ts-ignore
    scil.apply(webEditorApp.sequence.style, s);
    //@ts-ignore
    scil.apply(webEditorApp.notation.style, s);
    s = {width: sizes.rightwidth + 'px', height: (sizes.bottomheight + webEditorApp.toolbarheight) + 'px'};
    //@ts-ignore
    scil.apply(webEditorApp.properties.parent.style, s);
    webEditorApp.structureview.resize(sizes.rightwidth, sizes.bottomheight + webEditorApp.toolbarheight);
    webEditorApp.mex.resize(sizes.topheight - 80);
    setTimeout(() => {
      webEditorApp.canvas.helm.setSequence(helm, 'HELM');
    }, 200);
    return webEditorApp;
  }

  // -- Instance singleton --

  public static async getInstance(): Promise<IHelmHelper> {
    let res: HelmHelper = window.$helmHelper!;
    if (!res) {
      // create & init singleton
      window.$helmHelper = res = new HelmHelper();
      // await res.init();
    }
    return res;
  }
}
