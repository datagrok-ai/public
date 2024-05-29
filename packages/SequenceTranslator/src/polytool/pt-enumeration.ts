import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {IHelmHelper, getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';
import '@datagrok-libraries/bio/src/types/helm';
import '@datagrok-libraries/bio/src/types/jsdraw2';
import * as org from 'org';
import $ from 'cash-dom';
import {Unsubscribable, fromEvent} from 'rxjs';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {IMonomerLibFileManager, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

const PT_HELM_EXAMPLE = 'PEPTIDE1{[R].[F].[T].[G].[H].[F].[G].[A].[A].[Y].[P].[E].[NH2]}$$$$';

export async function getLibrariesList(): Promise<string[]> {
  const monomerLib: IMonomerLibHelper = await grok.functions.call('Bio:getMonomerLibHelper', {});

  const monomerFileManager: IMonomerLibFileManager = await monomerLib.getFileManager();

  return monomerFileManager.getValidLibraryPaths();
}


export class HelmInput {
  webEditorHost: HTMLDivElement | null;
  webEditorApp: org.helm.IWebEditorApp | null;
  helmHelper: IHelmHelper;
  editor: IHelmWebEditor;
  subs: Unsubscribable[];

  helmString: string = PT_HELM_EXAMPLE;
  helmSelection: number[];

  constructor(
    helmHelper: IHelmHelper, editor: IHelmWebEditor) {
    this.helmHelper = helmHelper;
    this.editor = editor;

    this.webEditorHost = null;
    this.webEditorApp = null;
    const subs: Unsubscribable[] = [];

    subs.push(fromEvent<MouseEvent>(editor.host, 'click').subscribe(() => {
      this.webEditorHost = ui.div();
      //TODO: use not hh, but anything from editor.getHelm(true)
      this.webEditorApp = helmHelper.createWebEditorApp(this.webEditorHost, PT_HELM_EXAMPLE);
      const dlg = ui.dialog({showHeader: false, showFooter: true})
        .add(this.webEditorHost!)
        .onOK(() => {
          try {
            const webEditorValue = this.webEditorApp!.canvas.getHelm(true)
              .replace(/<\/span>/g, '').replace(/<span style='background:#bbf;'>/g, '');
            editor.editor.setHelm(webEditorValue);
            this.helmString = webEditorValue;
            this.helmSelection = [];

            //@ts-ignore
            const selection = this.webEditorApp?.canvas.helm.jsd.m.atoms;
            for (let i = 0; i < selection.length; i++) {
              if (selection[i].selected)
                this.helmSelection.push(i);
            }
          } catch (err: any) {
            //this.logger.error(err);
          } finally {
            $(this.webEditorHost).empty();
            this.webEditorHost = null;
            this.webEditorApp = null;
          }
        })
        .onCancel(() => {
          $(this.webEditorHost).empty();
          this.webEditorHost = null;
          this.webEditorApp = null;
        })
        .show({modal: true, fullScreen: true});
    }));
  }

  static async init() {
    const helmHelper = await getHelmHelper();

    const editor = helmHelper.createHelmWebEditor();
    editor.host.style.width = '300px';
    editor.host.style.height = '60px';
    editor.editor.setHelm(PT_HELM_EXAMPLE);

    return new HelmInput(helmHelper, editor);
  }

  getHelmString(): string {
    return this.helmString;
  }

  getHelmSelections() {
    return this.helmSelection;
  }

  getDiv(): HTMLDivElement {
    const title = ui.divText('Macromolecule');
    return ui.divH([title, this.editor.host]);
  }
}
