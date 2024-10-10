import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';

import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {
  App, Editor, HelmType, IHelmWebEditor, HelmEditor, IHelmEditorOptions
} from '@datagrok-libraries/bio/src/helm/types';

import {JSDraw2Module, OrgHelmModule, ScilModule} from './types';

import {_package} from './package';

declare const scil: ScilModule;
declare const JSDraw2: JSDraw2Module;
declare const org: OrgHelmModule;

export class HelmWebEditor implements IHelmWebEditor {
  editor: HelmEditor;
  host: HTMLDivElement;

  w = 200;
  h = 100;
  private subs: Unsubscribable[];

  constructor(
    host?: HTMLDivElement,
    options?: Partial<IHelmEditorOptions>,
    private logger: ILogger = _package.logger
  ) {
    this.host = host ?? ui.div([], {style: {width: `${this.w}px`, height: `${this.h}px`}});

    const styleBackup = {
      width: this.host.style.width, height: this.host.style.height,
      overflow: this.host.style.overflow
    };
    this.editor = new JSDraw2.Editor(this.host,
      Object.assign({width: this.w, height: this.h, viewonly: true}, options));
    this.host.style.width = styleBackup.width;
    this.host.style.height = styleBackup.height;
    this.host.style.overflow = styleBackup.overflow;

    this.subs = [];
    this.subs.push(ui.onSizeChanged(this.host).subscribe(this.hostOnSizeChanged.bind(this)));
  }

  protected toLog(): string {
    return `Helm: HelmWebEditor<#>`;
  }

  private lastSize: { width: number, height: number } = {width: -1, height: -1};

  private hostOnSizeChanged(value: any): void {
    const size = {width: value.target.clientWidth, height: value.target.clientHeight};
    if (this.lastSize.width != size.width || this.lastSize.height != size.height) {
      const logPrefix = `${this.toLog()}.hostOnSizeChanged()`;
      this.logger.debug(`${logPrefix}`);
      this.editor.setSize(this.host.clientWidth, this.host.clientHeight);
    }
    this.lastSize = size;
  }

  resizeEditor(w: number, h: number) {
    this.editor.setSize(w, h);
  }
}
