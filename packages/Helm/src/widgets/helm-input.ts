import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {fromEvent, Unsubscribable} from 'rxjs';

import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';
import {IHelmHelper, IInputInitOptions} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmMol, IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';

import {defaultErrorHandler} from '../utils/err-info';
import {getHoveredMonomerFromEditorMol} from '../utils/get-hovered';

import {_package} from '../package';


export class HelmInput extends DG.JsInputBase<HelmMol> {
  /** Input type identifier (such as "Slider" for the slider input). See {@link InputType}. */
  get inputType(): string {
    return 'Macromolecule';
  }

  /** Data type this input can edit. See {@link DG.Type}. */
  get dataType(): string {
    return DG.TYPE.STRING;
  }

  getInput(): HTMLElement {
    return this.viewerHost;
  }

  getValue(): HelmMol {
    return this.viewer.editor.m;
  }

  setValue(value: HelmMol): void {
    this.viewer.editor.setMol(value);
  }

  getStringValue(): string {
    return this.viewer.editor.getHelm();
  }

  setStringValue(value: string): void {
    this.viewer.editor.setHelm(value);
  }

  viewerHost: HTMLDivElement;
  viewer: IHelmWebEditor;
  helmString: string = '';
  helmSelection: number[];

  private subs: Unsubscribable[];

  protected constructor(
    private readonly helmHelper: IHelmHelper,
    private readonly monomerLib: IMonomerLib,
    name?: string,
    private readonly logger = _package.logger,
  ) {
    super();

    if (name) this.captionLabel.innerText = name;

    this.viewerHost = ui.div([], {
      classes: 'ui-input-editor',
      style: {width: '100%', height: '100%', overflow: 'hidden'},
    });
    this.viewer = this.helmHelper.createHelmWebEditor(this.viewerHost);

    this.subs = [];
    this.subs.push(this.monomerLib.onChanged
      .subscribe(() => this.viewer.editor.redraw()));
    this.subs.push(fromEvent<MouseEvent>(this.viewer.host, 'mousemove')
      .subscribe(this.viewerOnMouseMove.bind(this)));
    this.subs.push(fromEvent<MouseEvent>(this.viewer.host, 'click')
      .subscribe(this.viewerOnClick.bind(this)));

    this.root.classList.add('ui-input-helm');
    this.root.append(this.viewerHost);
  }

  protected toLog(): string {
    return `Helm: HelmInput<#>`;
  }

  /** Inspired by {@link DG.MarkdownInput}, {@link ui._create} */
  static create(
    helmHelper: IHelmHelper, monomerLib: IMonomerLib, name?: string, options?: IInputInitOptions<HelmMol>
  ): HelmInput {
    const input = new HelmInput(helmHelper, monomerLib, name);
    // TODO: Apply options
    const value: HelmMol | string | undefined = options?.value;
    if (value !== undefined) {
      if (value instanceof String || typeof value === 'string')
        input.stringValue = value as any as string;
      else
        input.value = value;
    }

    return input;
  }

  // static async init(host?: HTMLElement, cell?: DG.Cell): Promise<HelmInput> {
  //   const helmHelper = await getHelmHelper();
  //   const libHelper = await getMonomerLibHelper();
  //
  //   const editor = helmHelper.createHelmWebEditor();
  //   editor.host.style.width = '270px';
  //   editor.host.style.height = '150px';
  //   editor.host.style.paddingLeft = '40px';
  //
  //   cell = cell ?? grok.shell.tv.dataFrame.currentCell;
  //
  //   if (cell.column.semType === DG.SEMTYPE.MACROMOLECULE && cell.column.tags[DG.TAGS.UNITS] === NOTATION.HELM)
  //     editor.editor.setHelm(cell.value);
  //   else
  //     editor.editor.setHelm(PT_HELM_EXAMPLE);
  //
  //   return new HelmInput(helmHelper, editor, libHelper);
  // }

  setHelmString(helm: string): void {
    this.helmString = helm;
    this.viewer.editor.setHelm(helm);
    this.helmSelection = [];
  }

  getDiv(): HTMLDivElement {
    const title = ui.divText('Macromolecule', {style: {paddingTop: '43px', color: 'var(--grey-4)'}});

    return ui.divH([title, this.viewer.host], {style: {paddingLeft: '48px'}});
  }

  // -- Handle events --

  private viewerOnClick(_event: MouseEvent): void {
    ui.tooltip.hide();

    const webEditorHost = ui.div();
    const webEditorApp = this.helmHelper.createWebEditorApp(webEditorHost, this.viewer.editor.getHelm());

    const destroyEditorDialog = () => {
      $(webEditorHost).empty();
    };

    // Edit the molecule and select monomers to enumerate
    const dlg = ui.dialog({showHeader: false, showFooter: true})
      .add(webEditorHost!)
      .onOK(() => {
        try {
          const webEditorValue = webEditorApp.canvas!.getHelm(true)
            .replace(/<\/span>/g, '').replace(/<span style='background:#bbf;'>/g, '');
          this.viewer.editor.setHelm(webEditorValue);
          this.helmString = webEditorValue;

          const selection = webEditorApp.canvas!.helm!.jsd.m.atoms;
          this.helmSelection = [];
          for (let i = 0; i < selection.length; i++) {
            if (selection[i].selected)
              this.helmSelection.push(i);

            this.viewer.editor.m.atoms[i].highlighted = selection[i].selected;
          }
          this.viewer.editor.redraw();
        } catch (err: any) {
          defaultErrorHandler(err);
        } finally {
          destroyEditorDialog();
        }
      })
      .onCancel(() => {
        destroyEditorDialog();
      })
      .show({modal: true, fullScreen: true});
  }

  private viewerOnMouseMove(event: MouseEvent): void {
    const logPrefix = `${this.toLog()}.viewerOnMouseMove()`;
    const dpr = window.devicePixelRatio;

    const argsX = event.offsetX;
    const argsY = event.offsetY;
    const mol = this.viewer.editor.m;
    const seqMonomer = getHoveredMonomerFromEditorMol(argsX, argsY, mol, this.viewer.editor.div.clientHeight);
    if (seqMonomer) {
      const monomerLib = _package.monomerLib;
      const tooltipEl = monomerLib ? monomerLib.getTooltip(seqMonomer.polymerType, seqMonomer.symbol) :
        ui.divText('Monomer library is not available');
      ui.tooltip.show(tooltipEl, event.x + 16, event.y + 16);
    } else {
      // Tooltip for missing monomers
      ui.tooltip.hide();
    }
    this.logger.debug(`${logPrefix}, x = ${event.x}, y = ${event.y}`);
  }
}
