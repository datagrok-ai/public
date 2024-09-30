import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {fromEvent, Observable, Unsubscribable} from 'rxjs';

import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/index';
import {HelmInputBase, IHelmHelper, IHelmInputInitOptions} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmAtom, HelmMol, HelmString, IHelmWebEditor} from '@datagrok-libraries/bio/src/helm/types';

import {defaultErrorHandler} from '../utils/err-info';
import {getHoveredMonomerFromEditorMol, getSeqMonomerFromHelmAtom} from '../utils/get-hovered';

import {MonomerNumberingTypes} from '@datagrok-libraries/bio/src/helm/consts';

import {_package} from '../package';

export class HelmInput extends HelmInputBase {
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

  getValue(): HelmString {
    return this.viewer.editor.getHelm();
  }

  setValue(value: HelmString): void {
    this.viewer.editor.setHelm(value);
  }

  getStringValue(): string {
    return this.viewer.editor.getHelm();
  }

  setStringValue(value: string): void {
    this.viewer.editor.setHelm(value);
  }

  // -- IHelmInput --

  get molValue(): HelmMol {
    return this.viewer.editor.m;
  }

  set molValue(value: HelmMol) {
    this.viewer.editor.setMol(value);
  }

  get onMouseMove(): Observable<MouseEvent> { return fromEvent<MouseEvent>(this.viewer.host, 'mousemove'); }

  get onClick(): Observable<MouseEvent> { return fromEvent<MouseEvent>(this.viewer.host, 'click'); }

  readonly viewerHost: HTMLDivElement;
  readonly viewer: IHelmWebEditor;
  readonly editHintDiv: HTMLDivElement;

  helmString: string = '';
  helmSelection: number[];

  private subs: Unsubscribable[];

  protected constructor(
    private readonly helmHelper: IHelmHelper,
    private readonly monomerLib: IMonomerLibBase,
    name?: string,
    public readonly options?: IHelmInputInitOptions,
    private readonly logger = _package.logger,
  ) {
    super();

    if (name) this.captionLabel.innerText = name;

    this.viewerHost = ui.div([], {
      classes: 'ui-input-editor',
      style: {width: '100%', height: '100%', overflow: 'hidden'},
    });
    this.viewer = this.helmHelper.createHelmWebEditor(this.viewerHost, {
      monomerNumbering: MonomerNumberingTypes.continuous
    });

    this.editHintDiv = ui.divH([
      ui.link('Click to edit', () => { this.viewer.host.click(); }, undefined, {}),
    ], {style: {display: 'none', position: 'absolute', top: '0', right: '0'}});

    this.subs = [];
    this.subs.push(this.monomerLib.onChanged
      .subscribe(() => this.viewer.editor.redraw()));
    /* eslint-disable rxjs/no-ignored-subscription */
    fromEvent<MouseEvent>(this.viewer.host, 'mousemove')
      .subscribe(this.viewerOnMouseMove.bind(this));
    fromEvent<MouseEvent>(this.viewer.host, 'mouseenter')
      .subscribe(this.viewerOnMouseEnter.bind(this));
    fromEvent<MouseEvent>(this.viewer.host, 'mouseleave')
      .subscribe(this.viewerOnMouseLeave.bind(this));
    fromEvent<MouseEvent>(this.viewer.host, 'click')
      .subscribe(this.viewerOnClick.bind(this));
    /* eslint-enable rxjs/no-ignored-subscription */

    this.subs.push(ui.onSizeChanged(this.root).subscribe(() => {
      const rootStyle = window.getComputedStyle(this.root);
      const h = this.root.clientHeight - parseFloat(rootStyle.paddingTop) - parseFloat(rootStyle.paddingBottom);
      this.input.style.height = `${h}px`;
    }));

    this.subs.push(ui.onSizeChanged(this.input).subscribe(() => {
      const rootStyle = window.getComputedStyle(this.root);
      const h = this.root.clientHeight - parseFloat(rootStyle.paddingTop) - parseFloat(rootStyle.paddingBottom);
      this.input.style.height = `{h}px`;
      const w = this.input.clientWidth;
      this.viewer.editor.setSize(w, h);
    }));

    this.root.classList.add('ui-input-helm');
    this.root.append(this.viewerHost, this.editHintDiv);
  }

  detach(): void {
    for (const sub of this.subs)
      sub.unsubscribe();
  }

  protected toLog(): string {
    return `Helm: HelmInput<#>`;
  }

  /** Inspired by {@link DG.MarkdownInput}, {@link ui._create} */
  static create(helmHelper: IHelmHelper, monomerLib: IMonomerLibBase,
    name?: string, options?: IHelmInputInitOptions
  ): HelmInput {
    const input = new HelmInput(helmHelper, monomerLib, name, options);
    // TODO: Apply options
    const value: HelmMol | string | undefined = options?.value;
    if (value !== undefined) {
      if (value instanceof String || typeof value === 'string')
        input.stringValue = value as any as string;
      else if (typeof value === 'object' /* && value.T === 'MOL'*/)
        input.molValue = value;
      else
        throw new Error(`Unsupported value of type '${typeof value}'.`);
    }

    return input;
  }

  redraw(): void {
    this.viewer.editor.redraw();
  }

  showTooltip(content: HTMLElement | string, a: HelmAtom): void {
    /** Monomer sign offset for tooltip inspired by {@link org.helm.webeditor.drawMonomer()} */
    const so = this.viewer.editor.getDrawOptions().fontsize * 1.2;
    const bcr = this.input.getBoundingClientRect();
    ui.tooltip.show(content, bcr.left + a.p.x + so, bcr.top + a.p.y + so);
  }

  // -- Handle events --

  private viewerOnClick(_event: MouseEvent): void {
    ui.tooltip.hide();

    if (this.options?.editable != false)
      this.showEditorDialog();
  }

  private viewerOnMouseMove(event: MouseEvent): void {
    const logPrefix = `${this.toLog()}.viewerOnMouseMove()`;

    const argsX = event.offsetX;
    const argsY = event.offsetY;
    const mol = this.viewer.editor.m;
    const hoveredAtom = getHoveredMonomerFromEditorMol(argsX, argsY, mol, this.viewer.editor.div.clientHeight);
    if (hoveredAtom) {
      const seqMonomer = getSeqMonomerFromHelmAtom(hoveredAtom);
      const monomerLib = _package.monomerLib;
      const tooltipEl = monomerLib ? monomerLib.getTooltip(seqMonomer.biotype, seqMonomer.symbol) :
        ui.divText('Monomer library is not available');
      this.showTooltip(tooltipEl, hoveredAtom);
      event.preventDefault();
      event.stopPropagation();
    } else
      ui.tooltip.hide();

    //this.logger.debug(`${logPrefix}, x = ${event.x}, y = ${event.y}`);
  }

  private viewerOnMouseEnter(_event: MouseEvent): void {
    //this.logger.debug(`${this.toLog()}.viewerOnMouseEnter()`);
    if (this.options?.editable != false)
      this.showEditHint();
  }

  private viewerOnMouseLeave(event: MouseEvent): void {
    if (this.editHintDiv.contains(event.relatedTarget as Node)) return; // prevents flickering

    //this.logger.debug(`${this.toLog()}.viewerOnMouseLeave()`);
    if (this.options?.editable != false)
      this.hideEditHint();
  }

  // editHint

  showEditHint(): void {
    this.editHintDiv.style.removeProperty('display');
  }

  hideEditHint(): void {
    this.editHintDiv.style.display = 'none';
  }

  // editor dialog

  showEditorDialog(): void {
    const webEditorHost = ui.div();
    const webEditorApp = this.helmHelper.createWebEditorApp(webEditorHost, this.viewer.editor.getHelm());

    const destroyEditorDialog = () => {
      $(webEditorHost).empty();
    };

    // Edit the molecule and select monomers to enumerate
    const _editorDialog = ui.dialog({showHeader: false, showFooter: true})
      .add(webEditorHost!)
      .onOK(() => {
        try {
          const webEditorValue = webEditorApp!.canvas!.getHelm(true)
            .replace(/<\/span>/g, '').replace(/<span style='background:#bbf;'>/g, '');
          this.viewer.editor.setHelm(webEditorValue);
          this.helmString = webEditorValue;

          const selection = webEditorApp!.canvas!.helm!.jsd.m.atoms;
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
}
