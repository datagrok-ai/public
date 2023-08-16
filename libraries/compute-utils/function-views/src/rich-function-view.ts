/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import html2canvas from 'html2canvas';
import wu from 'wu';
import $ from 'cash-dom';
import {Subject, BehaviorSubject, Observable} from 'rxjs';
import '../css/rich-function-view.css';
import {UiUtils} from '../../shared-components';
import {FunctionView} from './function-view';
import {startWith} from 'rxjs/operators';
import {EXPERIMENTAL_TAG, viewerTypesMapping} from './shared/consts';
import {boundImportFunction, getFuncRunLabel, getPropViewers} from './shared/utils';
import {FuncCallInput, SubscriptionLike} from '../../shared-components/src/FuncCallInput';

const FILE_INPUT_TYPE = 'file';

export type InputVariants = DG.InputBase | FuncCallInput;

function isInputBase(input: FuncCallInput): input is DG.InputBase {
  const inputAny = input as any;
  return (inputAny.dart && DG.toJs(inputAny.dart) instanceof DG.InputBase);
}

function getObservable<T>(onInput: (f: Function) => SubscriptionLike): Observable<T> {
  return new Observable((observer: any) => {
    const sub = onInput((val: T) => {
      observer.next(val);
    });
    return () => sub.unsubscribe();
  });
}

export interface AfterInputRenderPayload {
  prop: DG.Property;
  input: InputVariants;
}

export interface AfterOutputRenderPayload {
  prop: DG.Property;
  output: DG.Viewer;
}

enum SYNC_FIELD {
  INPUTS = 'inputs',
  OUTPUTS = 'outputs'
}

type SyncFields = SYNC_FIELD.INPUTS | SYNC_FIELD.OUTPUTS;
const syncParams = {
  [SYNC_FIELD.INPUTS]: 'inputParams',
  [SYNC_FIELD.OUTPUTS]: 'outputParams',
} as const;
export class RichFunctionView extends FunctionView {
  // emitted when runButton disability should be checked
  private checkDisability = new Subject();

  // stores the running state
  private isRunning = new BehaviorSubject(false);

  // stores simulation or upload mode flag
  private isUploadMode = new BehaviorSubject<boolean>(false);
  private inputsOverride: Record<string, FuncCallInput> = {};
  private inputsMap: Record<string, FuncCallInput> = {};

  static fromFuncCall(
    funcCall: DG.FuncCall,
    options: {historyEnabled: boolean, isTabbed: boolean} =
    {historyEnabled: true, isTabbed: false},
  ) {
    return new this(funcCall, options);
  }

  constructor(
    initValue: string | DG.FuncCall,
    public options: { historyEnabled: boolean, isTabbed: boolean} =
    {historyEnabled: true, isTabbed: false},
  ) {
    super(initValue, options);
  }

  protected async onFuncCallReady() {
    await this.loadInputsOverrides();
    await super.onFuncCallReady();
    this.basePath = `scripts/${this.funcCall.func.id}/view`;

    if (this.runningOnStart && this.isRunnable()) await this.doRun();
  }

  protected prevOpenedTab = null as DG.TabPane | null;
  /**
   * Saving previously opened tab
   * @param runFunc
   */
  public override onBeforeRun(): Promise<void> {
    this.prevOpenedTab = Object.keys(this.categoryToDfParamMap.inputs).includes(this.outputsTabsElem.currentPane.name) ? null: this.outputsTabsElem.currentPane;

    return Promise.resolve();
  }

  /**
   * Showing UI after completion of function call.
   * @param runFunc
   */
  public override onAfterRun(): Promise<void> {
    this.showOutputTabsElem();

    if (this.prevOpenedTab) {
      this.outputsTabsElem.currentPane = this.prevOpenedTab;
      return Promise.resolve();
    }

    const firstOutputTab = this.outputsTabsElem.panes
      .find((tab) => Object.keys(this.categoryToDfParamMap.outputs).includes(tab.name));
    if (firstOutputTab) this.outputsTabsElem.currentPane = firstOutputTab;

    return Promise.resolve();
  }

  // scripting api events
  // regular and experimental inputs
  public beforeInputPropertyRender = new Subject<DG.Property>();
  public afterInputPropertyRender = new Subject<AfterInputRenderPayload>();
  public afterOutputPropertyRender = new Subject<AfterOutputRenderPayload>();
  // output scalars table
  public afterOutputSacalarTableRender = new Subject<HTMLElement>();

  public getRunButton(name = 'Run') {
    const runButton = ui.bigButton(getFuncRunLabel(this.func) ?? name, async () => await this.doRun());
    const disabilitySub = this.checkDisability.subscribe(() => {
      const isValid = this.isRunnable();
      runButton.disabled = !isValid;
    });
    this.subs.push(disabilitySub);

    return runButton;
  }

  public async loadInputsOverrides() {
    const inputParams = [...this.funcCall.inputParams.values()] as DG.FuncCallParam[];
    await Promise.all(inputParams.map(async (param) => {
      if (param.property.options.input) {
        const func: DG.Func = await grok.functions.eval(param.property.options.input);
        const call = func.prepare({params: JSON.parse(param.property.options.inputOptions || '{}')});
        await call.call();
        this.inputsOverride[param.name] = call.outputs.input;
      }
    }));
  }

  private getSaveButton(name = 'Save') {
    const saveButton = ui.bigButton(name, async () => await this.saveExperimentalRun(this.funcCall), 'Save uploaded data');

    this.isUploadMode.subscribe((newValue) => {
      this.buildRibbonPanels();
      if (this.runningOnInput) return;

      if (newValue)
        $(saveButton).show();
      else
        $(saveButton).hide();
    });
    if (!this.runningOnInput || this.options.isTabbed) $(saveButton).hide();

    return saveButton;
  }

  private getStandardButtons(): HTMLElement[] {
    const runButton = this.getRunButton();
    const runButtonWrapper = ui.div([runButton]);
    ui.tooltip.bind(runButtonWrapper, () => runButton.disabled ? (this.isRunning.value ? 'Computations are in progress' : 'Some inputs are invalid') : '');
    const saveButton = this.getSaveButton();

    if (this.runningOnInput) $(runButtonWrapper).hide();

    return [saveButton, runButtonWrapper];
  }

  /**
   * Override to change additional buttons placed between navigation and run buttons.
   */
  protected additionalBtns = ui.divH([]) as HTMLElement;
  /**
   * Changes additional buttons to provided ones.
   * @param additionalBtns Array of HTML elements to place instead of the existing additional buttons.
   */
  public setAdditionalButtons(additionalBtns: HTMLElement[]) {
    const additionalBtnsContainer = ui.divH(additionalBtns);
    this.additionalBtns.replaceWith(additionalBtnsContainer);
    this.additionalBtns = additionalBtnsContainer;
  }

  /**
   * Override to change navigation buttons placed next to the additional buttons.
   */
  protected navBtns = ui.divH([]) as HTMLElement;
  /**
   * Changes navigation buttons to provided ones.
   * @param navBtns Array of HTML elements to place instead of the existing navigation buttons.
   */
  public setNavigationButtons(navBtns: HTMLElement[]) {
    const navBtnsContainer = ui.divH(navBtns);
    this.navBtns.replaceWith(navBtnsContainer);
    this.navBtns = navBtnsContainer;
  }

  /**
   * RichFunctionView has advanced automatic UI builder. It takes {@link this.funcCall} as a base and constructs flexible view.
   * This view is updated automatically when {@link this.funcCallReplaced} is emitted or any of input/output param changes.
   * @returns HTMLElement attached to the root of the view
   */
  public buildIO(): HTMLElement {
    const {inputBlock, inputForm, outputForm, controlsWrapper} = this.buildInputBlock();

    ui.tools.handleResize(inputBlock, () => {
      if (([
        ...Array.from(inputForm.childNodes),
        ...this.isUploadMode.value ? [Array.from(outputForm.childNodes)]: [],
      ]).some((child) => $(child).width() < 250) ||
      $(inputBlock).width() < 350) {
        $(inputForm).addClass('ui-form-condensed');
        $(outputForm).addClass('ui-form-condensed');
        $(controlsWrapper).addClass('ui-form-condensed');
      } else {
        $(inputForm).removeClass('ui-form-condensed');
        $(outputForm).removeClass('ui-form-condensed');
        $(controlsWrapper).removeClass('ui-form-condensed');
      }
    });

    const outputBlock = this.buildOutputBlock();
    outputBlock.style.height = '100%';
    outputBlock.style.width = '100%';
    $(this.outputsTabsElem.root).hide();

    if (Object.keys(this.categoryToDfParamMap.inputs).length > 0) {
      this.outputsTabsElem.panes.forEach((tab) => {
        $(tab.header).hide();
      });
    }

    const out = ui.splitH([inputBlock, ui.panel([outputBlock], {style: {'padding-top': '0px'}})], null, true);
    out.style.padding = '0 12px';

    inputBlock.style.maxWidth = '450px';

    return out;
  }

  public buildInputBlock() {
    const inputFormDiv = this.renderInputForm();
    const outputFormDiv = this.renderOutputForm();
    const standardButtons = this.getStandardButtons();

    const controllsDiv = ui.buttonsInput([
      this.navBtns as any,
      ui.divH([
        this.additionalBtns,
        ...standardButtons,
      ], {style: {'gap': '5px'}}),
    ]);
    $(controllsDiv).addClass('rfv-buttons');

    const controlsForm = ui.div(controllsDiv, 'ui-form ui-form-wide');
    $(controlsForm).css({
      'padding-left': '0px',
      'padding-bottom': '0px',
      'max-width': '100%',
      'min-height': '50px',
    });

    const experimentalDataSwitch = ui.switchInput('', this.isUploadMode.value, (v: boolean) => this.isUploadMode.next(v));
    this.isUploadMode.subscribe((newValue) => {
      experimentalDataSwitch.notify = false;
      experimentalDataSwitch.value = newValue,
      experimentalDataSwitch.notify = true;
    });

    const form = ui.divV([
      inputFormDiv,
      ...this.hasUploadMode ? [
        ui.divH([ui.h2('Experimental data'), experimentalDataSwitch.root], {style: {'flex-grow': '0'}}),
        outputFormDiv,
      ]: [],
      controlsForm,
    ], 'ui-box rfv-form');

    this.isUploadMode.subscribe((newValue) => {
      if (newValue)
        $(outputFormDiv).show();
      else
        $(outputFormDiv).hide();
    });

    return {
      inputBlock: form,
      inputForm: inputFormDiv,
      outputForm: outputFormDiv,
      controlsWrapper: controlsForm,
    };
  }

  buildRibbonPanels(): HTMLElement[][] {
    super.buildRibbonPanels();

    const play = ui.iconFA('play', async () => await this.doRun(), 'Run computations');
    play.classList.add('fas');

    const save = ui.iconFA('save', async () => {
      if (this.isUploadMode.value) {
        await this.saveExperimentalRun(this.funcCall);
        return;
      }

      if (this.lastCall)
        await this.saveRun(this.lastCall);
      else
        grok.shell.warning('Function was not called. Call it before saving');
    }, this.isUploadMode.value ? 'Save uploaded data': 'Save the last run');

    const toggleUploadMode = ui.iconFA('arrow-to-top', async () => {
      this.isUploadMode.next(!this.isUploadMode.value);

      if (boundImportFunction(this.func)) {
        const func = await grok.functions.eval(boundImportFunction(this.func)!) as DG.Func;
        func.prepare().edit();
        return;
      }

      toggleUploadMode.classList.toggle('d4-current');
    }, 'Upload experimental data');
    toggleUploadMode.classList.add(
      'd4-toggle-button',
      ...this.isUploadMode.value ? ['d4-current']: [],
    );

    const newRibbonPanels = [
      ...this.getRibbonPanels(),
      [
        ...this.runningOnInput || this.options.isTabbed ? []: [play],
        ...((this.hasUploadMode && this.isUploadMode.value) || this.runningOnInput) ? [save] : [],
        ...this.hasUploadMode ? [toggleUploadMode]: [],
      ],
    ];

    this.setRibbonPanels(newRibbonPanels);
    return newRibbonPanels;
  }

  // Main element of the output block. Stores all the tabs for the output and input
  private outputsTabsElem = ui.tabControl();

  private showOutputTabsElem() {
    $(this.outputsTabsElem.root).show();
    $(this.outputsTabsElem.root).css('display', 'flex');
  }

  public buildOutputBlock(): HTMLElement {
    this.outputsTabsElem.root.style.width = '100%';

    this.tabsLabels.forEach((tabLabel) => {
      const [tabParams, isInputTab] = this.categoryToDfParamMap.outputs[tabLabel] ? [this.categoryToDfParamMap.outputs[tabLabel], false] : [this.categoryToDfParamMap.inputs[tabLabel], true];

      const tabDfProps = tabParams.filter((p) => p.propertyType === DG.TYPE.DATA_FRAME);
      const tabScalarProps = tabParams.filter((p) => p.propertyType !== DG.TYPE.DATA_FRAME);

      const parsedTabDfProps = tabDfProps.map((dfProp) => getPropViewers(dfProp).config);

      let prevDfBlockTitle = '';
      const dfBlocks = tabDfProps.reduce((acc, dfProp, dfIndex) => {
        this.dfToViewerMapping[dfProp.name] = [];

        const promisedViewers: Promise<DG.Viewer>[] = parsedTabDfProps[dfIndex].map(async (viewerDesc: {[key: string]: string | boolean}, _) => {
          const initialValue: DG.DataFrame = this.funcCall.outputs[dfProp.name]?.value ?? this.funcCall.inputParams[dfProp.name]?.value ?? grok.data.demo.demog(1);

          const viewerType = viewerDesc['type'] as string;
          const viewer = Object.values(viewerTypesMapping).includes(viewerType) ? DG.Viewer.fromType(viewerType, initialValue): await initialValue.plot.fromType(viewerType) as DG.Viewer;
          viewer.setOptions(viewerDesc);

          this.dfToViewerMapping[dfProp.name].push(viewer);
          this.afterOutputPropertyRender.next({prop: dfProp, output: viewer});

          return viewer;
        });

        const reactiveViewers = promisedViewers.map((promisedViewer, viewerIdx) => promisedViewer.then((loadedViewer) => {
          const subscribeOnFcChanges = () => {
            const currentParam = this.funcCall.outputParams[dfProp.name] ?? this.funcCall.inputParams[dfProp.name];

            const updateViewerSource = async () => {
              const currentParam = this.funcCall.outputParams[dfProp.name] ?? this.funcCall.inputParams[dfProp.name];

              this.showOutputTabsElem();
              $(this.outputsTabsElem.getPane(tabLabel).header).show();

              if (Object.values(viewerTypesMapping).includes(loadedViewer.type)) {
                loadedViewer.dataFrame = currentParam.value;
                loadedViewer.setOptions(parsedTabDfProps[dfIndex][viewerIdx]);
              } else {
                // User-defined viewers (e.g. OutliersSelectionViewer) could created only asynchronously
                const newViewer = await currentParam.value.plot.fromType(loadedViewer.type) as DG.Viewer;
                newViewer.setOptions(parsedTabDfProps[dfIndex][viewerIdx]);
                loadedViewer.root.replaceWith(newViewer.root);
                loadedViewer = newViewer;
              }
              this.afterOutputPropertyRender.next({prop: dfProp, output: loadedViewer});
            };

            const paramSub = currentParam.onChanged.subscribe(async () => {
              await updateViewerSource();
            });

            this.funcCallReplaced.subscribe(async () => {
              await updateViewerSource();
            });

            this.subs.push(paramSub);
          };

          subscribeOnFcChanges();
          this.subs.push(
            this.funcCallReplaced.subscribe(subscribeOnFcChanges),
          );

          return loadedViewer;
        }));

        const dfBlockTitle: string = (prevDfBlockTitle !== (dfProp.options['caption'] ?? dfProp.name)) ? dfProp.options['caption'] ?? dfProp.name: '';
        prevDfBlockTitle = dfBlockTitle;

        if (isInputTab) {
          const subscribeOnFcChanges = () => {
            const currentParam = this.funcCall.outputParams[dfProp.name] ?? this.funcCall.inputParams[dfProp.name];

            const paramSub = currentParam.onChanged.subscribe(() => {
              this.showOutputTabsElem();
              Object.keys(this.categoryToDfParamMap.inputs).forEach((inputTabName) => {
                $(this.outputsTabsElem.getPane(inputTabName).header).show();
              });
            });

            this.subs.push(paramSub);
          };

          subscribeOnFcChanges();
          this.subs.push(
            this.funcCallReplaced.subscribe(subscribeOnFcChanges),
          );
        }

        const wrappedViewers = reactiveViewers.map((promisedViewer, viewerIndex) => {
          const blockWidth: string | boolean | undefined = parsedTabDfProps[dfIndex][viewerIndex]['block'];
          const viewerRoot = ui.wait(async () => (await promisedViewer).root);
          $(viewerRoot).css({
            'min-height': '300px',
            'flex-grow': '1',
          });

          return ui.divV([
            ...viewerIndex === 0 ? [ui.h2(dfBlockTitle)] : [ui.h2(' ', {style: {'white-space': 'pre'}})],
            viewerRoot,
          ], {style: {...blockWidth ? {
            'width': `${blockWidth}%`,
            'max-width': `${blockWidth}%`,
          } : {
            'flex-grow': '1',
          }}});
        });

        acc.append(...wrappedViewers);

        return acc;
      }, ui.divH([], {'style': {'flex-wrap': 'wrap', 'flex-grow': '1'}}));

      const generateScalarsTable = () => {
        const table = DG.HtmlTable.create(
          tabScalarProps,
          (scalarProp: DG.Property) => {
            const precision = scalarProp.options.precision;

            return [
              scalarProp.caption ?? scalarProp.name,
              precision && scalarProp.propertyType === DG.TYPE.FLOAT && this.funcCall.outputs[scalarProp.name]?
                this.funcCall.outputs[scalarProp.name].toPrecision(precision) : this.funcCall.outputs[scalarProp.name],
              scalarProp.options['units'],
            ];
          },
        ).root;
        $(table).css({
          'max-width': '400px',
        });
        this.afterOutputSacalarTableRender.next(table);
        return table;
      };

      let scalarsTable = generateScalarsTable();

      tabScalarProps.forEach((tabScalarProp) => {
        const subscribeOnFcChanges = () => {
          const paramSub = this.funcCall.outputParams[tabScalarProp.name].onChanged.subscribe(() => {
            const newScalarsTable = generateScalarsTable();
            scalarsTable.replaceWith(newScalarsTable);
            scalarsTable = newScalarsTable;


            $(this.outputsTabsElem.getPane(tabLabel).header).show();
          });

          this.funcCallReplaced.subscribe(() => {
            const newScalarsTable = generateScalarsTable();
            scalarsTable.replaceWith(newScalarsTable);
            scalarsTable = newScalarsTable;
          });

          this.subs.push(paramSub);
        };

        subscribeOnFcChanges();
        this.subs.push(
          this.funcCallReplaced.subscribe(subscribeOnFcChanges),
        );
      });

      this.outputsTabsElem.addPane(tabLabel, () => {
        return ui.divV([...tabDfProps.length ? [dfBlocks]: [], ...tabScalarProps.length ? [ui.h2('Scalar values'), scalarsTable]: []]);
      });
    });

    const outputBlock = ui.box();
    outputBlock.append(this.outputsTabsElem.root);

    return outputBlock;
  }

  public async onAfterLoadRun(loadedRun: DG.FuncCall) {
    this.showOutputTabsElem();
    this.outputsTabsElem.panes.forEach((tab) => {
      $(tab.header).show();
    });
  }

  // Stores mapping between DF and its' viewers
  private dfToViewerMapping: {[key: string]: DG.Viewer[]} = {};

  protected get tabsLabels() {
    return [
      ...Object.keys(this.categoryToDfParamMap.inputs),
      ...Object.keys(this.categoryToDfParamMap.outputs),
    ];
  }

  protected get categoryToDfParamMap() {
    const map = {
      inputs: {} as Record<string, DG.Property[]>,
      outputs: {} as Record<string, DG.Property[]>,
    };

    this.func.inputs
      .filter((inputProp) =>
        inputProp.propertyType === DG.TYPE.DATA_FRAME &&
        getPropViewers(inputProp).config.length !== 0,
      )
      .forEach((p) => {
        const category = p.category === 'Misc' ? 'Input': p.category;

        if (map.inputs[category])
          map.inputs[category].push(p);
        else
          map.inputs[category] = [p];
      });

    this.func.outputs
      .forEach((p) => {
        const category = p.category === 'Misc' ? 'Output': p.category;

        if (p.propertyType === DG.TYPE.DATA_FRAME &&
          getPropViewers(p).config.length === 0) return;

        if (map.outputs[category])
          map.outputs[category].push(p);
        else
          map.outputs[category] = [p];
      });

    return map;
  }

  public async doRun(): Promise<void> {
    this.isRunning.next(true);
    this.checkDisability.next();
    try {
      await this.run();
    } catch (e: any) {
      grok.shell.error(e.toString());
      console.log(e);
    } finally {
      this.isRunning.next(false);
      this.checkDisability.next();
    }
  }

  public getInput(name: string) {
    return this.inputsMap[name];
  }

  // TODO: implement warn if really needed
  public setInput(name: string, value: any, state: 'default' | 'disabled' | 'warn' = 'disabled') {
    const input = this.getInput(name);
    if (!input)
      throw new Error(`No input named ${name}`);

    // input value will not be synced, since doesn't trigger on
    // onInput for inputBase
    this.funcCall.inputs[name] = value;
    this.setInputState(input, state);
    if (state !== 'default') {
      const param: DG.FuncCallParam = this.funcCall.inputParams[name];
      param.aux.editState = state;
    }
  }

  private setInputState(input: FuncCallInput, state?: string) {
    if (state === 'disabled')
      input.enabled = false;

    if (state === 'default' || !state)
      input.enabled = true;
  }

  private renderOutputForm(): HTMLElement {
    return this.renderIOForm(SYNC_FIELD.OUTPUTS);
  }

  private renderInputForm(): HTMLElement {
    return this.renderIOForm(SYNC_FIELD.INPUTS);
  }

  private renderIOForm(field: SyncFields) {
    const inputs = ui.divH([], 'ui-form ui-form-wide');
    $(inputs).css({
      'flex-wrap': 'wrap',
      'flex-grow': '0',
    });

    let prevCategory = 'Misc';
    const params = this.funcCall[syncParams[field]].values() as DG.FuncCallParam[];
    wu(params as DG.FuncCallParam[])
      .filter((val) => !!val)
      .forEach((val) => {
        const prop = val.property;
        this.beforeInputPropertyRender.next(prop);
        const input = this.getInputForVal(val);
        if (!input) {
          prevCategory = prop.category;
          return;
        }
        this.inputsMap[val.property.name] = input;
        this.syncInput(val, input, field);
        this.disableInputsOnRun(input);
        if (field === SYNC_FIELD.INPUTS)
          this.bindOnHotkey(input);

        this.renderInput(inputs, val, input, prevCategory);
        this.afterInputPropertyRender.next({prop, input: input});
        prevCategory = prop.category;
      });

    inputs.classList.remove('ui-panel');
    inputs.style.paddingTop = '0px';
    inputs.style.paddingLeft = '0px';
    inputs.style.maxWidth = '100%';
    this.checkDisability.next();

    return inputs;
  }

  private disableInputsOnRun(t: InputVariants) {
    const disableOnRunSub = this.isRunning.subscribe((isRunning) => {
      if (isRunning)
        t.enabled = false;
      else
        t.enabled = true;
    });
    this.subs.push(disableOnRunSub);
  }

  private getInputForVal(val: DG.FuncCallParam): InputVariants | null {
    const prop = val.property;
    if (this.inputsOverride[val.property.name])
      return this.inputsOverride[val.property.name];

    if (prop.propertyType === DG.TYPE.STRING && prop.options.choices)
      return ui.choiceInput(prop.caption ?? prop.name, prop.defaultValue, JSON.parse(prop.options.choices));

    switch (prop.propertyType as any) {
    case DG.TYPE.DATA_FRAME:
      return ui.tableInput(prop.caption ?? prop.name, null, grok.shell.tables);
    case FILE_INPUT_TYPE:
      return UiUtils.fileInput(prop.caption ?? prop.name, null, null, null);
    case DG.TYPE.FLOAT:
      const floatInput = ui.input.forProperty(prop);
      const format = prop.options.format;
      if (format) floatInput.format = format;
      return floatInput;
    default:
      return ui.input.forProperty(prop);
    }
  }

  private bindOnHotkey(t: InputVariants) {
    if (isInputBase(t)) {
      t.input.onkeydown = async (ev) => {
        if (ev.key == 'Enter' && this.isRunnable()) this.doRun();
      };
    }
  }

  private renderInput(inputsDiv: HTMLDivElement, val: DG.FuncCallParam, t: InputVariants, prevCategory: string) {
    const prop = val.property;

    if (prop.category !== prevCategory)
      inputsDiv.append(ui.h2(prop.category, {style: {'width': '100%'}}));

    if (isInputBase(t))
      this.inputBaseAdditionalRenderHandler(val, t);

    inputsDiv.append(t.root);
  }

  private inputBaseAdditionalRenderHandler(val: DG.FuncCallParam, t: DG.InputBase) {
    const prop = val.property;

    $(t.root).css({
      'width': `${prop.options['block'] ?? '100'}%`,
      'box-sizing': 'border-box',
      'padding-right': '5px',
    });
    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13004
    t.captionLabel.firstChild!.replaceWith(ui.span([prop.caption ?? prop.name]));
    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13005
    if (prop.options['units']) t.addPostfix(prop.options['units']);
  }

  private syncInput(val: DG.FuncCallParam, t: InputVariants, fields: SyncFields) {
    this.syncFuncCallReplaced(t, val, fields);
    this.syncOnInput(t, val, fields);
  }

  private syncFuncCallReplaced(t: InputVariants, val: DG.FuncCallParam, field: SyncFields) {
    const prop = val.property;
    const name = val.name;
    // don't use val directly, get a fresh one, since it will be replaced with fc
    const sub = this.funcCallReplaced.pipe(startWith(true)).subscribe(() => {
      const newValue = this.funcCall[field][name] ?? prop.defaultValue ?? null;
      t.notify = false;
      t.value = newValue;
      t.notify = true;
      this.funcCall[field][name] = newValue;
      const newParam = this.funcCall[syncParams[field]][name];
      this.syncValOnChanged(t, newParam, field);
      this.setInputState(t, newParam.aux.editState);
    });
    this.subs.push(sub);
  }

  private syncValOnChanged(t: InputVariants, val: DG.FuncCallParam, field: SyncFields) {
    const sub = val.onChanged.subscribe(() => {
      const newValue = this.funcCall[field][val.name];
      // there is no notify for DG.FuncCallParam, so we need to
      // check if the value is not the same for floats, otherwise we
      // will overwrite a user input with a lower precicsion decimal
      // representation
      if (
        ((val.property.propertyType === DG.TYPE.FLOAT) && new Float32Array([t.value])[0] !== new Float32Array([newValue])[0]) ||
          val.property.propertyType !== DG.TYPE.FLOAT
      ) {
        t.notify = false;
        t.value = newValue;
        t.notify = true;
      }
      if (field === SYNC_FIELD.INPUTS) {
        this.hideOutdatedOutput();
        this.checkDisability.next();

        if (this.runningOnInput && this.isRunnable()) this.doRun();
      }
    });

    this.subs.push(sub);

    if (val.property.propertyType === DG.TYPE.DATA_FRAME) {
      const subscribeForInteriorMut = () => {
        if (!val.value) return;
        const sub = val.value.onDataChanged.subscribe(async () => {
          if (this.runningOnInput && this.isRunnable()) this.doRun();
        });
        this.subs.push(sub);
      };

      val.onChanged.subscribe(() => {
        subscribeForInteriorMut();
      });

      subscribeForInteriorMut();
    }
  }

  private syncOnInput(t: InputVariants, val: DG.FuncCallParam, field: SyncFields) {
    DG.debounce(getObservable(t.onInput.bind(t)), 350).subscribe(() => {
      if (this.isHistorical.value)
        this.isHistorical.next(false);

      this.funcCall[field][val.name] = t.value;
      if (isInputBase(t)) {
        if (t.value === null)
          setTimeout(() => t.input.classList.add('d4-invalid'), 100);
        else
          t.input.classList.remove('d4-invalid');
      }
      this.checkDisability.next();
      this.hideOutdatedOutput();
    });
  }

  private isRunnable() {
    return (wu(this.funcCall.inputs.values()).every((v) => v !== null && v !== undefined)) && !this.isRunning.value;
  }

  private async saveExperimentalRun(expFuncCall: DG.FuncCall) {
    // Dirty hack to set readonly 'started' field
    const tempCall = await(await grok.functions.eval('Sin')).prepare({x: 1}).call();
    expFuncCall.dart.r2 = tempCall.dart.r2;

    const tags = expFuncCall.options['tags'] || [];
    expFuncCall.options['tags'] = tags.includes(EXPERIMENTAL_TAG) ? tags: [...tags, EXPERIMENTAL_TAG];

    expFuncCall.newId();

    await this.saveRun(expFuncCall);
  }

  private hideOutdatedOutput() {
    this.outputsTabsElem.panes
      .filter((tab) => Object.keys(this.categoryToDfParamMap.outputs).includes(tab.name))
      .forEach((tab) => $(tab.header).hide());

    const firstInputTab = this.outputsTabsElem.panes
      .find((tab) => Object.keys(this.categoryToDfParamMap.inputs).includes(tab.name));
    if (firstInputTab)
      this.outputsTabsElem.currentPane = firstInputTab;
    else
      $(this.outputsTabsElem.root).hide();
  }

  private sheetNamesCache = {} as Record<string, string>;

  private getSheetName(initialName: string, wb: ExcelJS.Workbook) {
    if (this.sheetNamesCache[initialName]) return this.sheetNamesCache[initialName];

    let name = `${initialName}`;
    if (name.length > 31)
      name = `${name.slice(0, 31)}`;
    let i = 1;
    while (wb.worksheets.some((sheet) => sheet.name === name)) {
      let truncatedName = `${initialName}`;
      if (truncatedName.length > (31 - `-${i}`.length))
        truncatedName = `${initialName.slice(0, 31 - `-${i}`.length)}`;
      name = `${truncatedName}-${i}`;
      i++;
    }

    this.sheetNamesCache[initialName] = name;

    return name;
  };

  /**
   * RichFunctionView know everything about its UI, so it exports not only data, but also viewer screenshots.
   * This function iterates over all of the tabs and sequentally exports all dataframes, their viewers and scalars.
   * @param format format needed to export. See {@link this.defaultSupportedExportFormats} for available formats.
   * @returns Promise<Blob> with data ready for download
   */
  protected richFunctionExport = async (format: string) => {
    if (format === 'Excel') {
      try {
        const lastCall = this.lastCall;

        if (!lastCall) throw new Error(`Function was not called`);

        if (!this.exportConfig!.supportedFormats.includes(format)) throw new Error(`Format "${format}" is not supported.`);

        if (!this.func) throw new Error('The correspoding function is not specified');

        const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
        const exportWorkbook = new ExcelJS.Workbook();

        const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

        const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

        const dfInputs = this.func.inputs.filter((input) => isDataFrame(input));
        const scalarInputs = this.func.inputs.filter((input) => isScalarType(input.propertyType));
        const dfOutputs = this.func.outputs.filter((output) => isDataFrame(output));
        const scalarOutputs = this.func.outputs.filter((output) => isScalarType(output.propertyType));

        dfInputs.forEach((dfInput) => {
          const visibleTitle = dfInput.options.caption || dfInput.name;
          const currentDfSheet =
        exportWorkbook.worksheets.find((ws) => ws.name === this.getSheetName(visibleTitle, exportWorkbook)) ??
        exportWorkbook.addWorksheet(this.getSheetName(visibleTitle, exportWorkbook));

          const currentDf = lastCall.inputs[dfInput.name];
          dfToSheet(currentDfSheet, currentDf);
        });

        if (scalarInputs.length) {
          const inputScalarsSheet = exportWorkbook.addWorksheet('Input scalars');
          scalarsToSheet(inputScalarsSheet, scalarInputs.map((scalarInput) => ({
            caption: scalarInput.options['caption'] || scalarInput.name,
            value: lastCall.inputs[scalarInput.name],
            units: scalarInput.options['units'] || '',
          })));
        }

        dfOutputs.forEach((dfOutput) => {
          const visibleTitle = dfOutput.options.caption || dfOutput.name;
          const currentDfSheet =
        exportWorkbook.worksheets.find((ws) => ws.name === this.getSheetName(visibleTitle, exportWorkbook)) ??
        exportWorkbook.addWorksheet(this.getSheetName(visibleTitle, exportWorkbook));

          const currentDf = lastCall.outputs[dfOutput.name];
          dfToSheet(currentDfSheet, currentDf);
        });


        if (scalarOutputs.length) {
          const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
          scalarsToSheet(outputScalarsSheet, scalarOutputs.map((scalarOutput) => ({
            caption: scalarOutput.options['caption'] || scalarOutput.name,
            value: lastCall.outputs[scalarOutput.name],
            units: scalarOutput.options['units'] || '',
          })));
        }

        const tabControl = this.outputsTabsElem;

        for (const tabLabel of this.tabsLabels.filter((label) => Object.keys(this.categoryToDfParamMap.inputs).includes(label))) {
          for (const inputProp of this.categoryToDfParamMap.inputs[tabLabel].filter((prop) => isDataFrame(prop))) {
            const nonGridViewers = this.dfToViewerMapping[inputProp.name]
              .filter((viewer) => viewer.type !== DG.VIEWER.GRID)
              .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer.type));

            if (nonGridViewers.length === 0) continue;

            tabControl.currentPane = tabControl.getPane(tabLabel);
            await new Promise((r) => setTimeout(r, 100));

            const visibleTitle = inputProp.options.caption || inputProp.name;
            const currentDf = lastCall.inputs[inputProp.name];

            for (const [index, viewer] of nonGridViewers.entries()) {
              await plotToSheet(
                exportWorkbook,
                exportWorkbook.getWorksheet(this.getSheetName(visibleTitle, exportWorkbook)),
                viewer.root,
                currentDf.columns.length + 2,
                (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0,
              );
            };
          }
        }

        for (const tabLabel of this.tabsLabels.filter((label) => Object.keys(this.categoryToDfParamMap.outputs).includes(label))) {
          for (const outputProp of this.categoryToDfParamMap.outputs[tabLabel].filter((prop) => isDataFrame(prop))) {
            const nonGridViewers = this.dfToViewerMapping[outputProp.name]
              .filter((viewer) => viewer.type !== DG.VIEWER.GRID)
              .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer.type));

            if (nonGridViewers.length === 0) continue;

            tabControl.currentPane = tabControl.getPane(tabLabel);
            await new Promise((r) => setTimeout(r, 100));

            const visibleTitle = outputProp.options.caption || outputProp.name;
            const currentDf = lastCall.outputs[outputProp.name];

            for (const [index, viewer] of nonGridViewers.entries()) {
              if (viewer.type === DG.VIEWER.STATISTICS) {
                const length = currentDf.columns.length;
                const stats = DG.DataFrame.fromColumns([
                  DG.Column.string('Name', length).init((i: number) => currentDf.columns.byIndex(i).name),
                  DG.Column.int('Values', length).init((i: number) => currentDf.columns.byIndex(i).stats.valueCount),
                  DG.Column.int('Nulls', length).init((i: number) => currentDf.columns.byIndex(i).stats.missingValueCount),
                  DG.Column.float('Min', length).init((i: number) => currentDf.columns.byIndex(i).stats.min),
                  DG.Column.float('Max', length).init((i: number) => currentDf.columns.byIndex(i).stats.max),
                  DG.Column.float('Avg', length).init((i: number) => currentDf.columns.byIndex(i).stats.avg),
                  DG.Column.float('Stdev', length).init((i: number) => currentDf.columns.byIndex(i).stats.stdev),
                ]);
                dfToSheet(
                  exportWorkbook.getWorksheet(this.getSheetName(visibleTitle, exportWorkbook)),
                  stats,
                  currentDf.columns.length + 2,
                  (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0,
                );
              } else {
                await plotToSheet(
                  exportWorkbook,
                  exportWorkbook.getWorksheet(this.getSheetName(visibleTitle, exportWorkbook)),
                  viewer.root,
                  currentDf.columns.length + 2,
                  (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0,
                );
              }
            }
          }
        }

        const buffer = await exportWorkbook.xlsx.writeBuffer();

        return new Blob([buffer], {type: BLOB_TYPE});
      } catch (e) {
        console.log(e);
      }
    }

    if (format === 'DataUrl images') {
      const jsonText = {} as Record<string, Record<number, {dataUrl: string, width: number, height: number}>>;

      const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

      const tabControl = this.outputsTabsElem;

      for (const tabLabel of this.tabsLabels.filter((label) => Object.keys(this.categoryToDfParamMap.inputs).includes(label))) {
        for (const inputProp of this.categoryToDfParamMap.inputs[tabLabel].filter((prop) => isDataFrame(prop))) {
          const nonGridViewers = this.dfToViewerMapping[inputProp.name]
            .filter((viewer) => viewer.type !== DG.VIEWER.GRID && viewer.type !== DG.VIEWER.STATISTICS)
            .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer.type));

          if (nonGridViewers.length === 0) continue;

          tabControl.currentPane = tabControl.getPane(tabLabel);
          await new Promise((r) => setTimeout(r, 100));

          for (const [i, viewer] of nonGridViewers.entries()) {
            const dataUrl = (await html2canvas(viewer.root, {logging: false})).toDataURL();

            if (!jsonText[inputProp.name]) jsonText[inputProp.name] = {};

            jsonText[inputProp.name][i] = {dataUrl, width: viewer.root.clientWidth, height: viewer.root.clientHeight};
          }
        }
      }

      for (const tabLabel of this.tabsLabels.filter((label) => Object.keys(this.categoryToDfParamMap.outputs).includes(label))) {
        for (const outputProp of this.categoryToDfParamMap.outputs[tabLabel].filter((prop) => isDataFrame(prop))) {
          const nonGridViewers = this.dfToViewerMapping[outputProp.name]
            .filter((viewer) => viewer.type !== DG.VIEWER.GRID && viewer.type !== DG.VIEWER.STATISTICS)
            .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer.type));

          if (nonGridViewers.length === 0) continue;

          tabControl.currentPane = tabControl.getPane(tabLabel);
          await new Promise((r) => setTimeout(r, 100));

          for (const [i, viewer] of nonGridViewers.entries()) {
            const dataUrl = (await html2canvas(viewer.root, {logging: false})).toDataURL();

            if (!jsonText[outputProp.name]) jsonText[outputProp.name] = {};

            jsonText[outputProp.name][i] = {dataUrl, width: viewer.root.clientWidth, height: viewer.root.clientHeight};
          }
        };
      }
      return new Blob([JSON.stringify(jsonText)], {type: 'text/plain'});
    }

    throw new Error('Format is not supported');
  };

  richFunctionViewSupportedFormats() {
    return ['Excel', 'DataUrl images'];
  }

  richFunctionViewExportExtensions() {
    return {
      'Excel': 'xlsx',
      'DataUrl images': 'txt',
    };
  }

  exportConfig = {
    supportedExtensions: this.richFunctionViewExportExtensions(),
    supportedFormats: this.richFunctionViewSupportedFormats(),
    export: this.richFunctionExport,
    filename: this.defaultExportFilename,
  };
}

const scalarsToSheet = (sheet: ExcelJS.Worksheet, scalars: { caption: string, value: string, units: string }[]) => {
  sheet.addRow(['Parameter', 'Value', 'Units']).font = {bold: true};
  scalars.forEach((scalar) => {
    sheet.addRow([scalar.caption, scalar.value, scalar.units]);
  });

  sheet.getColumn(1).width = Math.max(
    ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length,
  ) * 1.2;
  sheet.getColumn(2).width = Math.max(...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
  sheet.getColumn(3).width = Math.max(...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;
};

let dfCounter = 0;
const dfToSheet = (sheet: ExcelJS.Worksheet, df: DG.DataFrame, column?: number, row?: number) => {
  const columnKey = sheet.getColumn(column ?? 1).letter;
  const tableConfig = {
    name: `ID_${dfCounter.toString()}`,
    ref: `${columnKey}${row ?? 1}`,
    columns: df.columns.toList().map((col) => ({name: col.name, filterButton: false})),
    rows: new Array(df.rowCount).fill(0).map((_, idx) => [...df.row(idx).cells].map((cell) => cell.value)),
  };
  sheet.addTable(tableConfig);
  sheet.columns.forEach((col) => {
    col.width = 25;
    col.alignment = {wrapText: true};
  });
  dfCounter++;
};

const plotToSheet = async (exportWb: ExcelJS.Workbook, sheet: ExcelJS.Worksheet, plot: HTMLElement, columnForImage: number, rowForImage: number = 0) => {
  const canvas = await html2canvas(plot as HTMLElement, {logging: false});
  const dataUrl = canvas.toDataURL('image/png');

  const imageId = exportWb.addImage({
    base64: dataUrl,
    extension: 'png',
  });
  sheet.addImage(imageId, {
    tl: {col: columnForImage, row: rowForImage},
    ext: {width: canvas.width, height: canvas.height},
  });
};
