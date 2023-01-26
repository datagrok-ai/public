/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView, INTERACTIVE_CSS_CLASS} from './function-view';
import ExcelJS from 'exceljs';
import html2canvas from 'html2canvas';
import wu from 'wu';
import $ from 'cash-dom';
import {Subject} from 'rxjs';
import {UiUtils} from './shared-components/ui-utils';

const viewerTypesMapping: {[key: string]: string} = {
  ['barchart']: DG.VIEWER.BAR_CHART,
  ['boxplot']: DG.VIEWER.BOX_PLOT,
  ['calendar']: DG.VIEWER.CALENDAR,
  ['corrplot']: DG.VIEWER.CORR_PLOT,
  ['densityplot']: DG.VIEWER.DENSITY_PLOT,
  ['filters']: DG.VIEWER.FILTERS,
  ['form']: DG.VIEWER.FORM,
  ['globe']: DG.VIEWER.GLOBE,
  ['googlemap']: DG.VIEWER.GOOGLE_MAP,
  ['grid']: DG.VIEWER.GRID,
  ['heatmap']: DG.VIEWER.HEAT_MAP,
  ['histogram']: DG.VIEWER.HISTOGRAM,
  ['linechart']: DG.VIEWER.LINE_CHART,
  ['markup']: DG.VIEWER.MARKUP,
  ['matrixplot']: DG.VIEWER.MATRIX_PLOT,
  ['networkdiagram']: DG.VIEWER.NETWORK_DIAGRAM,
  ['pcplot']: DG.VIEWER.PC_PLOT,
  ['piechart']: DG.VIEWER.PIE_CHART,
  ['scatterplot']: DG.VIEWER.SCATTER_PLOT,
  ['scatterplot3d']: DG.VIEWER.SCATTER_PLOT_3D,
  ['shapemap']: DG.VIEWER.SHAPE_MAP,
  ['statistics']: DG.VIEWER.STATISTICS,
  ['tileviewer']: DG.VIEWER.TILE_VIEWER,
  ['timelines']: DG.VIEWER.TIMELINES,
  ['treemap']: DG.VIEWER.TREE_MAP,
  ['trellisplot']: DG.VIEWER.TRELLIS_PLOT,
  ['wordcloud']: DG.VIEWER.WORD_CLOUD,
};

const FILE_INPUT_TYPE = 'file';

export class RichFunctionView extends FunctionView {
  // emitted when after a new FuncCall is linked
  private funcCallReplaced = new Subject<true>();

  // emitted when runButton disability should be checked
  private checkDisability = new Subject();

  // stores the running state
  private isRunning = false;

  constructor(funcCall: DG.FuncCall) {
    super();

    this.basePath = `scripts/${funcCall.func.id}/view`;

    this.exportConfig = {
      supportedExtensions: this.defaultSupportedExportExtensions(),
      supportedFormats: this.defaultSupportedExportFormats(),
      export: this.defaultExport,
      filename: this.defaultExportFilename,
    };

    this.linkFunccall(funcCall);
    this.init();
    this.build();

    this.funcCallReplaced.next();

    this.name = funcCall.func.friendlyName;
  }

  public override linkFunccall(funcCall: DG.FuncCall) {
    const isPreviousHistorical = this._funcCall?.options['isHistorical'];
    this._funcCall = funcCall;

    this.funcCallReplaced.next(true);

    if (funcCall.options['isHistorical']) {
      if (!isPreviousHistorical)
        this.name = `${this.name} — ${funcCall.options['title'] ?? new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`;
      else
        this.name = `${this.name.substring(0, this.name.indexOf(' — '))} — ${funcCall.options['title'] ?? new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`;


      // FIX ME: view name does not change in models
      document.querySelector('div.d4-ribbon-name')?.replaceChildren(ui.span([this.name]));
      this.path = `?id=${this._funcCall.id}`;
    } else {
      this.path = ``;

      if (isPreviousHistorical)
        this.name = `${this.name.substring(0, this.name.indexOf(' — '))}`;
    }
    this.buildRibbonPanels();
  }

  public override onAfterRun(runFunc: DG.FuncCall): Promise<void> {
    this.tabsElem.root.style.removeProperty('display');
    this.tabsElem.panes.forEach((tab) => {
      tab.header.style.removeProperty('display');
    });

    return Promise.resolve();
  }

  public buildIO(): HTMLElement {
    const inputBlock = this.buildInputBlock();

    ui.tools.handleResize(inputBlock, (width) => {
      if (width < 350)
        $(inputBlock.firstChild!).addClass('ui-form-condensed');
      else
        $(inputBlock.firstChild!).removeClass('ui-form-condensed');
    });

    const outputBlock = this.buildOutputBlock();
    outputBlock.style.height = '100%';
    outputBlock.style.width = '100%';
    this.tabsElem.root.style.display = 'none';

    if (!!this.tabsElem.getPane('Input')) {
      this.tabsElem.panes.forEach((tab) => {
        tab.header.style.display = 'none';
      });
    }

    const out = ui.splitH([inputBlock, ui.panel([outputBlock])], null, true);
    out.style.padding = '0 12px';

    inputBlock.parentElement!.style.maxWidth = '450px';

    return out;
  }

  public buildInputBlock(): HTMLElement {
    const formDiv = this.renderRunSection();

    return formDiv;
  }

  private tabsElem = ui.tabControl();

  public buildOutputBlock(): HTMLElement {
    const outputBlock = ui.box();

    this.tabsLabels.forEach((tabLabel) => {
      const tabDfProps = this.categoryToParamMap[tabLabel].filter((p) => p.propertyType === DG.TYPE.DATA_FRAME);
      const tabScalarProps = this.categoryToParamMap[tabLabel].filter((p) => p.propertyType !== DG.TYPE.DATA_FRAME);
      const parsedTabDfProps = tabDfProps.map((dfProp) => JSON.parse(dfProp.options['viewer']));

      const dfSections = tabDfProps.map((dfProp, dfIndex) => {
        const viewers: DG.Viewer[] = parsedTabDfProps[dfIndex].map((viewerDesc: {[key: string]: string}) => {
          const initialValue = this.funcCall?.outputs[dfProp.name]?.value ?? this.funcCall!.inputParams[dfProp.name]?.value ?? DG.DataFrame.fromCsv(``);

          const viewer = DG.Viewer.fromType(viewerTypesMapping[viewerDesc['type']], initialValue);
          viewer.setOptions(viewerDesc);

          if (this.dfToViewerMapping[dfProp.name]) this.dfToViewerMapping[dfProp.name].push(viewer); else this.dfToViewerMapping[dfProp.name] = [viewer];

          return viewer;
        });

        this.funcCallReplaced.subscribe(() => {
          const currentParam: DG.FuncCallParam = this.funcCall!.outputParams[dfProp.name] ?? this.funcCall!.inputParams[dfProp.name];

          currentParam.onChanged.subscribe(() => {
            viewers.forEach((viewer) => {
              viewer.dataFrame = currentParam.value;
            });
          });
        });

        const dfSectionTitle: string = dfProp.options['caption'] ?? dfProp.name;

        if (tabLabel !== 'Input') {
          return ui.divV([
            ui.h2(dfSectionTitle),
            ui.divH(viewers.map((viewer, viewerIndex) => {
              if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '25') return ui.block25([viewer.root]);
              if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '50') return ui.block50([viewer.root]);
              if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '75') return ui.block75([viewer.root]);
              if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '100') return ui.block([viewer.root]);

              return viewer.root;
            }), {style: {'flex-wrap': 'wrap'}})
          ]);
        } else {
          const inputSection = ui.divH(viewers.map((viewer, viewerIndex) => {
            if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '25') return ui.block25([viewer.root]);
            if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '50') return ui.block50([viewer.root]);
            if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '75') return ui.block75([viewer.root]);
            if (parsedTabDfProps[dfIndex][viewerIndex]['block'] === '100') return ui.block([viewer.root]);

            return viewer.root;
          }), {style: {'flex-wrap': 'wrap'}});

          this.funcCallReplaced.subscribe(() => {
            const currentParam: DG.FuncCallParam = this.funcCall!.outputParams[dfProp.name] ?? this.funcCall!.inputParams[dfProp.name];

            currentParam.onChanged.subscribe(() => {
              this.tabsElem.root.style.removeProperty('display');
            });
          });

          return ui.divV([
            ui.h2(dfSectionTitle),
            inputSection
          ]);
        }
      });


      const generateScalarsTable = () => ui.tableFromMap(
        tabScalarProps.reduce((acc, scalar) => {
          acc[scalar.name] = this.funcCall?.outputs[scalar.name];
          return acc;
        }, {} as {[key: string]: string})
      );

      let scalarsTable = generateScalarsTable();

      tabScalarProps.forEach((tabScalarProp) => {
        this.funcCallReplaced.subscribe(() => {
          (this.funcCall!.outputParams[tabScalarProp.name] as DG.FuncCallParam).onChanged.subscribe(() => {
            const newScalarsTable = generateScalarsTable();
            scalarsTable.replaceWith(newScalarsTable);
            scalarsTable = newScalarsTable;
          });
        });
      });

      this.tabsElem.addPane(tabLabel, () => {
        return ui.divV([...dfSections, ...tabScalarProps.length ? [ui.h2('Scalar values')]: [], scalarsTable]);
      });
    });
    outputBlock.append(this.tabsElem.root);

    this.tabsElem.header.classList.add(INTERACTIVE_CSS_CLASS);

    return outputBlock;
  }

  public async onAfterLoadRun(loadedRun: DG.FuncCall) {
    wu(this.funcCall!.outputParams.values() as DG.FuncCallParam[]).forEach((out) => {
      this.funcCall!.setParamValue(out.name, loadedRun.outputs[out.name]);
    });

    wu(this.funcCall!.inputParams.values() as DG.FuncCallParam[]).forEach((inp) => {
      this.funcCall!.setParamValue(inp.name, loadedRun.inputs[inp.name]);
    });

    this.tabsElem.root.style.removeProperty('display');
    this.tabsElem.panes.forEach((tab) => {
      tab.header.style.removeProperty('display');
    });
  }

  private dfToViewerMapping: {[key: string]: DG.Viewer[]} = {};

  protected get dfParams() {
    return [
      ...this.func!.inputs.filter((p) => p.propertyType == DG.TYPE.DATA_FRAME),
      ...this.func!.outputs.filter((p) => p.propertyType == DG.TYPE.DATA_FRAME),
    ];
  }

  protected get dfOutputParams() {
    return this.func!.outputs.filter((p) => p.propertyType == DG.TYPE.DATA_FRAME);
  }

  protected get isInputPanelRequired() {
    return this.func?.inputs.some((p) => p.propertyType == DG.TYPE.DATA_FRAME && p.options['viewer'] != null) || false;
  }

  protected get outUniqueParamCategories() {
    return [
      ...new Set(this.func!.outputs.map((p) => p.category)) // get all output params' categories
    ]; // keep only unique of them
  }

  protected get outputTabsLabels() {
    return [
      ...this.outUniqueParamCategories,
      ...this.outUniqueParamCategories.find((val) => val === 'Misc') ? ['Output'] : [], // if no categories are stated, the default category is added
    ];
  }

  protected get tabsLabels() {
    return Object.keys(this.categoryToParamMap);
  }

  protected get categoryToParamMap() {
    const map = {} as Record<string, DG.Property[]>;

    if (this.isInputPanelRequired) {
      this.func!.inputs.
        filter((p) => p.propertyType == DG.TYPE.DATA_FRAME && p.options['viewer'] != null).
        forEach((p) => map['Input'] ? map['Input'].push(p): map['Input'] = [p]);
    }

    this.func!.outputs.forEach((p) => {
      const category = p.category === 'Misc' ? 'Output': p.category;

      if (map[category])
        map[category].push(p);
      else
        map[category] = [p];
    });

    return map;
  }

  private renderRunSection(): HTMLElement {
    const runButton = ui.bigButton('Run', async () => {
      this.isRunning = true;
      this.checkDisability.next();
      try {
        await this.run();
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        this.isRunning = false;
        this.checkDisability.next();
      }
    });
    const buttonWrapper = ui.div(runButton);
    ui.tooltip.bind(buttonWrapper, () => runButton.disabled ? (this.isRunning ? 'Computations are in progress' : 'Some inputs are invalid') : '');

    this.checkDisability.subscribe(() => {
      const isValid = (wu(this.funcCall!.inputs.values()).every((v) => v !== null)) && !this.isRunning;
      runButton.disabled = isValid ? false : true;
    });

    const inputs = ui.divV([], 'ui-form');
    let prevCategory = 'Misc';
    wu(this.funcCall!.inputParams.values() as DG.FuncCallParam[])
      .filter((val) => !!val)
      .forEach((val) => {
        const prop = val.property;
        if (prop.propertyType.toString() === FILE_INPUT_TYPE) {
          const t = UiUtils.fileInput(prop.caption ?? prop.name, null, (file: File) => this.funcCall!.inputs[prop.name] = file);
          if (prop.category !== prevCategory)
            inputs.append(ui.h2(prop.category));

          inputs.append(t.root);
        } else {
          // FIX ME: .toDart added to prevent bug of tables initial non-rendering
          const t = prop.propertyType === DG.TYPE.DATA_FRAME ?
            ui.tableInput(prop.name, null, grok.shell.tables.map(DG.toDart)):
            ui.input.forProperty(prop);

          t.captionLabel.firstChild!.replaceWith(ui.span([prop.caption ?? prop.name]));
          if (prop.options['units']) t.addPostfix(prop.options['units']);

          this.funcCallReplaced.subscribe(() => {
            t.value = this.funcCall!.inputs[val.name] ?? prop.defaultValue ?? null;
          });
          t.onChanged(() => {
            this.funcCall!.inputs[val.name] = t.value;
            if (t.value === null) setTimeout(() => t.input.classList.add('d4-invalid'), 100); else t.input.classList.remove('d4-invalid'); ;

            this.checkDisability.next();
          });
          if (prop.category !== prevCategory)
            inputs.append(ui.h2(prop.category));

          if (t.value === null) runButton.disabled = true;
          inputs.append(t.root);
        }
        prevCategory = prop.category;
      });
    // @ts-ignore
    const buttons = ui.buttonsInput([buttonWrapper]);
    inputs.append(buttons);
    inputs.classList.remove('ui-panel');
    inputs.style.paddingTop = '0px';
    inputs.style.paddingLeft = '0px';

    return ui.div([inputs]);
  }

  protected defaultExport = async (format: string) => {
    const lastCall = this.lastCall;
    if (!lastCall) throw new Error(`Function was not called`);

    if (!this.exportConfig!.supportedFormats.includes(format)) throw new Error(`Format "${format}" is not supported.`);

    if (!this.func) throw new Error('The correspoding function is not specified');

    const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
    const exportWorkbook = new ExcelJS.Workbook();

    const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

    const isDataFrame = (type: DG.TYPE) => (type === DG.TYPE.DATA_FRAME);

    const dfInputs = this.func.inputs.filter((input) => isDataFrame(input.propertyType));
    const scalarInputs = this.func.inputs.filter((input) => isScalarType(input.propertyType));
    const dfOutputs = this.func.outputs.filter((output) => isDataFrame(output.propertyType));
    const scalarOutputs = this.func.outputs.filter((output) => isScalarType(output.propertyType));

    const inputParams = [...lastCall.inputParams.values()] as DG.FuncCallParam[];
    const outputParams = [...lastCall.outputParams.values()] as DG.FuncCallParam[];

    dfInputs.forEach((dfInput) => {
      const visibleTitle = dfInput.options.caption || dfInput.name;
      const currentDfSheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT));

      const currentDf = (lastCall.inputs[dfInput.name] as DG.DataFrame);
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
      const currentDfSheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT));

      const currentDf = (lastCall.outputs[dfOutput.name] as DG.DataFrame);
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

    const tabControl = this.tabsElem;
    for (const tabLabel of this.tabsLabels) {
      tabControl.currentPane = tabControl.getPane(tabLabel);
      await new Promise((r) => setTimeout(r, 100));
      if (tabLabel === 'Input') {
        for (const inputParam of inputParams.filter((inputParam) => inputParam.property.propertyType === DG.TYPE.DATA_FRAME)) {
          const nonGridViewers = this.dfToViewerMapping[inputParam.name].filter((viewer) => viewer.type !== DG.VIEWER.GRID);

          const dfInput = dfInputs.find((input) => input.name === inputParam.name);
          const visibleTitle = dfInput!.options.caption || inputParam.name;
          const currentDf = (lastCall.inputs[dfInput!.name] as DG.DataFrame);

          for (const [index, viewer] of nonGridViewers.entries()) {
            await plotToSheet(
              exportWorkbook,
              exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT)),
              viewer.root,
              currentDf.columns.length + 2,
              (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0
            );
          };
        }
      } else {
        for (const outputParam of outputParams.filter((outputParam) => outputParam.property.propertyType === DG.TYPE.DATA_FRAME && outputParam.property.category === tabLabel)) {
          const nonGridViewers = this.dfToViewerMapping[outputParam.name].filter((viewer) => viewer.type !== DG.VIEWER.GRID);

          const dfOutput = dfOutputs.find((input) => input.name === outputParam.name)!;
          const visibleTitle = dfOutput.options.caption || outputParam.name;
          const currentDf = (lastCall.outputs[dfOutput.name] as DG.DataFrame);

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
                exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT)),
                stats,
                currentDf.columns.length + 2,
                (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0
              );
            } else {
              await plotToSheet(
                exportWorkbook,
                exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT)),
                viewer.root,
                currentDf.columns.length + 2,
                (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0
              );
            }
          }
        };
      }
    };
    const buffer = await exportWorkbook.xlsx.writeBuffer();

    return new Blob([buffer], {type: BLOB_TYPE});
  };
}

const getSheetName = (name: string, direction: DIRECTION) => {
  const idealName = `${direction} - ${name}`;
  return (idealName.length > 31) ? name.substring(0, 32) : idealName;
};

enum DIRECTION {
  INPUT = 'Input',
  OUTPUT = 'Output'
}

const scalarsToSheet = (sheet: ExcelJS.Worksheet, scalars: { caption: string, value: string, units: string }[]) => {
  sheet.addRow(['Parameter', 'Value', 'Units']).font = {bold: true};
  scalars.forEach((scalar) => {
    sheet.addRow([scalar.caption, scalar.value, scalar.units]);
  });

  sheet.getColumn(1).width = Math.max(
    ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length
  ) * 1.2;
  sheet.getColumn(2).width = Math.max(...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
  sheet.getColumn(3).width = Math.max(...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;
};

const dfToSheet = (sheet: ExcelJS.Worksheet, df: DG.DataFrame, column: number = 0, row: number = 0) => {
  for (let i= 0; i < df.columns.names().length; i++) {
    sheet.getCell(1 + row, 1 + i + column).value = df.columns.byIndex(i).name;
    sheet.getColumn(1 + i + column).width = Math.max(
      ...df.columns.byIndex(i).categories.map((category) => category.toString().length),
      df.columns.byIndex(i).name.length
    ) * 1.2;
  }
  for (let dfColumn = 0; dfColumn < df.columns.length; dfColumn++) {
    for (let i = 0; i < df.rowCount; i++)
      sheet.getCell(i + 2 + row, 1 + column+dfColumn).value = df.columns.byIndex(dfColumn).get(i);
  }
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
