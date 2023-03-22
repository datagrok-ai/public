/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from './function-view';
import ExcelJS from 'exceljs';
import html2canvas from 'html2canvas';
import wu from 'wu';
import $ from 'cash-dom';
import {Subject} from 'rxjs';
import {UiUtils} from './shared-components/ui-utils';
import '../css/rich-function-view.css';

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
  // emitted when runButton disability should be checked
  private checkDisability = new Subject();

  // stores the running state
  private isRunning = false;

  static fromFuncCall(
    funcCall: DG.FuncCall,
    options: {historyEnabled: boolean, isTabbed: boolean} =
    {historyEnabled: true, isTabbed: false}
  ) {
    return new this(funcCall.func.nqName, options);
  }

  constructor(
    funcName: string,
    public options: { historyEnabled: boolean, isTabbed: boolean} =
    {historyEnabled: true, isTabbed: false}
  ) {
    super(funcName, options);

    this.basePath = `scripts/${this.funcCall.func.id}/view`;
  }

  public override onAfterRun(runFunc: DG.FuncCall): Promise<void> {
    this.tabsElem.root.style.removeProperty('display');
    this.tabsElem.panes.forEach((tab) => {
      tab.header.style.removeProperty('display');
    });
    const firstOutputTab = this.tabsElem.panes.find((tab) => tab.name !== 'Input');
    if (firstOutputTab) this.tabsElem.currentPane = firstOutputTab;

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

    const out = ui.splitH([inputBlock, ui.panel([outputBlock], {style: {'padding-top': '0px'}})], null, true);
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
    this.tabsElem.root.style.width = '100%';

    this.tabsLabels.forEach((tabLabel) => {
      const tabDfProps = this.categoryToParamMap[tabLabel].filter((p) => p.propertyType === DG.TYPE.DATA_FRAME);
      const tabScalarProps = this.categoryToParamMap[tabLabel].filter((p) => p.propertyType !== DG.TYPE.DATA_FRAME);

      // EXPLAIN: undefined and JSON.parse
      const parsedTabDfProps = tabDfProps.map((dfProp) => (dfProp.options['viewer'] !== undefined) ?
        JSON.parse(dfProp.options['viewer'], (_, v) => v === 'true' || v === 'false' ? JSON.parse(v): v):
        []);

      const dfBlocks = tabDfProps.map((dfProp, dfIndex) => {
        let fullBlockWidth = 0;
        parsedTabDfProps[dfIndex].forEach((_: any, viewerIndex: number) => {
          const blockSize = parsedTabDfProps[dfIndex][viewerIndex]['block'];

          fullBlockWidth += blockSize ? parseInt(blockSize): 0;
        });

        const viewers: Promise<DG.Viewer[]> = Promise.all(parsedTabDfProps[dfIndex].map(async (viewerDesc: {[key: string]: string}) => {
          const initialValue: DG.DataFrame = this.funcCall?.outputs[dfProp.name]?.value ?? this.funcCall!.inputParams[dfProp.name]?.value ?? grok.data.demo.demog(1);

          const viewer = viewerTypesMapping[viewerDesc['type']] ? DG.Viewer.fromType(viewerTypesMapping[viewerDesc['type']], initialValue): await initialValue.plot.fromType(viewerDesc['type']) as DG.Viewer;
          viewer.setOptions(viewerDesc);

          if (this.dfToViewerMapping[dfProp.name]) this.dfToViewerMapping[dfProp.name].push(viewer); else this.dfToViewerMapping[dfProp.name] = [viewer];

          return viewer;
        }));

        viewers.then((loadedViewers) => {
          const subscribeOnFcChanges = () => {
            const currentParam: DG.FuncCallParam = this.funcCall!.outputParams[dfProp.name] ?? this.funcCall!.inputParams[dfProp.name];

            currentParam.onChanged.subscribe(() => {
              loadedViewers.forEach(async (viewer) => {
                if (Object.values(viewerTypesMapping).includes(viewer.type)) {
                  viewer.dataFrame = currentParam.value;
                } else {
                  // EXPLAIN WHY
                  const newViewer = await currentParam.value.plot.fromType(viewer.type) as DG.Viewer;
                  viewer.root.replaceWith(newViewer.root);
                  viewer = newViewer;
                }
              });
            });
          };

          subscribeOnFcChanges();
          this.funcCallReplaced.subscribe(subscribeOnFcChanges);
        });

        const dfBlockTitle: string = dfProp.options['caption'] ?? dfProp.name;

        const dfBlock = ui.wait(async () => {
          const loadedViewers = await viewers;
          return ui.divH(
            await Promise.all(
              loadedViewers.map((viewer) => viewer.root, {style: {'flex-wrap': 'wrap'}})
            )
          );
        });

        if (tabLabel === 'Input') {
          const subscribeOnFcChanges = () => {
            const currentParam: DG.FuncCallParam = this.funcCall!.outputParams[dfProp.name] ?? this.funcCall!.inputParams[dfProp.name];

            currentParam.onChanged.subscribe(() => {
              this.tabsElem.root.style.removeProperty('display');
            });
          };

          subscribeOnFcChanges();
          this.funcCallReplaced.subscribe(subscribeOnFcChanges);
        }

        return ui.divV([
          // Removing title margin on top of block titles
          ui.h2(dfBlockTitle, {...(dfIndex > 0) ? {style: {'margin-top': '0px'}} : {}}),
          dfBlock
        ], {style: {...(fullBlockWidth > 0) ? {'width': `${fullBlockWidth}%`}: {'width': `100%`}}});
      });


      const generateScalarsTable = () => ui.tableFromMap(
        tabScalarProps.reduce((acc, scalar) => {
          acc[scalar.name] = this.funcCall?.outputs[scalar.name];
          return acc;
        }, {} as {[key: string]: string})
      );

      let scalarsTable = generateScalarsTable();

      tabScalarProps.forEach((tabScalarProp) => {
        const subscribeOnFcChanges = () => {
          (this.funcCall!.outputParams[tabScalarProp.name] as DG.FuncCallParam).onChanged.subscribe(() => {
            const newScalarsTable = generateScalarsTable();
            scalarsTable.replaceWith(newScalarsTable);
            scalarsTable = newScalarsTable;
          });
        };

        subscribeOnFcChanges();
        this.funcCallReplaced.subscribe(subscribeOnFcChanges);
      });

      const dfSections = [] as HTMLElement[];
      let currentWidth = 0;
      let currentSection = ui.divH([], {style: {'margin-bottom': '25px'}}); // Adding small margin between sections
      dfBlocks.forEach((dfBlock, dfIndex) => {
        const blockWidth = parseInt(dfBlock.style.width);

        if (currentWidth + blockWidth === 100) {
          currentSection.append(dfBlock);
          dfSections.push(currentSection);
          currentSection = ui.divH([]);
          currentWidth = 0;
        } else if (currentWidth + blockWidth < 100) {
          currentSection.append(dfBlock);
          currentWidth += blockWidth;
        } else {
          dfSections.push(currentSection);
          currentSection = ui.divH([]);
          currentSection.append(dfBlock);
          currentWidth = blockWidth;
        }

        if (dfIndex === dfBlocks.length - 1)
          dfSections.push(currentSection);
      });

      this.tabsElem.addPane(tabLabel, () => {
        return ui.divV([...dfSections, ...tabScalarProps.length ? [ui.h2('Scalar values')]: [], scalarsTable]);
      });
    });

    const outputBlock = ui.box();
    outputBlock.append(this.tabsElem.root);

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
      ...this.outUniqueParamCategories.find((val) => val === 'Misc' || val === 'Output') ? ['Output'] : [], // if no categories are stated, the default category is added
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

  public override async run(): Promise<void> {
    if (!this.funcCall) throw new Error('The correspoding function is not specified');

    await this.onBeforeRun(this.funcCall);
    const pi = DG.TaskBarProgressIndicator.create('Calculating...');
    this.funcCall.newId();
    await this.funcCall.call(); // mutates the funcCall field
    pi.close();
    await this.onAfterRun(this.funcCall);
    this.lastCall = this.options.isTabbed ? this.funcCall.clone() : await this.saveRun(this.funcCall);
  }

  private async doRun(): Promise<void> {
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
  }

  private renderRunSection(): HTMLElement {
    const runButton = ui.bigButton('Run', async () => await this.doRun());
    // REPLACE BY TYPE GUARD
    const isFuncScript = () => {
      //@ts-ignore
      return !!this.func.script;
    };
    // TO DO: move button somewhere
    const openScriptBtn = ui.button('Open script', async () => {
      window.open(`${window.location.origin}/script/${(this.func as DG.Script).id}`, '_blank');
    });
    // const buttonWrapper = ui.div([...this.options.isTabbed && isFuncScript() ? [openScriptBtn]: [], runButton]);
    const buttonWrapper = ui.div([runButton]);
    ui.tooltip.bind(buttonWrapper, () => runButton.disabled ? (this.isRunning ? 'Computations are in progress' : 'Some inputs are invalid') : '');

    this.checkDisability.subscribe(() => {
      const isValid = (wu(this.funcCall!.inputs.values()).every((v) => v !== null && v !== undefined)) && !this.isRunning;
      runButton.disabled = !isValid;
    });

    const inputs = ui.divV([], 'ui-form');
    let prevCategory = 'Misc';
    wu(this.funcCall!.inputParams.values() as DG.FuncCallParam[])
      .filter((val) => !!val)
      .forEach((val) => {
        const prop = val.property;
        if (prop.propertyType.toString() === FILE_INPUT_TYPE) {
          const t = UiUtils.fileInput(prop.caption ?? prop.name, null, (file: File) => {
            this.funcCall!.inputs[prop.name] = file;
            this.checkDisability.next();
          });
          if (prop.category !== prevCategory)
            inputs.append(ui.h2(prop.category));

          this.checkDisability.next();
          inputs.append(t.root);
        } else {
          let t = prop.propertyType === DG.TYPE.DATA_FRAME ?
            ui.tableInput(prop.caption ?? prop.name, null, grok.shell.tables):
            ui.input.forProperty(prop);

          t.input.onkeydown = async (ev) => {
            if (ev.key == 'Enter')
              await this.doRun();
          };

          t.captionLabel.firstChild!.replaceWith(ui.span([prop.caption ?? prop.name]));
          if (prop.options['units']) t.addPostfix(prop.options['units']);

          const subscribeOnFcChanges = () => {
            t.value = this.funcCall!.inputs[val.name] ?? prop.defaultValue ?? null;
            this.funcCall!.inputs[val.name] = this.funcCall!.inputs[val.name] ?? prop.defaultValue ?? null;
          };

          subscribeOnFcChanges();
          this.funcCallReplaced.subscribe(subscribeOnFcChanges);

          t.onChanged(() => {
            // Resolving case of propertyType = DATAFRAME
            if (!!t.property) return;

            this.funcCall!.inputs[val.name] = t.value;
          });

          t.onInput(() => {
            this.funcCall!.inputs[val.name] = t.value;
            if (t.value === null) setTimeout(() => t.input.classList.add('d4-invalid'), 100); else t.input.classList.remove('d4-invalid'); ;

            this.checkDisability.next();
          });

          val.onChanged.subscribe(() => {
            const newValue = this.funcCall!.inputs[val.name];

            if (val.property.propertyType === DG.TYPE.DATA_FRAME) {
              // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12223
              const newTableInput = ui.tableInput(prop.caption ?? prop.name, newValue, [...grok.shell.tables, newValue]);
              t.root.replaceWith(newTableInput.root);
              t = newTableInput;

              t.onInput(() => {
                this.funcCall!.inputs[val.name] = t.value;
                if (t.value === null) setTimeout(() => t.input.classList.add('d4-invalid'), 100); else t.input.classList.remove('d4-invalid'); ;

                this.checkDisability.next();
              });
            } else {
              t.notify = false;
              t.value = newValue;
              t.notify = true;
            }

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
          const nonGridViewers = this.dfToViewerMapping[inputParam.name]
            .filter((viewer) => viewer.type !== DG.VIEWER.GRID)
            .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer.type));

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
        for (const outputParam of outputParams.filter(
          (outputParam) => outputParam.property.propertyType === DG.TYPE.DATA_FRAME &&
          (
            (tabLabel === 'Output' && outputParam.property.category === 'Misc' || outputParam.property.category === 'Output') ||
            (tabLabel !== 'Output' && outputParam.property.category === tabLabel)
          )
        )) {
          const nonGridViewers = this.dfToViewerMapping[outputParam.property.name]
            .filter((viewer) => viewer.type !== DG.VIEWER.GRID)
            .filter((viewer) => Object.values(viewerTypesMapping).includes(viewer.type));

          const dfOutput = dfOutputs.find((output) => output.name === outputParam.property.name)!;
          const visibleTitle = dfOutput.options.caption || outputParam.property.name;
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

  exportConfig = {
    supportedExtensions: this.defaultSupportedExportExtensions(),
    supportedFormats: this.defaultSupportedExportFormats(),
    export: this.defaultExport,
    filename: this.defaultExportFilename,
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
