/* eslint-disable max-len */
//@ts-ignore
export * from './package.g';

/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {PlateCellHandler} from './plate/plate-cell-renderer';
import {Plate} from './plate/plate';
// import {PlateWidget} from './plate/plate-widget';
import {PlateReader} from './plate/plate-reader';
import {initPlatesAppTree, platesAppView} from './plates/plates-app';
import {initPlates} from './plates/plates-crud';
import {__createDummyPlateData} from './plates/plates-demo';
import {getPlatesFolderPreview} from './plate/plates-folder-preview';
import {PlateTemplateHandler} from './plates/objects/plate-template-handler';
import * as api from './package-api';
import {PlateWidget} from './plate/plate-widget/plate-widget';
import {DrcAnalysis} from './plate/analyses/drc/drc-analysis';
export const _package = new DG.Package();


//tags: autostart
export async function autostart(): Promise<void> {
  await PackageFunctions.createDummyPlateData();
}


export class PackageFunctions {
  @grok.decorators.func({
    name: 'Assay Plates',
    description: 'Assasy plates with concentration, layout and readout data',
    meta: {demoPath: 'Plates | Assay Plates'},
  })
  static async assayPlatesDemo(): Promise<void> {
    const plateFile = (await grok.dapi.files.list('System:DemoFiles/hts/xlsx_plates'))[0];
    grok.shell.addView(await PackageFunctions.previewPlateXlsx(plateFile) as DG.ViewBase);
  }

  @grok.decorators.init()
  static _initPlates(): void {
    DG.ObjectHandler.register(new PlateCellHandler());
    DG.ObjectHandler.register(new PlateTemplateHandler());
  }

  @grok.decorators.folderViewer({outputs: [{'name': 'result', 'type': 'dynamic'}]})
  static async platesFolderPreview(folder: DG.FileInfo, files: DG.FileInfo[]): Promise<DG.Widget | DG.ViewBase | undefined> {
    const nameLowerCase = folder.name?.toLowerCase();
    if (!nameLowerCase?.includes('plate'))
      return undefined;
    return getPlatesFolderPreview(files);
  }

  @grok.decorators.fileViewer({fileViewer: 'txt', fileViewerCheck: 'Plates:checkFileIsPlate'})
  static previewPlate(file: DG.FileInfo): DG.View {
    const view = DG.View.create();
    view.name = file.name;
    file.readAsString().then((content) => {
      const plate = PlateReader.read(content);
      if (plate !== null)
        view.root.appendChild(PlateWidget.fromPlate(plate).root);
    });
    return view;
  }

  @grok.decorators.fileHandler({ext: 'txt', fileViewerCheck: 'Plates:checkFileIsPlate'})
  static async importPlate(fileContent: string): Promise<DG.DataFrame[]> {
    const plate = PlateReader.read(fileContent);
    const view = DG.View.create();
    if (plate !== null)
      view.root.appendChild(PlateWidget.fromPlate(plate).root);

    view.name = 'Plate';
    grok.shell.addView(view);
    return [];
  }

  // NOTE: I commented this out for now because analyses classes have been refactored and PlateDrcAnalysis no longer exists in this form,
  // NOTE: but i am leaving this here to reimplent correctly when i'm adressing xlsx handling (focused solely on csv recently ).
  // @grok.decorators.fileHandler({outputs: [], ext: 'xlsx', fileViewerCheck: 'Curves:checkExcelIsPlate'})
  // static async importPlateXlsx(fileContent: Uint8Array): Promise<any[]> {
  //   const view = DG.View.create();
  //   const plate = await PackageFunctions.parseExcelPlate(fileContent);

  //   const plateWidget = PlateDrcAnalysis.analysisView(plate, {}, 'excel');

  //   if (plateWidget) {
  //     view.root.appendChild(plateWidget.root);
  //   } else {
  //     grok.shell.error('Failed to create plate analysis view. Please check data columns.');
  //     view.close();
  //   }
  //   view.name = 'Plate';
  //   grok.shell.addView(view);
  //   return [];
  // }

  // NOTE: I commented this out for now because analyses classes have been refactored and PlateDrcAnalysis no longer exists in this form,
  // NOTE: but i am leaving this here to reimplent correctly when i'm adressing xlsx handling (focused solely on csv recently ).
  // @grok.decorators.fileViewer({name: 'viewPlateXlsx', fileViewer: 'xlsx', fileViewerCheck: 'Curves:checkExcelIsPlate'})
  // static async previewPlateXlsx(file: DG.FileInfo): Promise<DG.View> {
  //   const view = DG.View.create();
  //   view.name = file.friendlyName;
  //   const plate = await PackageFunctions.parseExcelPlate(await file.readAsBytes());
  //   const plateWidget = PlateDrcAnalysis.analysisView(plate, {}, 'excel');

  //   if (plateWidget) {
  //     view.root.appendChild(plateWidget.root);
  //   } else {
  //     grok.shell.error('Failed to create plate analysis view. Please check data columns.');
  //     view.close();
  //   }
  //   return view;
  // }


@grok.decorators.func({
  name: 'checkExcelIsPlate',
  description: 'Checks if an Excel file contains plate data.'
})
  static async checkExcelIsPlate(content: Uint8Array): Promise<boolean> {
    try {
      if (content.length > 1_000_000) // haven't seen plate files larger than 1MB
        return false;
      const plate = await PackageFunctions.parseExcelPlate(content);
      return plate !== null;
    } catch (e) {
      return false;
    }
  }

static async parseExcelPlate(content: string | Uint8Array, name?: string): Promise<Plate> {
  if (typeof content === 'string') {
    const blob = new Blob([content], {type: 'application/octet-binary'});
    const buf = await blob.arrayBuffer();
    return await Plate.fromExcel(new Uint8Array(buf), name);
  } else {
    return await Plate.fromExcel(content, name);
  }
}

@grok.decorators.fileHandler({outputs: [], ext: 'xlsx', fileViewerCheck: 'Plates:checkExcelIsPlate'})
static async importPlateXlsx(fileContent: Uint8Array): Promise<any[]> {
  const view = DG.View.create();
  const plate = await PackageFunctions.parseExcelPlate(fileContent);

  const plateWidget = PlateWidget.fromPlate(plate);
  const initialMappings = PackageFunctions.autoDetectDrcMappings(plate);
  const drcAnalysis = new DrcAnalysis();
  const analysisView = drcAnalysis.createView(
    plate,
    plateWidget,
    initialMappings,
    (target: string, source: string) => {
      // Handle mapping changes
      initialMappings.set(target, source);
    },
    (target: string) => {
      initialMappings.delete(target);
    }
  );

  const container = ui.divV([
    plateWidget.root,
    analysisView
  ], 'xlsx-plate-container');

  view.root.appendChild(container);
  view.name = 'Plate';
  grok.shell.addView(view);
  return [];
}

// this is a basic solution to make the demo files work. Perhaps this should be a Plate-specific registry of bespoke "heuristics" in the future.
private static autoDetectDrcMappings(plate: Plate): Map<string, string> {
  const mappings = new Map<string, string>();
  const columnNames = plate.data.columns.names();
  const activityCandidates = ['activity', 'response', 'readout', 'value', 'signal', 'raw data'];
  for (const candidate of activityCandidates) {
    const found = columnNames.find((name) => name.toLowerCase().includes(candidate));
    if (found) {
      mappings.set('Activity', found);
      break;
    }
  }

  const concCandidates = ['concentration', 'conc', 'dose', 'concentrations'];
  for (const candidate of concCandidates) {
    const found = columnNames.find((name) => name.toLowerCase().includes(candidate));
    if (found) {
      mappings.set('Concentration', found);
      break;
    }
  }

  const sampleCandidates = ['sample', 'compound', 'id', 'name', 'layout', 'plate layout'];
  for (const candidate of sampleCandidates) {
    const found = columnNames.find((name) => name.toLowerCase().includes(candidate));
    if (found) {
      mappings.set('SampleID', found);
      break;
    }
  }

  return mappings;
}

@grok.decorators.fileViewer({name: 'viewPlateXlsx', fileViewer: 'xlsx', fileViewerCheck: 'Plates:checkExcelIsPlate'})
static async previewPlateXlsx(file: DG.FileInfo): Promise<DG.View> {
  const view = DG.View.create();
  view.name = file.friendlyName;
  const plate = await PackageFunctions.parseExcelPlate(await file.readAsBytes());

  const plateWidget = PlateWidget.fromPlate(plate);

  const initialMappings = PackageFunctions.autoDetectDrcMappings(plate);

  // Create DRC analysis
  const drcAnalysis = new DrcAnalysis();
  const analysisView = drcAnalysis.createView(
    plate,
    plateWidget,
    initialMappings,
    (target: string, source: string) => {
      initialMappings.set(target, source);
    },
    (target: string) => {
      initialMappings.delete(target);
    }
  );

  const container = ui.divV([
    plateWidget.root,
    analysisView
  ], 'xlsx-plate-container');

  view.root.appendChild(container);
  return view;
}

@grok.decorators.func({
  name: 'checkCsvIsPlate',
  description: 'Checks if a CSV file can be parsed as a plate.'
})
static async checkCsvIsPlate(file: DG.FileInfo): Promise<boolean> {
  try {
    const contentSample = await file.readAsString();
    const firstLine = contentSample.substring(0, contentSample.indexOf('\n')).toLowerCase();
    const commonHeaders = ['well', 'position', 'pos'];
    return commonHeaders.some((h) => firstLine.includes(h));
  } catch {
    return false;
  }
}


  @grok.decorators.func()
static checkFileIsPlate(content: string): boolean {
  if (content.length > 1_000_000)
    return false;
  return PlateReader.getReader(content) != null;
}

  @grok.decorators.app({name: 'Browse', browsePath: 'Plates'})
  static platesApp(): DG.View {
    return platesAppView();
  }

  @grok.decorators.func()
  static async platesAppTreeBrowser(treeNode: DG.TreeViewGroup) : Promise<void> {
    await initPlatesAppTree(treeNode);
  }

  @grok.decorators.func()
  static async getPlateByBarcode(barcode: string): Promise<Plate> {
    await initPlates();
    const df: DG.DataFrame = await api.queries.getWellValuesByBarcode(barcode);
    return Plate.fromDbDataFrame(df);
  }

  @grok.decorators.func()
  static async createDummyPlateData(): Promise<void> {
    try {
      await __createDummyPlateData();
    } catch (e) {
      throw e;
    }
  }
}
