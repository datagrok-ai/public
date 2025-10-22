// src/plate/plates-folder-preview.ts

/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Plate} from './plate';
import {_package} from '../package';
import {DrcAnalysis} from './analyses/drc/drc-analysis';
import {AnalysisManager} from './analyses/analysis-manager';
import {PlateWidget} from './plate-widget/plate-widget';
// import {PlateWidget} from './plate-widget';


export async function getPlatesFolderPreview(files: DG.FileInfo[]): Promise<DG.Widget | DG.ViewBase | undefined> {
  const csvFiles = files.filter((f) => f?.name?.toLowerCase()?.endsWith('.csv'));
  let csvView: DG.Widget | undefined = undefined;

  if (csvFiles.length > 2) {
    const plate = Plate.fromPlates(await Promise.all(csvFiles.map(async (f) => await Plate.fromCsvTableFile(f.fullPath, f.name.toLowerCase().substring(0, f.name.length - 4)))));

    const drcAnalysis = await AnalysisManager.instance.getAnalysis('DRC') as DrcAnalysis;
    if (drcAnalysis) {
      const plateWidget = PlateWidget.fromPlate(plate);
      const mappings = new Map<string, string>([

      ]);
      const analysisView = drcAnalysis.createView(plate, plateWidget, mappings, () => {}, () => {});
      plateWidget.addAnalysisTab('DRC', analysisView);
      csvView = plateWidget;
      plate.data.name = `${csvFiles.map((file) => file.name.slice(0, file.name.length - 4)).join('_')}.csv`;
      if (csvFiles.length === files.length)
        return csvView;
    }
  }
  const multiView = new DG.MultiView({viewFactories: {}});
  if (csvView)
    multiView.addView('Plate 1', () => DG.View.fromRoot(csvView!.root), true);
  const xlsxFiles = files.filter((f) => f?.name?.toLowerCase()?.endsWith('.xlsx') && f?.name?.toLowerCase().includes('plate'));
  if (xlsxFiles.length == 0 && !csvView)
    return undefined;
  else if (xlsxFiles.length == 0 && csvView)
    return csvView;
  for (const xlsxFile of xlsxFiles) {
    try {
      const plate = await Plate.fromExcelFile(xlsxFile);
      const drcAnalysis = await AnalysisManager.instance.getAnalysis('DRC') as DrcAnalysis;
      if (drcAnalysis) {
        const plateWidget = PlateWidget.fromPlate(plate);
        const mappings = new Map<string, string>([]);
        const analysisView = drcAnalysis.createView(plate, plateWidget, mappings, () => {}, () => {});
        plateWidget.addAnalysisTab('DRC', analysisView);

        const v = DG.View.fromRoot(plateWidget.root);
        v.name = xlsxFile.name.substring(0, xlsxFile.name.length - 5);
        multiView.addView(v.name, () => v, true);
      }
    } catch (e) {
      _package.logger.error(e);
    }
  }
  setTimeout(() => {
    const header: HTMLElement | null = multiView.root.querySelector('.d4-tab-header-stripe');
    if (!header)
      return;
    header.style.overflow = 'scroll';

    const headerItems = header.querySelectorAll('.d4-tab-header');
    headerItems?.forEach((el) => (el as HTMLElement).style.whiteSpace = 'nowrap');
  }, 300);

  return multiView;
}
