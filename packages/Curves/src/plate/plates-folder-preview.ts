import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Plate} from './plate';
//@ts-ignore
import * as jStat from 'jstat';
import {FIT_FUNCTION_4PL_REGRESSION, FIT_FUNCTION_SIGMOID, FitMarkerType, IFitPoint} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import { FitConstants } from '../fit/const';
import {FitCellOutlierToggleArgs, setOutlier} from '../fit/fit-renderer';
import { _package } from '../package';
import {savePlate} from "../plates/plates-crud";
import { PlateWidget } from './plate-widget';
import { PlateDrcAnalysis } from './plate-drc-analysis';



export async function getPlatesFolderPreview(files: DG.FileInfo[]): Promise<DG.Widget | DG.ViewBase | undefined> {

  const csvFiles = files.filter((f) => f?.name?.toLowerCase()?.endsWith('.csv'));
  let csvView: DG.Widget | undefined = undefined;
  if (csvFiles.length > 2) {

    const plate = Plate.fromPlates(await Promise.all(csvFiles.map(async (f) => await Plate.fromCsvTableFile(f.fullPath, f.name.toLowerCase().substring(0, f.name.length - 4)))));
    csvView = PlateDrcAnalysis.analysisView(plate, {submitAction: () => {grok.shell.info('Plate Submitted')}});
    plate.data.name = `${csvFiles.map((file) => file.name.slice(0, file.name.length - 4)).join('_')}.csv`;
    if (csvFiles.length === files.length)
      return csvView;
  }
  const multiView = new DG.MultiView({viewFactories: {}});

  if (csvView) {
    multiView.addView('Plate 1', () => DG.View.fromRoot(csvView!.root), true);
  }

  const xlsxFiles = files.filter((f) => f?.name?.toLowerCase()?.endsWith('.xlsx') && f?.name?.toLowerCase().includes('plate'));

  if (xlsxFiles.length == 0)
    return csvView;

  for (const xlsxFile of xlsxFiles) {
    try {
      const plate = await Plate.fromExcelFile(xlsxFile);
      const pw = PlateDrcAnalysis.analysisView(plate, {submitAction: () => {grok.shell.info('Plate Submitted')}});
      const v = DG.View.fromRoot(pw.root);
      v.name = xlsxFile.name.substring(0, xlsxFile.name.length - 5);
      multiView.addView(v.name, () => v, true);
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
  }, 300)

  return multiView;
}