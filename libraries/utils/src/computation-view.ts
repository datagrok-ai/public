/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from './function-view';
import dayjs from 'dayjs';

/**
 * Base class for handling Compute models (see https://github.com/datagrok-ai/public/blob/master/help/compute/compute.md).
 * In most cases, the computation is a simple {@link Func}
 * Extend it in cases where a behavior or UI not supported by the {@link FunctionView} is needed.
 *
 * It provides the following functionality out-of-the-box, where each section could be customized:
 * - a structured way to represent input and output parameters: {@link parameters}
 * - generic way to generate UI for inputs, outputs, and interactivity (running the model, etc)
 *   - persisting historical results to the db (via {@link parameters})
 * - export (to Excel and PDF): {@link export}
 * - easy loading of historical runs
 * - routing
 * - entering the real, measured (as opposed to predicted) values manually
 * - notifications for changing inputs, completion of computations, etc: {@link onInputChanged}
 * */
export abstract class ComputationView extends FunctionView {
  /**
   * Overrided to:
   * 1) add getting Help and ReportBug URLs from the package.
   * 2) change routing logic for models
   */
  constructor(
    funcName: string,
    public options: {
      historyEnabled: boolean,
      isTabbed: boolean,
      parentCall?: DG.FuncCall
    } = {historyEnabled: true, isTabbed: false}
  ) {
    super(funcName, options);

    if (!options.parentCall) options.parentCall = grok.functions.getCurrentCall();

    const parentCall = options.parentCall;
    this.parentCall = parentCall;
    this.parentView = parentCall?.parentCall.aux['view'];
    this.basePath = `/${parentCall?.func.name}`;

    this.onFuncCallReady.subscribe({
      complete: async () => {
        await this.getPackageUrls();
        this.buildRibbonMenu();
        this.changeViewName(parentCall.func.friendlyName);
      }
    });
  }

  /** Override to customize getting mocks
    * @stability Stable
  */
  getMocks: ({mockName: string, action: () => Promise<void>}[]) | null = null;

  /** Override to customize getting templates
    * @stability Stable
  */
  getTemplates: ({name: string, action: () => Promise<void>}[]) | null = null;

  /** Override to customize getting help feature. Called when "Help" is clicked.
    * @stability Stable
  */
  getHelp: (() => Promise<void>) | null = null;

  /** Override to customize bug reporting feature. Called when "Report a bug" is clicked.
    * @stability Stable
  */
  reportBug: (() => Promise<void>) | null = null;

  /** Override to customize feature request feature. Called when "Request a feature" is clicked.
    * @stability Stable
  */
  requestFeature: (() => Promise<void>) | null = null;

  /** Override to customize "about" info obtaining feature. Called when "About" is clicked.
    * Default implementation finds {@link this.funcCall}'s package and shows it's properties.
    * @stability Stable
  */
  getAbout: (() => Promise<string>) | null = async () => {
    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12306
    const pack = (await grok.dapi.packages.list()).find((pack) => pack.id === this.func?.package.id);
    return pack ? `${pack.friendlyName} v.${pack.version}.\nLast updated on ${dayjs(pack.updatedOn).format('YYYY MMM D, HH:mm')}`: `No package info was found`;
  };

  /**
   * Looks for
   * {@link getMocks}, {@link getTemplates}, {@link getHelp}, {@link reportBug}, {@link requestFeature}, {@link getAbout}, {@link exportConfig}
   * members and creates "Model" menu
   * @stability Stable
  */
  override buildRibbonMenu() {
    this.ribbonMenu.clear();

    super.buildRibbonMenu();

    if (!this.exportConfig && !this.reportBug && !this.requestFeature && !this.getHelp && !this.getMocks && !this.getTemplates) return;

    const ribbonMenu = this.ribbonMenu.group('Model');

    if (this.getMocks && this.getMocks.length > 0) {
      if (this.getMocks.length === 1) {
        ribbonMenu.item('Input data mock', this.getMocks[0].action);
      } else {
        const dataGroup = ribbonMenu.group('Input data mocks');
        this.getMocks.forEach((val) => {
          dataGroup.item(val.mockName, val.action);
        });
        ribbonMenu.endGroup();
      }
    }

    if (this.getTemplates && this.getTemplates.length > 0) {
      if (this.getTemplates.length === 1) {
        ribbonMenu.item('Input data template', this.getTemplates[0].action);
      } else {
        const dataGroup = ribbonMenu.group('Input data templates');
        this.getTemplates.forEach((val) => {
          dataGroup.item(val.name, val.action);
        });
        ribbonMenu.endGroup();
      }
    }

    if (this.exportConfig && this.exportConfig.supportedFormats.length > 0) {
      ribbonMenu
        .group('Export')
        .items(this.exportConfig.supportedFormats, async (format: string) => DG.Utils.download(this.exportConfig!.filename(format), await this.exportConfig!.export(format)))
        .endGroup();
    }

    if (this.reportBug)
      ribbonMenu.item('Report a bug', () => this.reportBug!());

    if (this.requestFeature)
      ribbonMenu.item('Request a feature', () => this.requestFeature!());

    if (this.getHelp)
      ribbonMenu.item('Help', () => this.getHelp!());

    if (this.getAbout) {
      ribbonMenu.item('About', async () => {
        const dialog = ui.dialog('Current version');
        (await this.getAbout!()).split('\n').forEach((line) => dialog.add(ui.label(line)));
        dialog.onOK(() => {});
        dialog.getButton('CANCEL').style.display = 'none';
        dialog.show({center: true});
      });
    }
  }

  /**
   * Finds {@link this.funcCall}'s package and retrieves it's variables.
   * @stability Stable
  */
  private async getPackageUrls() {
    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12306
    const pack = (await grok.dapi.packages.list()).find((pack) => pack.id === this.parentCall?.func.package.id);

    const reportBugUrl = (await pack?.getProperties() as any).REPORT_BUG_URL;
    if (reportBugUrl && !this.reportBug)
      this.reportBug = async () => { window.open(reportBugUrl, '_blank'); };

    const reqFeatureUrl = (await pack?.getProperties() as any).REQUEST_FEATURE_URL;
    if (reqFeatureUrl && !this.requestFeature)
      this.requestFeature = async () => { window.open(reqFeatureUrl, '_blank'); };
  }
}
