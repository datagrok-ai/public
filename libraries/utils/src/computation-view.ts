/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from './function-view';
import '../css/computation-view.css';
import dayjs from 'dayjs';
import {historyUtils} from './history-utils';

const url = new URL(grok.shell.startUri);
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
export class ComputationView extends FunctionView {
  /** Find the function by fully specified name, link it to the view and constructs the view.
    * If function name is not specified, calls {@link init} and {@link build} without FuncCall linkage.
    * @stability Stable
  */
  constructor(funcName?: string) {
    const runId = url.searchParams.get('id');
    super();

    this.parentCall = grok.functions.getCurrentCall();
    this.parentView = grok.functions.getCurrentCall().parentCall.aux['view'];
    this.basePath = `/${grok.functions.getCurrentCall()?.func.name}`;

    ui.setUpdateIndicator(this.root, true);
    if (runId) {
      setTimeout(async () => {
        await this.init();
        await this.getPackageUrls();
        const urlRun = await historyUtils.loadRun(runId);
        this.linkFunccall(urlRun);
        this.build();
        await this.onBeforeLoadRun();
        this.linkFunccall(urlRun);
        await this.onAfterLoadRun(urlRun);
        this.setRunViewReadonly();
        ui.setUpdateIndicator(this.root, false);
        url.searchParams.delete('id');
      }, 0);
    } else {
      if (funcName) {
        (grok.functions.eval(funcName) as Promise<DG.Func>).then(async (func: DG.Func) => {
          const funccall = func.prepare({});
          this.linkFunccall(funccall);
          await this.init();
          await this.getPackageUrls();
          this.build();
        }).finally(() => {
          ui.setUpdateIndicator(this.root, false);
        });
      } else {
        setTimeout(async () => {
          await this.init();
          await this.getPackageUrls();
          this.build();
          this.buildRibbonPanels();
          ui.setUpdateIndicator(this.root, false);
        }, 0);
      }
    }
  }

  /** Override to customize getting mocks
    * @stability Experimental
  */
  getMocks: ({mockName: string, action: () => Promise<void>}[]) | null = null;

  /** Override to customize getting templates
    * @stability Experimental
  */
  getTemplates: ({name: string, action: () => Promise<void>}[]) | null = null;

  /** Override to customize getting help feature
    * @stability Stable
  */
  getHelp: (() => Promise<void>) | null = null;

  /** Override to customize bug reporting feature
    * @stability Stable
  */
  reportBug: (() => Promise<void>) | null = null;

  /** Override to customize feature request feature
    * @stability Stable
  */
  requestFeature: (() => Promise<void>) | null = null;

  /** Override to customize "about" info obtaining feature.
    * @stability Experimental
  */
  getAbout: (() => Promise<string>) | null = async () => {
    const pack = (await grok.dapi.packages.list()).find((pack) => pack.id === this.func?.package.id);
    return pack ? `${pack.friendlyName} v.${pack.version}.\nLast updated on ${dayjs(pack.updatedOn).format('YYYY MMM D, HH:mm')}`: `No package info was found`;
  };

  /**
   * Looks for {@link reportBug}, {@link getHelp} and {@link exportConfig} members and creates model menus
   * @stability Stable
  */
  override buildRibbonMenu() {
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

  private async getPackageUrls() {
    const pack = (await grok.dapi.packages.list()).find((pack) => pack.id === this.func?.package.id);
    const reportBugUrl = (await pack?.getProperties() as any).REPORT_BUG_URL;
    if (reportBugUrl && !this.reportBug)
      this.reportBug = async () => { window.open(reportBugUrl, '_blank'); };

    const reqFeatureUrl = (await pack?.getProperties() as any).REQUEST_FEATURE_URL;
    if (reqFeatureUrl && !this.requestFeature)
      this.requestFeature = async () => { window.open(reqFeatureUrl, '_blank'); };
  }
}
