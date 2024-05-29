/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {ILineChartSettings} from "datagrok-api/dg";

export class UsageWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    const uniqueUsersDiv = ui.box(null, {style: {margin: '0 12px 0 12px'}});
    const userErrorsDiv = ui.box(null, {style: {margin: '0 12px 0 12px'}});
    const services = ui.box(null, {style: {margin: '0 12px 0 12px'}});
    const link = ui.link('Open Usage Analysis', async () => {
      const progress = DG.TaskBarProgressIndicator.create('Opening Usage Analysis...');
      try {
        grok.shell.addView(await grok.functions.eval('UsageAnalysis:usageAnalysisApp()'));
      } catch (e) {
        console.error(e);
      } finally {
        progress.close();
      }
    });
    const linkDiv = ui.box( ui.div([link],
      {style: {display: 'flex', justifyContent: 'end', alignItems: 'center', paddingRight: '8px'}}), {style: {maxHeight: '40px'}});
    super(ui.box(ui.splitV([linkDiv, uniqueUsersDiv, userErrorsDiv, services], {classes: 'ua-widget'})));

    uniqueUsersDiv.appendChild(ui.waitBox(async () => {
      return ui.splitH([ui.box(ui.divText('Users'),
        {style: {maxWidth: '70px'}}), ui.box(DG.Viewer.fromType('Line chart',
        await grok.functions.call('UsageAnalysis:UniqueUsersSummary'), uniqueUsersChartStyle).root, {style: {paddingRight: '12px'}})]);
    }));

    services.appendChild(ui.waitBox(async () => {
      const infos = (await grok.dapi.admin.getServiceInfos()).map((info) => {
        const d = ui.span([info.name]);
        d.setAttribute('data', info.status);
        d.classList.add('ui-service-info');
        return ui.div([ui.tooltip.bind(d, () => info.status)]);
      });
      const d = ui.div(infos, 'ui-service-info-wrapper');
      return ui.splitH([ui.box(ui.divText('System'),
        {style: {maxWidth: '70px'}}), d]);
    }));

    userErrorsDiv.appendChild(ui.waitBox(async () => {
      return ui.splitH([ui.box(ui.divText('Errors'),
        {style: {maxWidth: '70px'}}), ui.box(DG.Viewer.fromType('Line chart',
        await grok.functions.call('UsageAnalysis:UsersErrorsSummary'), userErrorsChartStyle).root, {style: {paddingRight: '12px'}})]);
    }));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Usage');
    this.order = super.addProperty('order', DG.TYPE.STRING, '3');
  }
}

const uniqueUsersChartStyle: Partial<ILineChartSettings> = {
  'aggrType': 'count',
  'innerChartMarginTop': 0,
  'innerChartMarginBottom': 0,
  'outerChartMarginTop': 0,
  'outerChartMarginBottom': 0,
  'outerChartMarginLeft': 0,
  'outerChartMarginRight': 0,
  'yGlobalScale': false,
  'showTopPanel': false,
  'showMouseOverRowLine': false,
  'showXSelector': false,
  'showYSelectors': false,
  'showAggrSelectors': false,
  'showSplitSelector': false,
  'showYAxis': false,
  'showMarkers': 'Never',
  'autoLayout': false,
};

const userErrorsChartStyle: Partial<ILineChartSettings> = {
  'aggrType': 'count',
  'innerChartMarginTop': 0,
  'innerChartMarginBottom': 0,
  'outerChartMarginTop': 0,
  'outerChartMarginBottom': 0,
  'outerChartMarginLeft': 0,
  'outerChartMarginRight': 0,
  'yGlobalScale': false,
  'showTopPanel': false,
  'showMouseOverRowLine': false,
  'showXSelector': false,
  'showYSelectors': false,
  'showAggrSelectors': false,
  'showSplitSelector': false,
  'showYAxis': false,
  'legendVisibility': DG.VisibilityMode.Never,
  'showMarkers': 'Never',
  'autoLayout': false,
};
