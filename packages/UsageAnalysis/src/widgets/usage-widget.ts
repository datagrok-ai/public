/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import {usageAnalysisApp} from './package';

export class UsageWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    const uniqueUsersDiv = ui.box(null, {style: {margin: '0 12px 0 12px'}});
    const userEventsDiv = ui.box(null, {style: {margin: '0 12px 0 12px'}});
    const userErrorsDiv = ui.box(null, {style: {margin: '0 12px 0 12px'}});
    const link = ui.link('Open Usage Analysis', () => grok.functions.eval('UsageAnalysis:usageAnalysisApp()'));
    const linkDiv = ui.box( ui.div([link],
      {style: {display: 'flex', justifyContent: 'end', alignItems: 'center', paddingRight: '8px'}}), {style: {maxHeight: '40px'}});
    super(ui.box(ui.splitV([linkDiv, uniqueUsersDiv, userEventsDiv, userErrorsDiv], {classes: 'ua-widget'})));

    uniqueUsersDiv.appendChild(ui.waitBox(async () => {
      return ui.splitH([ui.box(ui.divText('Users'),
        {style: {maxWidth: '70px'}}), ui.box(DG.Viewer.fromType('Line chart',
        await grok.functions.call('UsageAnalysis:UniqueUsersSummary'), uniqueUsersChartStyle).root, {style: {paddingRight: '12px'}})]);
    }));

    userEventsDiv.appendChild(ui.waitBox(async () => {
      return ui.splitH([ui.box(ui.divText('Events'),
        {style: {maxWidth: '70px'}}), ui.box(DG.Viewer.fromType('Line chart',
        await grok.functions.call('UsageAnalysis:UsersEventsSummary'), userEventsChartStyle).root, {style: {paddingRight: '12px'}})]);
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

const uniqueUsersChartStyle = {
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

const userEventsChartStyle = {
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
  'legendVisibility': 'Never',
  'showMarkers': 'Never',
  'autoLayout': false,
};


const userErrorsChartStyle = {
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
  'legendVisibility': 'Never',
  'showMarkers': 'Never',
  'autoLayout': false,
};
