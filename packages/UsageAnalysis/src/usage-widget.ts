/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {usageAnalysisApp} from './package';

export class UsageWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor(header: HTMLDivElement) {
    const icon = ui.iconFA('external-link', (_) => usageAnalysisApp());
    icon.style.marginRight = '6px';
    header.appendChild(icon);
    const uniqueUsersDiv = ui.box();
    const userEventsDiv = ui.box();
    super(ui.box(ui.splitV([uniqueUsersDiv, userEventsDiv], {classes: 'ua-widget'})));

    uniqueUsersDiv.appendChild(ui.wait(async () => {
      return DG.Viewer.fromType('Line chart',
        await grok.data.query('UsageAnalysis:UniqueUsersSummary'), uniqueUsersChartStyle).root;
    }));

    userEventsDiv.appendChild(ui.wait(async () => {
      return DG.Viewer.fromType('Line chart',
        await grok.data.query('UsageAnalysis:UsersEventsSummary'), userEventsChartStyle).root;
    }));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Usage');
    this.order = super.addProperty('order', DG.TYPE.STRING, '2');
  }
}

const uniqueUsersChartStyle = {
  'aggrType': 'count',
  'innerChartMarginTop': 0,
  'innerChartMarginBottom': 0,
  'outerChartMarginTop': 5,
  'outerChartMarginBottom': 0,
  'yGlobalScale': false,
  'showTopPanel': false,
  'showMouseOverRowLine': false,
  'showXSelector': false,
  'showYSelectors': false,
  'showAggrSelectors': false,
  'showSplitSelector': false,
  'showYAxis': false,
  'showMarkers': 'Never',
  'title': 'Unique users',
};

const userEventsChartStyle = {
  'aggrType': 'count',
  'innerChartMarginTop': 0,
  'innerChartMarginBottom': 0,
  'outerChartMarginTop': 5,
  'outerChartMarginBottom': 0,
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
  'title': 'Total events',
};
