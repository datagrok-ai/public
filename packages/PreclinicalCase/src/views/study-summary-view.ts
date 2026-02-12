import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BW_RES_N, SEX, SPECIES, SUBJECT_ID,
  VISIT_DAY} from '../constants/columns-constants';
import {removeExtension, studyConfigToMap, updateDivInnerHTML} from '../utils/utils';
import {addDomainAsTableView, studies} from '../utils/app-utils';
import {ACT_TRT_ARM, PLANNED_TRT_ARM} from '../constants/columns-constants';


export class StudySummaryView extends DG.ViewBase {
  validationView: DG.ViewBase | null = null;
  errorsByDomain: any;
  errorsByDomainWithLinks: any;
  studyId: string;
  centralChart = ui.box();
  armChart = ui.box();
  sexChart = ui.box();
  speciesChart = ui.box();
  weightChart = ui.box();
  viewerTitle = {
    style: {
      'color': 'var(--grey-6)',
      'margin': '12px 0px 6px 12px',
      'font-size': '16px',
    },
  };
  validationErrorLinkHandler;
  loaded = false;

  constructor(name: string, studyId: string, errorLinkHandler?: () => void) {
    super();
    this.name = name;
    this.studyId = studyId;
    this.path = '/summary';
    if (errorLinkHandler)
      this.validationErrorLinkHandler = errorLinkHandler;
  }

  load() {
    if (this.loaded)
      return;
    this.loaded = true;
    this.errorsByDomain = studies[this.studyId].errorsByDomain;
    this.errorsByDomainWithLinks = this.createErrorsMapWithLinks();
    this.buildView();
  }

  buildView() {
    const configMap = studyConfigToMap(studies[this.studyId].config);
    const summaryDiv = ui.divV([]);
    summaryDiv.append(ui.h1('Study summary',
      {style: {fontSize: '15px', marginLeft: '12px', position: 'absolute', width: '98%', background: 'white'}}));
    const summaryTableDiv = ui.div(ui.tableFromMap(configMap), {style: {marginTop: '15px'}});
    summaryDiv.append(summaryTableDiv);
    updateDivInnerHTML(this.centralChart, summaryDiv);

    // Treatment arm chart
    const dmDomain = studies[this.studyId].domains.dm;
    const armColName = dmDomain?.col(PLANNED_TRT_ARM) ? PLANNED_TRT_ARM :
      (dmDomain?.col(ACT_TRT_ARM) ? ACT_TRT_ARM : null);
    if (dmDomain && armColName) {
      const arm = DG.Viewer.barChart(dmDomain, //@ts-ignore
        {split: armColName, style: 'dashboard', barColor: DG.Color.lightBlue});
      arm.root.prepend(ui.divText('Treatment arm', this.viewerTitle));
      updateDivInnerHTML(this.armChart, arm.root);
    }

    // Sex chart
    if (dmDomain?.col(SEX)) { //@ts-ignore
      const sex = DG.Viewer.barChart(dmDomain, {split: SEX, style: 'dashboard'});
      sex.root.prepend(ui.divText('Sex', this.viewerTitle));
      updateDivInnerHTML(this.sexChart, sex.root);
    }

    // Species chart
    if (dmDomain?.col(SPECIES)) { //@ts-ignore
      const species = DG.Viewer.barChart(dmDomain, {split: SPECIES, style: 'dashboard'});
      species.root.prepend(ui.divText('Species', this.viewerTitle));
      updateDivInnerHTML(this.speciesChart, species.root);
    }

    // Initial BW distribution
    const visitDayCol = studies[this.studyId].domains.bw?.col(VISIT_DAY) ??
      studies[this.studyId].domains.bw?.col('BWDY');
    const visitDayName = visitDayCol ? visitDayCol.name : VISIT_DAY;
    if (studies[this.studyId].domains.bw?.col(visitDayName) &&
        studies[this.studyId].domains.bw?.col(SUBJECT_ID) &&
        studies[this.studyId].domains.bw?.col(BW_RES_N)) {
      const firstVisit = studies[this.studyId].domains.bw!.col(visitDayName)!.stats.min;
      const df = studies[this.studyId].domains.bw!.groupBy([SUBJECT_ID, visitDayName, BW_RES_N])
        .where(`${visitDayName} = ${firstVisit}`).aggregate();
      const bwHist = DG.Viewer.histogram(df, {value: BW_RES_N});
      bwHist.root.prepend(ui.divText('Initial BW', this.viewerTitle));
      updateDivInnerHTML(this.weightChart, bwHist.root);
    }

    const addDomainToWorkspace = ui.icons.add(() => {
      const menu = DG.Menu.popup();
      for (const domain of studies[this.studyId].domains.all()) {
        const domainName = removeExtension(domain.name);
        menu.item(domainName, () => {
          addDomainAsTableView(this.studyId, domain, false);
        });
      }
      menu.show();
    }, 'Add domain to workspace');
    addDomainToWorkspace.classList.add('preclinical-case-add-domain-to-workspace-icon');

    this.root.className = 'grok-view ui-box';
    const bottomCharts = ui.splitH([
      this.armChart,
      this.sexChart,
      this.speciesChart,
      this.weightChart,
    ]);

    const viewDiv = ui.splitV([]);
    viewDiv.append(ui.div(addDomainToWorkspace,
      {style: {maxHeight: '15px', display: 'flex', alignSelf: 'flex-end'}}));
    viewDiv.append(this.centralChart);
    viewDiv.append(bottomCharts);
    this.root.append(viewDiv);
  }

  private createErrorsMapWithLinks(): {[key: string]: HTMLAnchorElement} | null {
    const errorsByDomain = studies[this.studyId].errorsByDomain;
    if (errorsByDomain) {
      const errorsByDomainWithLinks: {[key: string]: HTMLAnchorElement} = {};
      Object.keys(errorsByDomain).forEach((domain) => {
        const link = ui.link(errorsByDomain[domain].toString(), {}, '', {id: domain});
        link.addEventListener('click', (event) => {
          if (this.validationErrorLinkHandler)
            this.validationErrorLinkHandler();
          else
            grok.shell.v = this.validationView ?? DG.View.create();
          event.stopPropagation();
        });
        errorsByDomainWithLinks[domain] = link;
      });
      return errorsByDomainWithLinks;
    }
    return null;
  }
}
