import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {SEX, SUBJECT_ID, TSPARMCD, TSVAL, TSPARM} from '../constants/columns-constants';
import {removeExtension} from '../utils/utils';
import {studies} from '../utils/app-utils';

const STUDY_OVERVIEW_CODES = new Set([
  'SPECIES', 'STRAIN', 'SSTYP', 'SDESIGN', 'PLESSION',
]);

const TREATMENT_CODES = new Set([
  'TRT', 'TCNTRL', 'TRTV', 'ROUTE', 'SROTEFQ', 'APTS',
]);

const ADMINISTRATION_CODES = new Set([
  'SPONSOR', 'TESTFAC', 'SSTDDR', 'SREGID', 'SENDVER', 'SNDIGVER', 'PCLAS',
]);

const SKIP_CODES = new Set(['STSTDTC', 'STENDTC']);


export class StudySummaryView extends DG.ViewBase {
  studyId: string;
  loaded = false;
  validationErrorLinkHandler;

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
    this.buildView();
  }

  buildView() {
    const studyDetailsContent = this.buildStudyDetails();
    const reportContent = this.buildReportTab();
    const tabs = ui.tabControl({
      'Study Details': studyDetailsContent,
      'Report': reportContent,
    });
    tabs.root.style.width = '100%';
    tabs.root.style.height = '100%';

    this.root.className = 'grok-view ui-box';
    this.root.append(tabs.root);

    // const generateBtn = ui.button([ui.iconFA('file-alt'), ' Generate Report'], async () => {
    //   grok.shell.info('Report generation is not yet implemented');
    // });
    // generateBtn.classList.add('preclinical-generate-report-btn');
    // this.setRibbonPanels([[generateBtn]]);
  }

  private buildStudyDetails(): HTMLElement {
    const study = studies[this.studyId];
    const config = study.config;
    const container = ui.div([], 'preclinical-study-details');

    const title = document.createElement('h2');
    title.textContent = `Study: ${config.name || this.studyId}`;
    title.style.marginBottom = '4px';
    container.append(title);

    if (config.description) {
      const desc = ui.divText(config.description);
      desc.classList.add('preclinical-study-subtitle');
      container.append(desc);
    }

    const sections = this.parseTsParameters();

    const subjectsStr = this.getSubjectsString();
    if (subjectsStr)
      sections.overview.push(['Subjects', subjectsStr]);

    const startDate = config.startDate ? config.startDate.toString() : null;
    const endDate = config.endDate ? config.endDate.toString() : null;
    if (startDate)
      sections.overview.push(['Start date', startDate]);
    if (endDate)
      sections.overview.push(['End date', endDate]);

    const duration = this.getDuration(sections.overview, startDate, endDate);
    if (duration)
      sections.overview.push(['Duration', duration]);

    if (sections.overview.length)
      this.appendSection(container, 'STUDY OVERVIEW', sections.overview);
    if (sections.treatment.length)
      this.appendSection(container, 'TREATMENT', sections.treatment);
    if (sections.administration.length)
      this.appendSection(container, 'ADMINISTRATION', sections.administration);

    this.appendDomainsSection(container);

    const validationDiv = ui.div([]);
    container.append(validationDiv);
    this.loadValidationSummary(validationDiv);

    return container;
  }

  private parseTsParameters(): {
    overview: [string, string][],
    treatment: [string, string][],
    administration: [string, string][]
  } {
    const overview: [string, string][] = [];
    const treatment: [string, string][] = [];
    const administration: [string, string][] = [];

    const tsDf = studies[this.studyId].domains.ts;
    if (!tsDf) {
      const other = studies[this.studyId].config.other;
      if (other) {
        for (const [key, value] of Object.entries(other))
          overview.push([key, value]);
      }
      return {overview, treatment, administration};
    }

    const codeCol = tsDf.col(TSPARMCD);
    const labelCol = tsDf.col(TSPARM);
    const valCol = tsDf.col(TSVAL);
    if (!codeCol || !valCol)
      return {overview, treatment, administration};

    const unitValues: {[code: string]: string} = {};
    const params: {code: string, label: string, value: string}[] = [];

    for (let i = 0; i < tsDf.rowCount; i++) {
      const code = codeCol.get(i) ?? '';
      const label = labelCol ? (labelCol.get(i) ?? code) : code;
      const value = valCol.get(i) ?? '';

      if (SKIP_CODES.has(code))
        continue;

      if (code.endsWith('U') && code.length > 1) {
        unitValues[code.slice(0, -1)] = value;
        continue;
      }

      params.push({code, label, value});
    }

    for (const param of params) {
      if (unitValues[param.code])
        param.value = `${param.value} ${unitValues[param.code]}`;
    }

    for (const param of params) {
      if (STUDY_OVERVIEW_CODES.has(param.code))
        overview.push([param.label, param.value]);
      else if (TREATMENT_CODES.has(param.code))
        treatment.push([param.label, param.value]);
      else if (ADMINISTRATION_CODES.has(param.code))
        administration.push([param.label, param.value]);
    }

    return {overview, treatment, administration};
  }

  private getSubjectsString(): string {
    const study = studies[this.studyId];
    const dm = study.domains.dm;
    const total = study.config.totalSubjects ||
      (dm ? dm.col(SUBJECT_ID)?.categories?.length ?? dm.rowCount : 0);
    if (!total)
      return '';

    if (dm?.col(SEX)) {
      const sexCol = dm.col(SEX)!;
      const counts: {[key: string]: number} = {};
      for (let i = 0; i < dm.rowCount; i++) {
        const sex = sexCol.get(i);
        if (sex)
          counts[sex] = (counts[sex] || 0) + 1;
      }
      if (Object.keys(counts).length) {
        const parts = Object.entries(counts).map(([sex, count]) => `${count}${sex}`);
        return `${total} (${parts.join(', ')})`;
      }
    }
    return total.toString();
  }

  private getDuration(
    overviewParams: [string, string][],
    startDate: string | null,
    endDate: string | null,
  ): string | null {
    const existing = overviewParams.find(([label]) =>
      label.toLowerCase().includes('length') || label.toLowerCase().includes('duration'));
    if (existing)
      return null;

    if (startDate && endDate) {
      const start = dayjs(startDate);
      const end = dayjs(endDate);
      if (start.isValid() && end.isValid()) {
        const weeks = end.diff(start, 'week');
        if (weeks > 0)
          return `${weeks} weeks`;
        const days = end.diff(start, 'day');
        if (days > 0)
          return `${days} days`;
      }
    }
    return null;
  }

  private appendSection(container: HTMLElement, title: string, items: [string, string][]) {
    const header = ui.divText(title);
    header.classList.add('preclinical-section-header');
    container.append(header);

    const table = document.createElement('table');
    table.classList.add('preclinical-details-table');
    for (const [label, value] of items) {
      const row = table.insertRow();
      const labelCell = row.insertCell();
      labelCell.textContent = label;
      labelCell.classList.add('preclinical-details-label');
      const valueCell = row.insertCell();
      valueCell.textContent = value;
      valueCell.classList.add('preclinical-details-value');
    }
    container.append(table);
  }

  private appendDomainsSection(container: HTMLElement) {
    const domains = studies[this.studyId].domains.all();
    if (!domains.length)
      return;

    const header = ui.divText(`DOMAINS (${domains.length})`);
    header.classList.add('preclinical-section-header');
    container.append(header);

    const tagsContainer = ui.divH([]);
    tagsContainer.classList.add('preclinical-domains-container');
    for (const domain of domains) {
      const name = removeExtension(domain.name);
      const tag = ui.divText(name);
      tag.classList.add('preclinical-domain-tag');
      tagsContainer.append(tag);
    }
    container.append(tagsContainer);
  }

  private async loadValidationSummary(container: HTMLElement) {
    try {
      const path = `System:AppData/Preclinicalcase/SEND/${this.studyId}/validation_summary.json`;
      if (!await grok.dapi.files.exists(path))
        return;
      const json = JSON.parse(await grok.dapi.files.readAsText(path));
      const items: [string, string][] = [];
      if (json.totalIssues != null)
        items.push(['Total issues', `${json.totalIssues}`]);
      if (json.rulesPassed)
        items.push(['Rules passed', json.rulesPassed]);
      if (json.topRule)
        items.push(['Top rule', `${json.topRule.id} (${json.topRule.count} issues, ${json.topRule.percentage}%) \u2014 ${json.topRule.message}`]);
      if (json.mostAffected?.length)
        items.push(['Most affected', json.mostAffected.map((d: any) => `${d.domain} (${d.issues})`).join(', ')]);
      if (json.description)
        items.push(['Description', json.description]);
      if (items.length)
        this.appendSection(container, 'VALIDATION SUMMARY', items);
    }
    catch (_) {
      // silently skip if unavailable
    }
  }

  private buildReportTab(): HTMLElement {
    const container = ui.div([]);
    container.style.height = '100%';
    container.style.width = '100%';
    this.loadReport(container);
    return container;
  }

  private async loadReport(container: HTMLElement) {
    ui.setUpdateIndicator(container, true);
    try {
      const reportPath = `System:AppData/Preclinicalcase/SEND/${this.studyId}/study_report.html`;
      if (await grok.dapi.files.exists(reportPath)) {
        const html = await grok.dapi.files.readAsText(reportPath);
        const blob = new Blob([html], {type: 'text/html'});
        const url = URL.createObjectURL(blob);
        const iframe = document.createElement('iframe');
        iframe.src = url;
        iframe.style.width = '100%';
        iframe.style.height = '100%';
        iframe.style.border = 'none';
        container.append(iframe);
      }
      else
        container.append(ui.divText('No report available.',
          {style: {padding: '20px', color: 'var(--grey-4)'}}));
    }
    catch (e: any) {
      container.append(ui.divText(`Error loading report: ${e?.message ?? e}`,
        {style: {padding: '20px', color: 'var(--red-3)'}}));
    }
    finally {
      ui.setUpdateIndicator(container, false);
    }
  }
}
