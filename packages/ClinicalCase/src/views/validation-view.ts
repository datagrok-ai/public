import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {_package} from '../package';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import {ValidationResult} from '../types/validation-result';
import {studies} from '../utils/app-utils';

export class ValidationView extends ClinicalCaseViewBase {
  resultsDataframes: DG.DataFrame[];
  rulesDataframe: DG.DataFrame;
  resultsGrid: DG.Grid;
  rulesGrid: DG.Grid;
  domains: any;

  constructor(name: string, studyId: string) {
    super(name, studyId);
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/validation.md`;
  }

  createView(): void {
    if (!studies[this.studyId].validated) {
      ui.setUpdateIndicator(this.root, true, 'Validation in progress...');
      studies[this.studyId].validationCompleted.subscribe(() => {
        ui.setUpdateIndicator(this.root, false);
        this.generateUI();
      });
    } else
      this.generateUI();
  }

  private generateUI() {
    const validationResults: ValidationResult = studies[this.studyId].validationResults;

    if (!validationResults) {
      ui.setUpdateIndicator(this.root, false);
      this.root.appendChild(ui.divText('No validation results available'));
      return;
    }

    const tabControl = ui.tabControl(null, false);

    // Conformance_Details tab - use ui.tableFromMap
    if (validationResults.Conformance_Details) {
      const conformanceMap: {[key: string]: any} = {};
      Object.keys(validationResults.Conformance_Details).forEach((key) => {
        const value = validationResults.Conformance_Details[key as keyof typeof validationResults.Conformance_Details];
        conformanceMap[key] = value !== null && value !== undefined ? String(value) : '';
      });
      const conformanceTable = ui.tableFromMap(conformanceMap);
      conformanceTable.style.width = 'fit-content';
      tabControl.addPane('Conformance_Details', () => conformanceTable);
    }

    // Dataset_Details tab - grid from DataFrame
    if (validationResults.Dataset_Details && validationResults.Dataset_Details.length > 0) {
      const datasetDf = DG.DataFrame.fromObjects(validationResults.Dataset_Details);
      const datasetGrid = datasetDf.plot.grid();
      tabControl.addPane('Dataset_Details', () => datasetGrid.root);
    }

    // Issue_Summary tab - grid from DataFrame
    if (validationResults.Issue_Summary && validationResults.Issue_Summary.length > 0) {
      const issueSummaryDf = DG.DataFrame.fromObjects(validationResults.Issue_Summary);
      const issueSummaryGrid = issueSummaryDf.plot.grid();
      tabControl.addPane('Issue_Summary', () => issueSummaryGrid.root);
    }

    // Issue_Details tab - grid from DataFrame
    if (validationResults.Issue_Details && validationResults.Issue_Details.length > 0) {
      // Convert arrays (variables, values) to strings for display
      const issueDetailsForDf = validationResults.Issue_Details.map((issue) => ({
        ...issue,
        variables: Array.isArray(issue.variables) ? issue.variables.join(', ') : String(issue.variables || ''),
        values: Array.isArray(issue.values) ? issue.values.join(', ') : String(issue.values || ''),
      }));
      const issueDetailsDf = DG.DataFrame.fromObjects(issueDetailsForDf);
      const issueDetailsGrid = issueDetailsDf.plot.grid();
      tabControl.addPane('Issue_Details', () => issueDetailsGrid.root);
    }

    // Rules_Report tab - grid from DataFrame
    if (validationResults.Rules_Report && validationResults.Rules_Report.length > 0) {
      const rulesReportDf = DG.DataFrame.fromObjects(validationResults.Rules_Report);
      const rulesReportGrid = rulesReportDf.plot.grid();
      tabControl.addPane('Rules_Report', () => rulesReportGrid.root);
    }

    tabControl.root.style.width = '100%';

    this.root.appendChild(tabControl.root);
  }
}
