import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ADVERSE_EVENTS_VIEW_NAME, AE_BROWSER_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME,
  ANIMAL_PROFILE_VIEW_NAME, CORRELATIONS_VIEW_NAME,
  DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME,
  QUESTIONNAIRES_VIEW_NAME, SUMMARY_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME,
  TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME, VISITS_VIEW_NAME} from '../constants/view-names-constants';
import * as sdtmCols from '../constants/columns-constants';
import {AE_END_DAY_FIELD, AE_START_DAY_FIELD, AE_TERM_FIELD, CON_MED_END_DAY_FIELD, CON_MED_NAME_FIELD,
  CON_MED_START_DAY_FIELD, INV_DRUG_END_DAY_FIELD, INV_DRUG_NAME_FIELD, INV_DRUG_START_DAY_FIELD,
  TRT_ARM_FIELD} from '../views-config';
import {updateDivInnerHTML} from './utils';
import {CDISC_STANDARD} from './types';
import {VISIT} from '../constants/columns-constants';
import {studies} from './app-utils';
import {IssueDetail} from '../types/validation-result';
import {Subscription} from 'rxjs';

const ERROR_ICON_SIZE = 9;
const ERROR_ICON_MARGIN = 2;
const COLUMNS_WITH_VALIDATION_ERRORS_TAG = 'columnsWithErrors';

export function setupValidationErrorColumns(df: DG.DataFrame) {
  if (df.getTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG))
    return;

  // Find all columns that have corresponding _hasErrors columns
  const columnsWithErrors: {[colName: string]: {hasErrorsCol: string, errorsCol: string}} = {};

  for (const colName of df.columns.names()) {
    if (colName.endsWith(sdtmCols.COL_HAS_ERRORS_POSTFIX)) {
      const baseColName = colName.replace(sdtmCols.COL_HAS_ERRORS_POSTFIX, '');
      const baseCol = df.columns.byName(baseColName);
      const errorsCol = df.columns.byName(`${baseColName}${sdtmCols.ERRORS_POSTFIX}`);

      if (baseCol && errorsCol) {
        // Set semantic type on the errors column to enable custom renderer
        errorsCol.semType = 'sdisc-rule-violation';
        columnsWithErrors[baseColName] = {
          hasErrorsCol: colName,
          errorsCol: `${baseColName}${sdtmCols.ERRORS_POSTFIX}`,
        };
      }
    }
  }

  if (Object.keys(columnsWithErrors).length === 0)
    df.setTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG, '');

  // Store columnsWithErrors on df for tooltip handler access
  df.setTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG, JSON.stringify(columnsWithErrors));
}

// Track current cell with errors state
let currentErrorCell: {
    tableColName: string,
    tableRowIndex: number,
    errors: Array<{ruleID: string, message: string, value: string}>,
    iconX: number,
    iconY: number,
    iconRight: number,
    iconBottom: number
  } | null = null;

export function setupValidationErrorIndicators(grid: DG.Grid, df: DG.DataFrame, ruleId?: string): Subscription[] {
  const columnsWithErrorsString = df.getTag(COLUMNS_WITH_VALIDATION_ERRORS_TAG);
  if (!columnsWithErrorsString)
    return [];

  const columnsWithErrors = JSON.parse(columnsWithErrorsString);


  const onCellRenderedSub = grid.onCellRendered.subscribe((args: DG.GridCellRenderArgs) => {
    const cell = args.cell;
    const tableColName = cell.tableColumn?.name;
    const tableRowIndex = cell.tableRowIndex!;

    if (!cell.isTableCell || !tableColName || tableRowIndex === undefined || tableRowIndex === null)
      return;

    const errorColsNames: {hasErrorsCol: string, errorsCol: string} = columnsWithErrors[tableColName];
    if (!errorColsNames)
      return;
    const hasErrorsCol = df.col(errorColsNames.hasErrorsCol);
    const errorsCol = df.col(errorColsNames.errorsCol);
    if (!hasErrorsCol || !errorsCol)
      return;

    // Check if row index is valid
    if (tableRowIndex < 0 || tableRowIndex >= hasErrorsCol.length)
      return;

    const hasErrors = hasErrorsCol.get(tableRowIndex);

    if (hasErrors) {
      // Parse errors for this specific column
      const errorsStr = errorsCol.get(tableRowIndex);
      let errors: Array<{ruleID: string, message: string, value: string}> = [];
      if (errorsStr) {
        try {
          errors = JSON.parse(errorsStr);
        } catch (e) {
          console.error('Failed to parse column errors:', e);
        }
      }

      if (ruleId) {
        const filteredRule = errors.filter((it) => it.ruleID === ruleId);
        if (!filteredRule.length)
          return;
        errors = filteredRule;
      }

      // Get bounds and canvas context
      const bounds = args.bounds;
      const g = args.g;

      if (!g || !bounds)
        return;

      // Draw error indicator on canvas (similar to Peptides approach)
      g.save();
      g.beginPath();
      g.rect(bounds.x, bounds.y, bounds.width, bounds.height);
      g.clip();

      // Draw error icon in top-right corner with 2px margins
      const iconSize = ERROR_ICON_SIZE;
      const margin = ERROR_ICON_MARGIN;
      const iconX = bounds.x + bounds.width - iconSize - margin;
      const iconY = bounds.y + margin;

      // Draw red circle border (not filled)
      g.strokeStyle = '#dc3545';
      g.lineWidth = 1;
      g.beginPath();
      g.arc(iconX + iconSize / 2, iconY + iconSize / 2, iconSize / 2 - 0.5, 0, 2 * Math.PI);
      g.stroke();

      // Draw red exclamation mark inside
      g.fillStyle = '#dc3545';
      g.font = 'bold 6px Arial';
      g.textAlign = 'center';
      g.textBaseline = 'middle';
      g.fillText('!', iconX + iconSize / 2, iconY + iconSize / 2 + 0.5);

      //args.preventDefault();
      g.restore();
    }
  });

  // Subscribe to cell mouse enter - use hitTest to identify cell and check if it has errors
  const onCellMouseEnterSub = grid.onCellMouseEnter.subscribe((cell: DG.GridCell) => {
    // Use hitTest to identify which cell we're over at the moment
    //const cell = grid.hitTest(lastMouseX, lastMouseY);
    if (!cell || !cell.isTableCell || !cell.tableColumn?.name) {
      currentErrorCell = null;
      return;
    }

    const errorColsNames: {hasErrorsCol: string, errorsCol: string} = columnsWithErrors[cell.tableColumn.name];
    if (!errorColsNames)
      return;
    const hasErrorsCol = df.col(errorColsNames.hasErrorsCol);
    const errorsCol = df.col(errorColsNames.errorsCol);
    if (!hasErrorsCol || !errorsCol) {
      currentErrorCell = null;
      return;
    }

    // Parse errors for this specific column
    const errorsStr = errorsCol.get(cell.tableRowIndex);
    let errors: Array<{ruleID: string, message: string, value: string}> = [];
    if (errorsStr) {
      try {
        errors = JSON.parse(errorsStr);
      } catch (e) {
        console.error('Failed to parse column errors:', e);
        currentErrorCell = null;
        return;
      }
    }

    if (ruleId) {
      const filteredRule = errors.filter((it) => it.ruleID === ruleId);
      if (!filteredRule.length)
        return;
      errors = filteredRule;
    }

    // Get cell bounds and calculate icon position (same as in rendering code)
    const cellBounds = cell.bounds;
    if (!cellBounds) {
      currentErrorCell = null;
      return;
    }

    const iconSize = ERROR_ICON_SIZE;
    const margin = ERROR_ICON_MARGIN;
    const iconX = cellBounds.x + cellBounds.width - iconSize - margin;
    const iconY = cellBounds.y + margin;
    const iconRight = iconX + iconSize;
    const iconBottom = iconY + iconSize;

    // Store current cell with errors and icon bounds
    currentErrorCell = {tableColName: cell.tableColumn!.name,
      tableRowIndex: cell.tableRowIndex, errors, iconX, iconY, iconRight, iconBottom};
  });

  // Subscribe to mouse move on overlay - track position and check coordinates
  grid.overlay.addEventListener('mousemove', handleMouseMoveOverErrorCell);

  // Clear state when mouse leaves a cell
  const onCellMouseLeave = grid.onCellMouseLeave.subscribe(() => {
    currentErrorCell = null;
    ui.tooltip.hide();
  });

  return [onCellRenderedSub, onCellMouseEnterSub, onCellMouseLeave];
}

export function handleMouseMoveOverErrorCell(e: MouseEvent) {
  // Check if we're in a cell with errors
  if (!currentErrorCell) {
    ui.tooltip.hide();
    return;
  }

  // Check if mouse is within icon area (with tolerance)
  const tolerance = 3;
  if (e.offsetX >= currentErrorCell.iconX - tolerance && e.offsetX <= currentErrorCell.iconRight + tolerance &&
        e.offsetY >= currentErrorCell.iconY - tolerance && e.offsetY <= currentErrorCell.iconBottom + tolerance) {
    // Show tooltip
    const tooltipContent = createColumnValidationTooltip(currentErrorCell.errors);
    ui.tooltip.show(tooltipContent, e.clientX, e.clientY);
  } else {
    // Mouse moved outside icon area, hide tooltip
    ui.tooltip.hide();
  }
}

export function createColumnValidationTooltip(
  errors: Array<{ruleID: string, message: string, value: string}>,
): HTMLElement | string {
  if (!errors || errors.length === 0)
    return 'No validation errors';

  const tooltipDiv = ui.div([], {style: {maxWidth: '400px'}});
  errors.forEach((error, index) => {
    const errorDiv = ui.div([],
      {style: {
        marginBottom: '8px',
        paddingBottom: '8px',
        borderBottom: index < errors.length - 1 ? '1px solid var(--grey-2)' : 'none',
      }});

    const ruleIdDiv = ui.div([
      ui.span(['Rule ID: '], {style: {fontWeight: 'bold'}}),
      ui.span([error.ruleID]),
    ]);
    errorDiv.append(ruleIdDiv);

    const messageDiv = ui.div([
      ui.span(['Message: '], {style: {fontWeight: 'bold'}}),
      ui.span([error.message]),
    ]);
    errorDiv.append(messageDiv);

    if (error.value) {
      const valueDiv = ui.div([
        ui.span(['Value: '], {style: {fontWeight: 'bold'}}),
        ui.span([error.value]),
      ]);
      errorDiv.append(valueDiv);
    }

    tooltipDiv.append(errorDiv);
  });

  return tooltipDiv;
}


export function createErrorsByDomainMap(validationResults: DG.DataFrame): {[key: string]: number} {
  const errorsMapWithCount: {[key: string]: number} = {};
  if (validationResults.rowCount) {
    const validationSummary = validationResults.groupBy(['Domain']).count().aggregate();
    for (let i = 0; i < validationSummary.rowCount; ++i)
      errorsMapWithCount[validationSummary.get('Domain', i)] = validationSummary.get('count', i);
    return errorsMapWithCount;
  }
  return null;
}


export function checkRequiredColumns(df: DG.DataFrame, columns: string[], viwerName: string) {
  const missingCols = columns.filter((it) => !df.columns.names().includes(it));
  if (missingCols.length)
    // eslint-disable-next-line max-len
    return `The following columns are required for ${viwerName} viewer: ${columns.join(',')}. Missing ${missingCols.join(',')}`;

  return null;
}

export function checkColumnsAndCreateViewer(df: DG.DataFrame, columns: string[],
  div: HTMLDivElement, createViewer: () => any, viewerName: string): boolean {
  const message = checkRequiredColumns(df, columns, viewerName);
  message ? updateDivInnerHTML(div, ui.info(`${message}`)) : createViewer();
  return !(!!message);
}

export function createValidationErrorsDiv(missingDomains: string[],
  missingColumnsInReqDomains: any, missingColumnsInOptDomains: any) {
  const errorsDiv = ui.divV([], {style: {margin: 'auto', textAlign: 'center'}});
  if (missingDomains.length)
    createMissingDataDiv(errorsDiv, missingDomains, 'Missing domains:');

  createMissingColumnsDiv(missingColumnsInReqDomains, errorsDiv);
  createMissingColumnsDiv(missingColumnsInOptDomains, errorsDiv);
  return errorsDiv;
}


export function createMissingColumnsDiv(domainsWithMissingCols: any, div: HTMLDivElement) {
  Object.keys(domainsWithMissingCols).forEach((domain) => {
    if (domainsWithMissingCols[domain].length)
      createMissingDataDiv(div, domainsWithMissingCols[domain], `Missing columns in ${domain}:`);
  });
}

export function createMissingDataDiv(div: HTMLDivElement, missingDomainsOrCols: string[], header: string) {
  const domainsDiv = ui.div();
  missingDomainsOrCols.forEach((it) => {domainsDiv.append(ui.divText(it));});
  div.append(ui.div([
    ui.h2(header),
    domainsDiv,
  ]));
}

export function getRequiredColumnsByView(studyId: string) {
  // req - all coulmns must be present, opt - at least one of the columns must be present
  // req_domains - all domains must be present, opt_domains - at least one of domains must be present
  return {
    [SUMMARY_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.STUDY_ID,
            sdtmCols.SUBJECT_ID,
          ],
        },
      },
    },
    [TIMELINES_VIEW_NAME]: {
      'opt_domains': {
        'ae': {
          'req': [
            sdtmCols.DOMAIN,
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][AE_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][AE_END_DAY_FIELD],
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][AE_TERM_FIELD],
          ],
        },
        'cm': {
          'req': [
            sdtmCols.DOMAIN,
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][CON_MED_NAME_FIELD],
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][CON_MED_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][CON_MED_END_DAY_FIELD],
          ],
        },
        'ex': {
          'req': [
            sdtmCols.DOMAIN,
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][INV_DRUG_NAME_FIELD],
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][INV_DRUG_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[TIMELINES_VIEW_NAME][INV_DRUG_END_DAY_FIELD],
          ],
        },
      },
    },
    [PATIENT_PROFILE_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
        },
      },
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.LAB_DAY,
            sdtmCols.LAB_TEST,
            sdtmCols.LAB_RES_N,
            sdtmCols.LAB_LO_LIM_N,
            sdtmCols.LAB_HI_LIM_N,
          ],
        },
        'ae': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][AE_TERM_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][AE_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][AE_END_DAY_FIELD],
          ],
        },
        'ex': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][INV_DRUG_NAME_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][INV_DRUG_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][INV_DRUG_END_DAY_FIELD],
          ],
        },
        'cm': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][CON_MED_NAME_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][CON_MED_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][CON_MED_END_DAY_FIELD],
          ],
        },
      },
    },
    [ANIMAL_PROFILE_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
        },
      },
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.LAB_DAY,
            sdtmCols.LAB_TEST,
            sdtmCols.LAB_RES_N,
            sdtmCols.LAB_LO_LIM_N,
            sdtmCols.LAB_HI_LIM_N,
          ],
        },
        'ae': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][AE_TERM_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][AE_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][AE_END_DAY_FIELD],
          ],
        },
        'ex': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][INV_DRUG_NAME_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][INV_DRUG_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][INV_DRUG_END_DAY_FIELD],
          ],
        },
        'cm': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][CON_MED_NAME_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][CON_MED_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[PATIENT_PROFILE_VIEW_NAME][CON_MED_END_DAY_FIELD],
          ],
        },
      },
    },
    [ADVERSE_EVENTS_VIEW_NAME]: {
      'req_domains': {
        'ae': {
          'req': [
          ],
        },
      },
    },
    [LABORATORY_VIEW_NAME]: {
      'req_domains': {
        'lb': {
          'req': [
          ],
        },
      },
    },
    [AE_RISK_ASSESSMENT_VIEW_NAME]: {
      'req_domains': {
        'ae': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[AE_RISK_ASSESSMENT_VIEW_NAME][AE_TERM_FIELD],
          ],
        },
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
            studies[studyId].viewsConfig.config[AE_RISK_ASSESSMENT_VIEW_NAME][TRT_ARM_FIELD],
          ],
        },
      },
    },
    [SURVIVAL_ANALYSIS_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.SUBJ_REF_STDT,
            sdtmCols.SUBJ_REF_ENDT,
          ],
        },
      },
    },
    [DISTRIBUTIONS_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
          'opt': [
            sdtmCols.ETHNIC,
            sdtmCols.SEX,
            sdtmCols.RACE,
            studies[studyId].viewsConfig.config[DISTRIBUTIONS_VIEW_NAME][TRT_ARM_FIELD],
          ],
        },
      },
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            //  sdtmCols.VISIT_DAY,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
            sdtmCols.LAB_RES_N,
            sdtmCols.LAB_TEST,
          ],
        },
        'vs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            //  sdtmCols.VISIT_DAY,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
            sdtmCols.VS_RES_N,
            sdtmCols.VS_TEST,
          ],
        },
      },
    },
    [CORRELATIONS_VIEW_NAME]: {
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.LAB_TEST,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
            sdtmCols.LAB_RES_N,
          ],
        },
        'vs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VS_TEST,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
            sdtmCols.VS_RES_N,
          ],
        },
      },
    },
    [TIME_PROFILE_VIEW_NAME]: {
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.LAB_TEST,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
            // sdtmCols.VISIT_DAY,
            sdtmCols.LAB_RES_N,
          ],
        },
        'vs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VS_TEST,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
            //   sdtmCols.VISIT_DAY,
            sdtmCols.VS_RES_N,
          ],
        },
      },
    },
    [TREE_MAP_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
          'opt': [
            sdtmCols.ETHNIC,
            sdtmCols.SEX,
            sdtmCols.RACE,
            studies[studyId].viewsConfig.config[TREE_MAP_VIEW_NAME][TRT_ARM_FIELD],
          ],
        },
        'ae': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
        },

      },
    },
    [MEDICAL_HISTORY_VIEW_NAME]: {
      'req_domains': {
        'mh': {
          'req': [
          ],
        },

      },
    },
    [VISITS_VIEW_NAME]: {
      'req_domains': {
        'sv': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VISIT_START_DATE,
            sdtmCols.VISIT_DAY,
            studies[studyId].config.standard === CDISC_STANDARD.SEND ? sdtmCols.VISIT_DAY_STR : VISIT,
          ],
        },
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
        },
      },
    },
    [QUESTIONNAIRES_VIEW_NAME]: {
      'req_domains': {
        'qs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.QS_CATEGORY,
            sdtmCols.QS_SUB_CATEGORY,
            sdtmCols.QS_TEST,
            sdtmCols.QS_RES,
          ],
        },
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
          ],
          'opt': [
            sdtmCols.ETHNIC,
            sdtmCols.SEX,
            sdtmCols.RACE,
            studies[studyId].viewsConfig.config[QUESTIONNAIRES_VIEW_NAME][TRT_ARM_FIELD],
          ],
        },
      },
    },
    [AE_BROWSER_VIEW_NAME]: {
      'req_domains': {
        'ae': {
          'req': [
            studies[studyId].viewsConfig.config[AE_BROWSER_VIEW_NAME][AE_TERM_FIELD],
            studies[studyId].viewsConfig.config[AE_BROWSER_VIEW_NAME][AE_START_DAY_FIELD],
            studies[studyId].viewsConfig.config[AE_BROWSER_VIEW_NAME][AE_END_DAY_FIELD],
          ],
        },
      },
    },
  };
}
