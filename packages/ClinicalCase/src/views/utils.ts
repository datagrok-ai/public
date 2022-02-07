import { study } from "../clinical-study";
import * as ui from "datagrok-api/ui";
import * as DG from 'datagrok-api/dg';
import { ADVERSE_EVENTS_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME, CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME, SUMMARY_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME, VISITS_VIEW_NAME } from "../view-names-constants";
import * as sdtmCols from "../columns-constants";
import { AE_TERM_FIELD, CON_MED_NAME_FIELD, INV_DRUG_NAME_FIELD, TRT_ARM_FIELD, VIEWS_CONFIG } from "../views-config";

export function updateDivInnerHTML(div: HTMLDivElement, content: any) {
  div.innerHTML = '';
  div.append(content);
}

export function checkRequiredColumns(df: DG.DataFrame, columns: string[], viwerName: string) {
  if (columns.filter(it => !df.columns.names().includes(it)).length) {
    return `The following columns are required for ${viwerName} viewer: ${columns.join(',')}`;
  }
  return null;
}

export function checkColumnsAndCreateViewer(df: DG.DataFrame, columns: string[], div: HTMLDivElement, createViewer: () => any, viewerName: string) {
  const message = checkRequiredColumns(df, columns, viewerName);
  message ? updateDivInnerHTML(div, ui.info(`${message}`)) : createViewer();
}

export function checkMissingDomains(requiredDomainsAndCols: any, obj: any) {
  let loadObject = (obj) => {
    obj.createView();
    obj.loaded = true;
  }

  if (!requiredDomainsAndCols) {
    loadObject(obj);
    return;
  }
  let reqDomains = requiredDomainsAndCols['req_domains'] ? Object.keys(requiredDomainsAndCols['req_domains']) : [];
  let optDomains = requiredDomainsAndCols['opt_domains'] ? Object.keys(requiredDomainsAndCols['opt_domains']) : [];
  let missingReqDomains = reqDomains.filter(it => study.domains[it] === null);
  let missingOptDomains = optDomains.some(it => study.domains[it] !== null) ? [] : optDomains;
  let presentOptDomains = optDomains.filter(it => study.domains[it] !== null);
  let totalMissingDomains = missingReqDomains.concat(missingOptDomains);
  let requiredColumns = {};
  reqDomains.forEach(domain => {
    requiredColumns[domain] = requiredDomainsAndCols['req_domains'][domain];
  });
  optDomains.forEach(domain => {
    requiredColumns[domain] = requiredDomainsAndCols['opt_domains'][domain];
  });
  if (!totalMissingDomains.length) {
    if (checkMissingColumns(obj, reqDomains.concat(presentOptDomains), requiredColumns)) {
      loadObject(obj);
    }
  } else {
    const errorsDiv = ui.divV([], { style: { margin: 'auto', textAlign: 'center' } });
    createMissingDataDiv(errorsDiv, totalMissingDomains, 'Missing domains:');
    checkMissingColumns(errorsDiv, reqDomains.concat(optDomains), requiredColumns, true);
    updateDivInnerHTML(obj.root, errorsDiv);
  }
}

export function checkMissingColumns(obj: any, reqDomains: string[], requiredDomainsAndCols: any, append?: boolean) {
  const errorsDiv = ui.divV([], { style: { margin: 'auto', textAlign: 'center' } });
  let noMissingCols = true;
  reqDomains.forEach(domain => {
    const domainColumns = study.domains[domain] ? study.domains[domain].columns.names() : [];
    const reqCols = requiredDomainsAndCols[domain]['req'] ?? [];
    const optCols = requiredDomainsAndCols[domain]['opt'] ?? []; //at least one of optional columns should exist in domain
    const missingReqColumns = reqCols.filter(it => !domainColumns.includes(it));
    const missingOptColumns = optCols.some(it => study.domains[it] !== null) ? [] : optCols.filter(it => !domainColumns.includes(it));
    const missingColumns = missingReqColumns.concat(missingOptColumns);
    if (missingColumns.length) {
      noMissingCols = false;
      createMissingDataDiv(errorsDiv, missingColumns, `Missing columns in ${domain}:`)
    }
  })
  if (!noMissingCols) {
    let root = obj.root ?? obj;
    append ? root.append(errorsDiv) : updateDivInnerHTML(root, errorsDiv);
  }
  return noMissingCols;
}

export function createMissingDataDiv(div: HTMLDivElement, missingDomainsOrCols: string[], header: string) {
  let domainsDiv = ui.div();
  missingDomainsOrCols.forEach(it => { domainsDiv.append(ui.divText(it)) })
  div.append(ui.div([
    ui.h2(header),
    domainsDiv
  ]));
}

export function getRequiredColumnsByView() {
  // req - all coulmns must be present, opt - at least one of the columns must be present
  // req_domains - all domains must be present, opt_domains - at least one of domains must be present
  return {
    [SUMMARY_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.STUDY_ID,
            sdtmCols.SUBJECT_ID
          ]
        }
      }
    },
    [TIMELINES_VIEW_NAME]: {
      'opt_domains': {
        'ae': {
          'req': [
            sdtmCols.DOMAIN,
            sdtmCols.SUBJECT_ID,
            sdtmCols.AE_START_DAY,
            sdtmCols.AE_END_DAY,
            VIEWS_CONFIG[TIMELINES_VIEW_NAME][AE_TERM_FIELD]
          ],
        },
        'cm': {
          'req': [
            sdtmCols.DOMAIN,
            sdtmCols.SUBJECT_ID,
            VIEWS_CONFIG[TIMELINES_VIEW_NAME][CON_MED_NAME_FIELD],
            sdtmCols.CON_MED_START_DAY,
            sdtmCols.CON_MED_END_DAY
          ]
        },
        'ex': {
          'req': [
            sdtmCols.DOMAIN,
            sdtmCols.SUBJECT_ID,
            VIEWS_CONFIG[TIMELINES_VIEW_NAME][INV_DRUG_NAME_FIELD],
            sdtmCols.INV_DRUG_START_DAY,
            sdtmCols.INV_DRUG_END_DAY
          ]
        }
      }
    },
    [PATIENT_PROFILE_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID
          ]
        }
      }
    },
    [ADVERSE_EVENTS_VIEW_NAME]: {
      'req_domains': {
        'ae': {
          'req': [
          ]
        }
      }
    },
    [LABORATORY_VIEW_NAME]: {
      'req_domains': {
        'lb': {
          'req': [

          ]
        }
      }
    },
    [AE_RISK_ASSESSMENT_VIEW_NAME]: {
      'req_domains': {
        'ae': {
          'req': [
            sdtmCols.SUBJECT_ID,
            VIEWS_CONFIG[AE_RISK_ASSESSMENT_VIEW_NAME][AE_TERM_FIELD]
          ]
        },
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
            VIEWS_CONFIG[AE_RISK_ASSESSMENT_VIEW_NAME][TRT_ARM_FIELD]
          ]
        }
      }
    },
    [SURVIVAL_ANALYSIS_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.SUBJ_REF_STDT,
            sdtmCols.SUBJ_REF_ENDT,
          ]
        }
      }
    },
    [DISTRIBUTIONS_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID
          ],
          'opt': [
            sdtmCols.ETHNIC,
            sdtmCols.SEX,
            sdtmCols.RACE,
            VIEWS_CONFIG[DISTRIBUTIONS_VIEW_NAME][TRT_ARM_FIELD]
          ]
        },
      },
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VISIT_DAY,
            sdtmCols.VISIT_NAME,
            sdtmCols.LAB_RES_N,
            sdtmCols.LAB_TEST
          ]
        },
        'vs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VISIT_DAY,
            sdtmCols.VISIT_NAME,
            sdtmCols.VS_RES_N,
            sdtmCols.VS_TEST
          ]
        }
      }
    },
    [CORRELATIONS_VIEW_NAME]: {
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.LAB_TEST,
            sdtmCols.VISIT_NAME,
            sdtmCols.LAB_RES_N
          ]
        },
        'vs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VS_TEST,
            sdtmCols.VISIT_NAME,
            sdtmCols.VS_RES_N
          ]
        },
      }
    },
    [TIME_PROFILE_VIEW_NAME]: {
      'opt_domains': {
        'lb': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.LAB_TEST,
            sdtmCols.VISIT_NAME,
            sdtmCols.VISIT_DAY,
            sdtmCols.LAB_RES_N
          ]
        },
        'vs': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VS_TEST,
            sdtmCols.VISIT_NAME,
            sdtmCols.VISIT_DAY,
            sdtmCols.VS_RES_N
          ]
        }
      }
    },
    [TREE_MAP_VIEW_NAME]: {
      'req_domains': {
        'dm': {
          'req': [
            sdtmCols.SUBJECT_ID
          ],
          'opt': [
            sdtmCols.ETHNIC,
            sdtmCols.SEX,
            sdtmCols.RACE,
            VIEWS_CONFIG[TREE_MAP_VIEW_NAME][TRT_ARM_FIELD]
          ]
        },
        'ae': {
          'req': [
            sdtmCols.SUBJECT_ID
          ]
        }

      }
    },
    [MEDICAL_HISTORY_VIEW_NAME]: {
      'req_domains': {
        'mh': {
          'req': [
          ]
        }

      }
    },
    [VISITS_VIEW_NAME]: {
      'req_domains': {
        'sv': {
          'req': [
            sdtmCols.SUBJECT_ID,
            sdtmCols.VISIT_START_DATE,
            sdtmCols.VISIT_DAY,
            sdtmCols.VISIT_NAME
          ]
        },

      }
    },
  }
}