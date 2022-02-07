import * as sdtmCols from "./columns-constants";
import { ADVERSE_EVENTS_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME, CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME, PATIENT_PROFILE_VIEW_NAME, SUMMARY_VIEW_NAME, SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME, VISITS_VIEW_NAME } from "./view-names-constants";
import { AE_TERM_FIELD, CON_MED_NAME_FIELD, INV_DRUG_NAME_FIELD, TRT_ARM_FIELD, VIEWS_CONFIG } from "./views-config";

export const ALT = 'ALT';
export const AST = 'AST';
export const BILIRUBIN = 'BLN';
export const AP = 'AP';
export const colorsForSurvivalChart = ['red', 'green', 'blue', 'yellow'];
export const SURVIVAL_ANALYSIS_GUIDE = `1. Select dataset parameters and click 'Create dataset'.
2. Filter data if needed.
3. Go to 'Survival Chart' tab to see Kaplan-Meier curve. You can modify the survival chart parameters - confidence interval and strata.
4. Go to 'Covariates' tab to perform covariates analysis. Check one or more covariates to see graphs.
`
export const SEVERITY_COLOR_DICT = {
        'MILD': 'var(--green-2)', 
        'MODERATE': 'var(--orange-2)', 
        'SEVERE': 'var(--red-3)', 
        'U': 'var(--grey-3)', 
        'UNK': 'var(--grey-3)', 
        'Unknown': 'var(--grey-3)'}

export const CLINICAL_TRIAL_GOV_FIELDS = {
        'NCTId': 'NCT ID',
        'Condition': 'Condition',
        'BriefTitle': 'Brief Title',
        'OrgFullName': 'Organization',
        'Phase': 'Phase',
        'StartDate': 'Start Date',
        'CompletionDate':'Completion Date',
        'BriefSummary': 'Brief Summary'
}

export const ACTIVE_ARM_POSTTFIX = 'ACTIVE';
export const PLACEBO_ARM_POSTTFIX = 'PLACEBO';
export const RELATIVE_RISK = 'RR';
export const RISK_DIFFERENCE = 'RD';
export const ODDS_RATIO = 'OR';
export const P_VALUE = 'p-value';
export const NEG_LOG10_P_VALUE = '-log10(p-value)';
export const AE_PERCENT = 'AE_PERCENT';
export const STANDARD_ERROR_RD = 'SE_RD';
export const SE_RD_WITH_SIGN_LEVEL = 'SE_RD_SIGN_LEVEL';

// req - all coulmns must be present, opt - at least one of the columns must be present
// req_domains - all domains must be present, opt_domains - at least one of domains must be present

export const requiredColumnsByView = {
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