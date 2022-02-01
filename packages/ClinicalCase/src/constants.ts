import * as sdtmCols from "./columns-constants";
import { sdtmSummaryPanel } from "./package";

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
        'Summary': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.STUDY_ID,
                                        sdtmCols.SUBJECT_ID
                                ]
                        }
                }    
        },
        'Timelines': {
                'opt_domains': {
                        'ae': {
                                'req': [
                                        sdtmCols.DOMAIN,
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.AE_START_DAY,
                                        sdtmCols.AE_END_DAY,
                                        sdtmCols.AE_TERM
                                ],
                        },
                        'cm': {
                                'req': [
                                        sdtmCols.DOMAIN,
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.CON_MED_NAME,
                                        sdtmCols.CON_MED_START_DAY,
                                        sdtmCols.CON_MED_END_DAY
                                ]
                        },
                        'ex': {
                                'req': [
                                        sdtmCols.DOMAIN,
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.INV_DRUG_NAME,
                                        sdtmCols.INV_DRUG_START_DAY,
                                        sdtmCols.INV_DRUG_END_DAY
                                ]
                        }
                }
        },
        'Patient Profile': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID  
                                ]
                        },
                        'ae': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.AE_START_DAY,
                                        sdtmCols.AE_END_DAY,
                                        sdtmCols.AE_TERM
                                ]
                        },
                        'cm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.CON_MED_NAME,
                                        sdtmCols.CON_MED_START_DAY,
                                        sdtmCols.CON_MED_END_DAY
                                ]
                        },
                        'ex': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.INV_DRUG_NAME,
                                        sdtmCols.INV_DRUG_START_DAY,
                                        sdtmCols.INV_DRUG_END_DAY,
                                        sdtmCols.INV_DRUG_DOSE,
                                        sdtmCols.INV_DRUG_DOSE_UNITS
                                ]
                        },
                        'lb': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.LAB_DAY,
                                        sdtmCols.LAB_HI_LIM_N,
                                        sdtmCols.LAB_LO_LIM_N,
                                        sdtmCols.LAB_RES_N,
                                        sdtmCols.LAB_TEST
                                ]
                        }
                }
        },
        'Adverse Events': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.TREATMENT_ARM,
                                ]
                        },
                        'ae': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.AE_DECOD_TERM,
                                        sdtmCols.AE_BODY_SYSTEM,
                                        sdtmCols.AE_CAUSALITY,
                                        sdtmCols.AE_OUTCOME,
                                        sdtmCols.AE_START_DAY,
                                        sdtmCols.AE_SEVERITY
                                ]
                        }
                }
        },
        'Laboratory': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.TREATMENT_ARM
                                ]
                        },
                        'lb': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.LAB_TEST,
                                        sdtmCols.VISIT_NAME,
                                        sdtmCols.VISIT_DAY,
                                        sdtmCols.LAB_RES_N,
                                        sdtmCols.LAB_LO_LIM_N,
                                        sdtmCols.LAB_HI_LIM_N
                                ]
                        }
                }
        },
        'AE Risk Assessment': {
                'req_domains': {
                        'ae': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.AE_TERM
                                ]
                        },
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.TREATMENT_ARM
                                ]
                        }
                }
        },
        'Survival Analysis': {
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
        'Distributions': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID
                                ],
                                'opt': [
                                        sdtmCols.ETHNIC,
                                        sdtmCols.SEX,
                                        sdtmCols.RACE,
                                        sdtmCols.TREATMENT_ARM,
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
        'Correlations': {
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
        'Time Profile': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                ],
                                'opt': [
                                        sdtmCols.ETHNIC,
                                        sdtmCols.SEX,
                                        sdtmCols.RACE,
                                        sdtmCols.TREATMENT_ARM
                                ]
        
                        },
                },
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
        'Tree map': {
                'req_domains': {
                        'dm': {
                                'req': [
                                        sdtmCols.SUBJECT_ID
                                ],
                                'opt': [
                                        sdtmCols.ETHNIC,
                                        sdtmCols.SEX,
                                        sdtmCols.RACE,
                                        sdtmCols.TREATMENT_ARM,
                                ]
                        },
                        'ae': {
                                'req': [
                                        sdtmCols.SUBJECT_ID
                                ]
                        }

                }
        },
        'Medical History': {
                'req_domains': {
                        'mh': {
                                'req': [
                                        sdtmCols.SUBJECT_ID,
                                        sdtmCols.MH_DECOD_TERM,
                                        sdtmCols.MH_CATEGORY,
                                        sdtmCols.MH_BODY_SYSTEM,
                                        sdtmCols.MH_TERM
                                ]
                        }

                }
        },
        'Visits': {
                'req_domains': {
                        'tv': {
                                'req': [
                                        sdtmCols.VISIT_DAY,
                                        sdtmCols.VISIT_NAME
                                ]
                        },
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