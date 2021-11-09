import * as sdtmCols from "./columns-constants";
import { sdtmSummaryPanel } from "./package";

export const ALT = 'ALT';
export const AST = 'AST';
export const BILIRUBIN = 'BLN';
export const AP = 'AP';
export const colorsForSurvivalChart = ['red', 'green', 'blue', 'yellow']
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

export const requiredColumnsByView = {
        'Summary': {
                'dm': [
                        sdtmCols.STUDY_ID,
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AGE,
                        sdtmCols.SEX,
                        sdtmCols.RACE,
                        sdtmCols.TREATMENT_ARM,
                        sdtmCols.SUBJ_REF_STDT
                ]
        },
        'Timelines': {
                'ae': [
                        sdtmCols.DOMAIN,
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AE_START_DAY,
                        sdtmCols.AE_END_DAY,
                        sdtmCols.AE_TERM,
                        sdtmCols.AE_SEVERITY,
                        sdtmCols.AE_BODY_SYSTEM
                ],
                'cm': [
                        sdtmCols.DOMAIN,
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.CON_MED_NAME,
                        sdtmCols.CON_MED_START_DAY,
                        sdtmCols.CON_MED_END_DAY
                ],
                'ex': [
                        sdtmCols.DOMAIN,
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.INV_DRUG_NAME,
                        sdtmCols.INV_DRUG_START_DAY,
                        sdtmCols.INV_DRUG_END_DAY
                ]
        },
        'Patient Profile': {
                'dm': [
                        sdtmCols.SUBJECT_ID  
                ],
                'ae': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AE_START_DAY,
                        sdtmCols.AE_END_DAY,
                        sdtmCols.AE_TERM
                ],
                'cm': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.CON_MED_NAME,
                        sdtmCols.CON_MED_START_DAY,
                        sdtmCols.CON_MED_END_DAY
                ],
                'ex': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.INV_DRUG_NAME,
                        sdtmCols.INV_DRUG_START_DAY,
                        sdtmCols.INV_DRUG_END_DAY,
                        sdtmCols.INV_DRUG_DOSE,
                        sdtmCols.INV_DRUG_DOSE_UNITS
                ],
                'lb': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.LAB_DAY,
                        sdtmCols.LAB_HI_LIM_N,
                        sdtmCols.LAB_LO_LIM_N,
                        sdtmCols.LAB_RES_N,
                        sdtmCols.LAB_TEST,
                        sdtmCols.LAB_VISIT_DAY
                ]
        },
        'Adverse Events': {
                'dm': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.TREATMENT_ARM,
                ],
                'ae': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AE_DECOD_TERM,
                        sdtmCols.AE_BODY_SYSTEM,
                        sdtmCols.AE_CAUSALITY,
                        sdtmCols.AE_OUTCOME,
                        sdtmCols.AE_START_DAY,
                        sdtmCols.AE_SEVERITY
                ]
        },
        'Laboratory': {
                'dm': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.TREATMENT_ARM
                ],
                'lb': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.LAB_TEST,
                        sdtmCols.LAB_VISIT_NAME,
                        sdtmCols.LAB_VISIT_DAY,
                        sdtmCols.LAB_RES_N,
                        sdtmCols.LAB_LO_LIM_N,
                        sdtmCols.LAB_HI_LIM_N
                ]
        },
        'AE Risk Assessment': {
                'ae': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AE_TERM
                ],
                'ex': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.INV_DRUG_NAME
                ]
        },
        'Survival Analysis': {
                'dm': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AGE,
                        sdtmCols.SEX,
                        sdtmCols.RACE,
                        sdtmCols.TREATMENT_ARM,
                        sdtmCols.DEATH_DATE 
                ],
                'ae': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.AE_START_DATE,
                        sdtmCols.AE_REQ_HOSP,
                        sdtmCols.AE_SEVERITY,
                        sdtmCols.AE_CAUSALITY,
                        sdtmCols.AE_SEQ
                ]
        },
        'Biomarkers': {
                'dm': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.ETHNIC,
                        sdtmCols.SEX,
                        sdtmCols.RACE,
                        sdtmCols.TREATMENT_ARM,
                ],
                'lb': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.LAB_VISIT_DAY,
                        sdtmCols.LAB_VISIT_NAME,
                        sdtmCols.LAB_RES_N,
                        sdtmCols.LAB_TEST
                ]
        },
        'Correlations': {
                'lb': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.LAB_TEST,
                        sdtmCols.LAB_VISIT_NAME,
                        sdtmCols.LAB_RES_N
                ]
        },
        'Time Profile': {
                'dm': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.ETHNIC,
                        sdtmCols.SEX,
                        sdtmCols.RACE,
                        sdtmCols.TREATMENT_ARM,
                ],
                'lb': [
                        sdtmCols.SUBJECT_ID,
                        sdtmCols.LAB_TEST,
                        sdtmCols.LAB_VISIT_NAME,
                        sdtmCols.LAB_VISIT_DAY,
                        sdtmCols.LAB_TEST,
                        sdtmCols.LAB_RES_N
                ]
        }
}