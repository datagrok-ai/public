export const ALT = 'ALT';
export const AST = 'AST';
export const BILIRUBIN = 'BLN';
export const AP = 'AP';
export const TREATMENT_ARM = 'ACTARM';
export const SUBJECT_ID = 'USUBJID';
export const STUDY_ID = 'STUDYID';
export const AGE = 'AGE';
export const SEX = 'SEX';
export const RACE = 'RACE';
export const ETHNIC = 'ETHNIC';
export const colorsForSurvivalChart = ['red', 'green', 'blue', 'yellow']
export const SURVIVAL_ANALYSIS_GUIDE = `1. Select dataset paramenters and click 'Create dataset'
2. Filter data if neded,
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