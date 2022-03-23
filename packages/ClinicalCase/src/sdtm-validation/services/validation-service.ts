import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { identicalValues, isNullVariable, negativeValue, nonIdenticalValuesWithExternal, notISO8601, seriousnessCriteriaNotIndicated, startAfterEnd } from '../rules';
import { getColumnNamesWithDomain, getListOfVariablesToValidate, validateColumns } from '../validation-utils';
import { isRequired } from '../../sdtm-meta';
import { seriousnessCriteriaVariables } from '../constants';
import { validatePairOfVariables } from './common-validation-functions';
import { endsWithAny } from '../utils';
import { AE_SERIOUS, DOMAIN } from '../../constants/columns-constants';


export function validateRequiredVariables(df: DG.DataFrame, domain: string, ruleId: string, validationResults: DG.DataFrame){
    const filter = ([k,v]) => v['req'] === isRequired.REQUIRED;
    const presentRequiredVariables = getListOfVariablesToValidate(df, domain, filter);
    presentRequiredVariables.forEach(item => validateColumns([df.getCol(item)], df.rowCount, isNullVariable, 
    validationResults, ruleId, domain, false));
}

export function validateStartEndVariables(df: DG.DataFrame, columnPostfixes: string[], domain: string, ruleId: string, validationResults: DG.DataFrame){
    validatePairOfVariables(df, columnPostfixes, startAfterEnd, validationResults, domain, ruleId)
}

export function validateNonIdenticalVariables(df: DG.DataFrame, columnPostfixes: string[], domain: string, ruleId: string, validationResults: DG.DataFrame){
    validatePairOfVariables(df, columnPostfixes, identicalValues, validationResults, domain, ruleId)
}

export function validateNonNegativeVariables(df: DG.DataFrame, columnPostfixes: string[], domain: string, ruleId: string, validationResults: DG.DataFrame){
    const columnsToValidate = getListOfVariablesToValidate(df, domain, null, getColumnNamesWithDomain(df, columnPostfixes, domain));
    if (columnsToValidate.length){
        columnsToValidate.forEach(item => validateColumns([df.getCol(item)], df.rowCount, negativeValue, 
        validationResults, ruleId, domain, true))
    }    
}

export function validateISO8601Variables(df: DG.DataFrame, columnPostfixes: string[], domain: string, ruleId: string, validationResults: DG.DataFrame){
    const filter = ([k,v]) => endsWithAny(columnPostfixes, k);
    const presentISO8601Variables = getListOfVariablesToValidate(df, domain, filter);
    if(presentISO8601Variables.length){
        presentISO8601Variables.forEach(item => validateColumns([df.getCol(item)], df.rowCount, notISO8601, 
        validationResults, ruleId, domain, true));
    }   
}

export function validateSeriousnesCriteriaVariables(df: DG.DataFrame, domain: string, ruleId: string, validationResults: DG.DataFrame){
    const seriousnesVariable = getListOfVariablesToValidate(df, domain, null, [AE_SERIOUS]);
    const criteriaVariables = getListOfVariablesToValidate(df, domain, null, seriousnessCriteriaVariables);
    if (seriousnesVariable.length){
        validateColumns(seriousnesVariable.concat(criteriaVariables).map(it => df.getCol(it)), df.rowCount, seriousnessCriteriaNotIndicated, 
        validationResults, ruleId, domain, true);
    }
}

export function validateDomainName(df: DG.DataFrame, domain: string, ruleId: string, validationResults: DG.DataFrame){
    const domainVariable = getListOfVariablesToValidate(df, domain, null, [DOMAIN]);
    if(domainVariable.length){
        validateColumns([df.getCol(domainVariable[0])], df.rowCount, nonIdenticalValuesWithExternal, validationResults, ruleId, domain, true, domain);
    }  
}

/* export function validateSubjectIdsInDM(dfDomain: DG.DataFrame, dfDM: DG.DataFrame, domain: string, ruleId: string){
    const subjIdDomainVariable = DG.DataFrame.fromColumns([dfDomain.getCol('USUBJID')]);
    const subjIdDMVariables = DG.DataFrame.fromColumns([dfDomain.getCol('USUBJID')]).columns.addNewInt('dm').init(1);
    const joinedWithDM = grok.data.joinTables(subjIdDomainVariable, subjIdDMVariables, ['USUBJID'], ['USUBJID'], ['USUBJID'], ['USUBJID'], DG.JOIN_TYPE.LEFT, false);

} */

export function vaidateAEDomain(df: DG.DataFrame, validationResults: DG.DataFrame){
    validateRequiredVariables(df, 'ae', 'SD0002', validationResults);
    validateStartEndVariables(df, ['STDY', 'ENDY'], 'ae','SD0012', validationResults);
    validateNonNegativeVariables(df, ['DUR'], 'ae','SD0015', validationResults);
    validateISO8601Variables(df, ['DTC'], 'ae', 'SD0003', validationResults);
    validateISO8601Variables(df, ['DUR', 'ELTM', 'EVLINT', 'STINT', 'ENINT', 'TDSTOFF', 'TDTGTPAI', 'TDMINPAI', 'TDMAXPAI'], 'ae', 'SD1011', validationResults);
    validateSeriousnesCriteriaVariables(df, 'ae', 'SD0009', validationResults);
    validateNonIdenticalVariables(df, ['CAT', 'SCAT'],'ae', 'SD1041', validationResults);
    validateDomainName(df, 'ae', 'SD0004', validationResults);
}

export function vaidateDMDomain(df: DG.DataFrame, validationResults: DG.DataFrame){
    validateRequiredVariables(df, 'dm', 'SD0002', validationResults);
    validateStartEndVariables(df, ['STDY', 'ENDY'], 'dm','SD0012', validationResults);
    validateNonNegativeVariables(df, ['DUR'], 'dm','SD0015', validationResults);
    validateISO8601Variables(df, ['DTC'], 'dm', 'SD0003', validationResults);
    validateISO8601Variables(df, ['DUR', 'ELTM', 'EVLINT', 'STINT', 'ENINT', 'TDSTOFF', 'TDTGTPAI', 'TDMINPAI', 'TDMAXPAI'], 'dm', 'SD1011', validationResults);
    validateDomainName(df, 'dm', 'SD0004', validationResults);
}

