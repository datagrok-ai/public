import * as DG from 'datagrok-api/dg';
import { notISO8601 } from '../rules';
import { getColumnNamesWithDomain, getListOfVariablesToValidate, validateColumns } from '../validation-utils';

export function validatePairOfVariables(df: DG.DataFrame, columnPostfixes: string[], rule: any, validationResults: DG.DataFrame, domain: string, ruleId: string){
    const columnsToCompare = getListOfVariablesToValidate(df, domain, null, getColumnNamesWithDomain(df, columnPostfixes, domain));
    if (columnsToCompare.length === 2) {
        validateColumns(columnsToCompare.map(it => df.getCol(it)), df.rowCount, rule, validationResults, ruleId, domain, true);
    }
}