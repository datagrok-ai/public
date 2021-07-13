import * as DG from 'datagrok-api/dg';
import * as meta from '../sdtm-meta';


export function createValidationDataFrame(){
    const validationResults = DG.DataFrame.create();
    validationResults.columns.addNewString('Domain');
    validationResults.columns.addNewString('Column');
    validationResults.columns.addNewInt('Row number');
    validationResults.columns.addNewString('Violated rule id');
    return validationResults;
}


export function getColumnNamesWithDomain(df: DG.DataFrame, columns: string[]) {
    const domain = df.getTag('sdtm-domain').toUpperCase();
    return columns.map(col => `${domain}${col}`);
}


export function validateColumns(columns: DG.Column[], 
    rowCount: number, 
    filter: any,
    validationResults: DG.DataFrame, 
    ruleId: string,
    domainId: string,
    getValue: boolean,
    externalConditionVariable?: any){
    for (let i = 0; i < rowCount; ++i) {
        const values = getValue? columns.map((item) => item.get(i)): {column: columns[0], index: i};
        const condition = externalConditionVariable? filter(values, externalConditionVariable): filter(values);
        if(condition){
            updateValidationResult(i, columns[0].name, ruleId, domainId, validationResults);
        }
      }
}


export function updateValidationResult(i, colName, ruleId, domainId, validationResults: DG.DataFrame){
    validationResults.rows.addNew([domainId, colName, i, ruleId]);
}


export function getListOfVariablesToValidate(df: DG.DataFrame, domainId: string, filter: any, columnNames?: string[]){
    const variablesToValidate = filter? 
        Object.entries(meta.domains[domainId]).filter(filter).map(([k,v]) => k):
        columnNames;
    return variablesToValidate.filter(value => df.columns.contains(value));
}