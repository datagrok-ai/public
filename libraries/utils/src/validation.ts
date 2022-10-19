/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import ExcelJS from 'exceljs';

type ValidationSuccess = {isValid: true, chainedFuncs?: ValidatingFunc[], result?: any};
type ValidationFailure = {isValid: false, caption: string, content: HTMLElement | string};
type ValidationResult = ValidationSuccess | ValidationFailure;
export type ValidatingFunc = ((arg: any) => Promise<ValidationResult>);
type InitialValidatingFunc<InitialArgType> = ((arg: InitialArgType) => Promise<ValidationResult>);

class Validation<InitialArgType> {
  constructor(
    // Validation functions to be stored and applied
    public validators: InitialValidatingFunc<InitialArgType>[] = []
  ) {

  }

  public async validate(initialArg: InitialArgType) {
    const failedValidations = [] as ValidationFailure[];
    const checkCustomValidators = async (validators: ValidatingFunc[], result: any) => {
      for (const validator of validators) {
        const validationRes = await validator(result);
        if (!validationRes.isValid)
          failedValidations.push(validationRes);
        else if (validationRes.chainedFuncs) await checkCustomValidators(validationRes.chainedFuncs, validationRes.result);
      }
    };

    for (const validator of this.validators) {
      const validationRes = await validator(initialArg);
      if (!validationRes.isValid)
        failedValidations.push(validationRes);
      else if (validationRes.chainedFuncs) await checkCustomValidators(validationRes.chainedFuncs, validationRes.result);
    }

    if (failedValidations.length == 0) return true;

    const wizard = new DG.Wizard({title: 'Data validation failed'});

    failedValidations.forEach((validation) => {
      wizard.page({
        caption: validation.caption,
        root: (validation.content instanceof HTMLElement ? validation.content : ui.markdown(validation.content)) as HTMLDivElement,
      });
    });
    wizard.onOK(() => {});

    const btnCancel = wizard.getButton('CANCEL');
    const btnBack = wizard.getButton('<<');
    const btnNext = wizard.getButton('>>');
    const btnOk = wizard.getButton('OK');

    const pagesLength = wizard.pages.length;

    const stepCounter = ui.divText('1/'+pagesLength);
    stepCounter.style.position = 'absolute';
    stepCounter.style.left = '12px';
    stepCounter.style.color = 'var(--grey-5)';

    wizard.captionHost.style.fontSize = '20px';
    wizard.captionHost.style.fontWeight = '400';
    wizard.captionHost.style.margin = '6px 0';

    $(wizard.root).find('.d4-command-bar').append(stepCounter);
    $(btnBack).html('BACK');
    $(btnNext).html('NEXT');
    $(btnOk).addClass('disabled');
    btnCancel.style.display = 'none';

    btnBack.addEventListener('click', ()=>{
      const curIndex = wizard.pageIndex+1;
      const pageIndex = wizard.pages.length;

      wizard.captionHost.style.fontSize = '20px';
      wizard.captionHost.style.fontWeight = '400';
      wizard.captionHost.style.margin = '6px 0';

      if (!wizard.completable)
        $(btnOk).addClass('disabled');

      stepCounter.textContent = curIndex+'/'+pageIndex;
    });

    btnNext.addEventListener('click', ()=>{
      const curIndex = wizard.pageIndex+1;
      const pageIndex = wizard.pages.length;

      wizard.captionHost.style.fontSize = '20px';
      wizard.captionHost.style.fontWeight = '400';
      wizard.captionHost.style.margin = '6px 0';

      if (wizard.completable)
        $(btnOk).removeClass('disabled');
      stepCounter.textContent = curIndex+'/'+pageIndex;
    });

    wizard.show();
    return false;
  }
}

export namespace v8n {
  export const isSheetExist = (sheetName: string, chainedFuncs?: ValidatingFunc[]) => async (importedWb: ExcelJS.Workbook) => {
    const datasetWs = importedWb.worksheets.find((ws) => ws.name === sheetName);
    if (!datasetWs) {
      return {
        isValid: false,
        caption: `Sheet named '${sheetName}' is not found`,
        content: `Sheet name should be exactly \`${sheetName}\` to be parsed. 
\n   
Please double-check:
- No extra spaces are in the sheet name
- Capitalization is the same as in the template
\n
${importedWb.worksheets.length > 0 ? `Found sheets: \`${importedWb.worksheets.map((ws) => ws.name).join('`,`')}\``: `No columns have found.`} 
        `
      };
    }
    return {isValid: true as true, chainedFuncs, result: datasetWs};
  };

  export const isColumnExist = (columnName: string, chainedFuncs?: ValidatingFunc[]) => async (sheet: ExcelJS.Worksheet) => {
    const vals = sheet.getRow(1).values as ExcelJS.CellValue[];
    if (!vals.includes(columnName)) {
      return {
        isValid: false,
        caption: `No column "${columnName}" in "${sheet.name}" sheet`,
        content: `Column name should be exactly \`${columnName}\` to be parsed. 
\n
Please double-check:
- Column names are in the first row of the sheet
- No extra spaces are in the cell
- Capitalization is the same as in the template
\n
${vals.length > 0 ? `Found columns: \`${vals.filter((val) => val !== '').join('`,`')}\``: `No columns have found.`} 
        `
      };
    }
    return {isValid: true as true, chainedFuncs};
  };

  export const isExcelParsed = (chainedFuncs?: ValidatingFunc[]) => async (file: File) => {
    const importedWb = new ExcelJS.Workbook();
    return importedWb.xlsx.load(await file.arrayBuffer())
      .then(async (workbook) => {
        return {isValid: true as true, chainedFuncs, result: workbook};
      })
      .catch(() => {
        return {
          isValid: false,
          caption: 'File is failed to parse',
          content: `Please make sure:
- XLSX format is used (not XLS or XLSM)
- File should not be corrupted
\n
XLS could be converted to XLSX using Excel's export functionality. \n
If the file is corrupted, it will not open in Excel. 
          `,
        };
      });
  };
}

export default Validation;
