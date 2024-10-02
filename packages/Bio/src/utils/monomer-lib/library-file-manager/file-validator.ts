import {JSONSchemaType, ValidateFunction} from 'ajv';
import Ajv2020 from 'ajv/dist/2020';
import addErrors from 'ajv-errors';

const NA_CODE: string = '#N/A';

import {
  HELM_REQUIRED_FIELD as REQ,
} from '@datagrok-libraries/bio/src/utils/const';

export class MonomerLibFileValidator {
  private validateMonomerSchema: ValidateFunction<any>;

  constructor(
    private helmMonomerSchema: JSONSchemaType<any>
  ) {
    // HELMMonomerSchema.json / #/properties/id uses a union type (string added by Maria Dolotova)
    const ajv = new Ajv2020({allErrors: true, strictTuples: false, allowUnionTypes: true});
    addErrors(ajv);
    this.validateMonomerSchema = ajv.compile(this.helmMonomerSchema);
  }

  validateFile(fileContent: string, fileName: string): boolean {
    const jsonContent = this.parseJson(fileContent, fileName);
    if (jsonContent === null)
      return false;

    if (!Array.isArray(jsonContent)) {
      console.warn(`Bio: Monomer Library File Validator file '${fileName}': Invalid JSON format: ` +
        'The file must contain an array of monomers.');
      return false;
    }

    return this.validateJsonContent(jsonContent, fileName);
  }

  private parseJson(fileContent: string, fileName: string): any[] | null {
    try {
      return JSON.parse(fileContent);
    } catch (e) {
      console.error(`Bio: Monomer Library File Validator file '${fileName}': Invalid JSON format:`, e);
      return null;
    }
  }

  private validateJsonContent(jsonContent: any[], fileName: string): boolean {
    let isValid = true;
    const existingMonomerSymbols = new Set<string>();
    for (const monomer of jsonContent) {
      const name = monomer[REQ.SYMBOL] ?? monomer[REQ.ID] ?? monomer[REQ.NAME] ?? NA_CODE;
      isValid = this.validateMonomerSchema(monomer);
      if (!isValid) {
        console.warn(
          `Bio: Monomer Library File Validator file ${fileName}, monomer '${name}' violating JSON schema:`,
          monomer,
          '\nError reason: ',
          this.validateMonomerSchema.errors,
          `\nThere may be other errors in ${fileName} since the validation is stopped after the first error.`,
          ' Please, verify that the monomer library file satisfies the JSON schema'
        );
        break;
      }
      const key = `${(monomer[REQ.POLYMER_TYPE] ?? '')}-${name}`;
      if (existingMonomerSymbols.has(key)) {
        console.warn(`Bio: Monomer Library File Validator file ${fileName}, monomer '${name}' is duplicated.`,
          'Please, verify that the monomer library file does not contain duplicated monomer symbols.'
        );
      }
      existingMonomerSymbols.add(key);
    }
    return isValid;
  }
}
