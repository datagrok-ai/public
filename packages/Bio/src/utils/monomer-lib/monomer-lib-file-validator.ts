import Ajv, {JSONSchemaType} from 'ajv-draft-04';
// import draft04 from 'ajv/dist/refs/json-schema-draft-04.json';

export class MonomerLibFileValidator {
  private validateMonomerSchema: (data: any) => boolean;

  constructor(private helmMonomerSchema: JSONSchemaType<any>) {
    const ajv = new Ajv();
    this.validateMonomerSchema = ajv.compile(this.helmMonomerSchema);
  }

  validateFile(fileContent: string): boolean {
    const jsonContent = this.parseJson(fileContent);
    if (jsonContent === null)
      return false;

    return this.validateJsonContent(jsonContent);
  }

  private parseJson(fileContent: string): any[] | null {
    try {
      return JSON.parse(fileContent);
    } catch (e) {
      console.error('Bio: Monomer Library File Validator: Invalid JSON format:', e);
      return null;
    }
  }

  private validateJsonContent(jsonContent: any[]): boolean {
    if (!Array.isArray(jsonContent)) {
      console.error('Bio: Monomer Library File Validator: Invalid JSON format: Expected an array');
      return false;
    }

    return jsonContent.every((monomer) => this.validateMonomerSchema(monomer));
  }
}
