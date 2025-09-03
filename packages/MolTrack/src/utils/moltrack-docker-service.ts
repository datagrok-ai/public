import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, ResultOutput, scopeToUrl } from './constants';
import { MolTrackProperty, MolTrackSearchQuery, MolTrackSearchResponse} from './types';

export class MolTrackDockerService {
  private static _container: DG.DockerContainer;

  static async init(): Promise<void> {
    this._container = await grok.dapi.docker.dockerContainers
      .filter('name = "moltrack"')
      .first();
  }

  private static get container(): DG.DockerContainer {
    if (!this._container)
      throw new Error('MolTrackDockerService not initialized. Call init() first.');

    return this._container;
  }

  static async postToEndpoint(endpoint: string, jsonPayload: string): Promise<string> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, endpoint, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: jsonPayload,
    });
    return response.text();
  }

  static async checkHealth(): Promise<string> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, '/v1/health');
    return response.text();
  }

  static async fetchCompoundProperties(): Promise<string> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, '/v1/schema/compounds');
    return response.text();
  }

  static async fetchBatchProperties(): Promise<string> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, '/v1/schema/batches');
    return response.text();
  }

  static async fetchSchema(): Promise<MolTrackProperty[]> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, '/v1/schema/');
    if (!response.ok)
      throw new Error(`HTTP error!: ${response.status}`, { cause: response.status });
    const res = await response.json();
    return res;
  }


  static async fetchDirectSchema(): Promise<Partial<MolTrackProperty>[]> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id, '/v1/schema-direct/');
    const data = await response.json();
    return data.flat();
  }

  static async search(query: MolTrackSearchQuery, endpointLevel: string): Promise<MolTrackSearchResponse> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(this.container.id,
      `/v1/search/${endpointLevel}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(query),
      });
    if (!response.ok)
      throw new Error(`HTTP error!: ${response.status}`, { cause: response.status });
    return await response.json();
  }

  static async updateSchema(jsonPayload: string): Promise<string> {
    return this.postToEndpoint('/v1/schema/', jsonPayload);
  }

  static async registerAssay(jsonPayload: string): Promise<string> {
    return this.postToEndpoint('/v1/assays', jsonPayload);
  }

  static async retrieveEntity(scope: string): Promise<any> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(
      this.container.id,
      `/v1/${scope}/`,
      { method: 'GET'},
    );
    if (!response.ok)
      throw new Error(`HTTP error!: ${response.status}`, { cause: response.status });
    return response.json();
  }

  static async getCompoundById(id: number): Promise<any> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(
      this.container.id,
      `/v1/compounds/id/${id}`,
      { method: 'GET' },
    );
    return response.json();
  }

  static async getCompoundByCorporateId(corporateCompoundId: string): Promise<any> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(
      this.container.id,
      `/v1/compounds/${corporateCompoundId}`,
      { method: 'GET' },
    );
    return response.json();
  }

  static async getBatchByCorporateId(corporateBatchId: string): Promise<any> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(
      this.container.id,
      `/v1/batches/${corporateBatchId}`,
      { method: 'GET' },
    );
    return response.json();
  }

  static async checkCompoundExists(smiles: string): Promise<boolean> {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(
      this.container.id,
      `/v1/compounds/check-duplicate?smiles=${smiles}`,
      {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
      },
    );
    if (response.ok) {
      const data = await response.json();
      return data['exists'];
    }
    return false;
  }

  static async registerBulk(
    file: DG.FileInfo,
    scope: string,
    mapping: string,
    errorHandling: string,
  ): Promise<DG.DataFrame> {
    try {
      const formData = await buildFormData(file, mapping, errorHandling);
      const response = await grok.dapi.docker.dockerContainers.fetchProxy(
        this.container.id,
        scopeToUrl[scope],
        { method: 'POST', body: formData },
      );

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`MolTrack API error: ${response.status} ${errorText}`);
      }

      const json = await response.json();
      return DG.DataFrame.fromJson(JSON.stringify(json));
    } catch (e) {
      grok.shell.error(String(e));
      throw e;
    }
  }
}

async function buildFormData(
  file: DG.FileInfo,
  mapping: string,
  errorHandling: string,
): Promise<FormData> {
  const formData = new FormData();

  const rawBytes = await file.readAsBytes();
  const bytes = new Uint8Array(rawBytes.length);
  bytes.set(rawBytes);
  const blob = new Blob([bytes], { type: 'text/csv' });
  const preparedFile = new File([blob], file.fileName, { type: 'text/csv' });

  formData.append('file', preparedFile, file.name);
  formData.append('error_handling', errorHandling || ErrorHandling.REJECT_ROW);
  formData.append('mapping', mapping || '');
  formData.append('output_format', ResultOutput.JSON);

  return formData;
}
