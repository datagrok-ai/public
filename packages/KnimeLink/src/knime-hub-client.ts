import * as grok from 'datagrok-api/grok';
import {IKnimeClient} from './knime-client';
import {KnimeInputParam, KnimeWorkflowSpec, KnimeParamType, KnimeExecutionInput, KnimeExecutionResult, KnimeJobStatus, KnimeDeployment, KnimeOutputResource} from './types';
import {KnimeJobState} from './constants';
import {getCredentials, getBasicAuthHeader} from './credentials';

/**
 * Client for KNIME Business Hub.
 *
 * API lives on a separate subdomain: if the Hub is at `hub.knime.com`,
 * the REST API is at `api.hub.knime.com`.  The constructor accepts
 * either form and normalises it to the `api.` variant.
 *
 * Key endpoints:
 *   GET  /accounts/identity                     — current user (includes teams[])
 *   GET  /deployments/{scope}                   — list all deployments for a team
 *   GET  /deployments/{id}                      — single deployment details
 *   POST /deployments/{id}/execution            — synchronous execution (id = rest:{uuid})
 *   POST /deployments/{id}/jobs                 — create async job
 *   GET  /jobs/{uuid}                           — poll job status
 *   DELETE /jobs/{uuid}                         — discard job
 *   DELETE /jobs/{uuid}/execution               — cancel running job
 *
 * Auth: HTTP Basic with an Application Password (id + secret).
 */
export class KnimeHubClient implements IKnimeClient {
  private apiUrl: string;
  private teamScopes: string[] | null = null;

  constructor(baseUrl: string) {
    this.apiUrl = toApiUrl(baseUrl);
  }

  // ------------------------------------------------------------------ HTTP
  private async request<T>(method: string, path: string, body?: any): Promise<T> {
    const creds = await getCredentials();
    const headers: Record<string, string> = {
      'Authorization': getBasicAuthHeader(creds.id, creds.secret),
      'Accept': 'application/json',
    };
    if (body !== undefined)
      headers['Content-Type'] = 'application/json';

    const url = `${this.apiUrl}${path}`;
    const response = await grok.dapi.fetchProxy(url, {
      method,
      headers,
      body: body !== undefined ? JSON.stringify(body) : undefined,
    });

    if (!response.ok) {
      const text = await response.text();
      throw new Error(`KNIME Hub API error ${response.status}: ${text}`);
    }

    const ct = response.headers.get('content-type') ?? '';
    if (ct.includes('application/json'))
      return await response.json() as T;
    return await response.text() as unknown as T;
  }

  private async requestMultipart<T>(method: string, path: string, formData: FormData): Promise<T> {
    const creds = await getCredentials();
    const url = `${this.apiUrl}${path}`;
    const response = await grok.dapi.fetchProxy(url, {
      method,
      headers: {
        'Authorization': getBasicAuthHeader(creds.id, creds.secret),
        'Accept': 'application/json',
      },
      body: formData,
    });

    if (!response.ok) {
      const text = await response.text();
      throw new Error(`KNIME Hub API error ${response.status}: ${text}`);
    }

    const ct = response.headers.get('content-type') ?? '';
    if (ct.includes('application/json'))
      return await response.json() as T;
    return await response.text() as unknown as T;
  }

  // ----------------------------------------------------------- deployments
  async listDeployments(): Promise<KnimeDeployment[]> {
    const scopes = await this.getTeamScopes();
    const list: KnimeDeployment[] = [];
    for (const scope of scopes) {
      try {
        const data = await this.request<any>('GET', `/deployments/${scope}`);
        const items = Array.isArray(data) ? data : data?.deployments ?? data?.items ?? [];
        for (const d of items)
          list.push({
            id: d.id,
            name: d.name ?? d.workflowPath ?? d.id,
            type: d.type ?? 'REST',
            state: d.state,
            workflowPath: d.workflowPath,
          });
      }
      catch {
        // scope may not have deployments or user may lack permission
      }
    }
    return list;
  }

  private async getTeamScopes(): Promise<string[]> {
    if (this.teamScopes)
      return this.teamScopes;

    const identity = await this.request<any>('GET', '/accounts/identity');
    const scopes: string[] = [];

    // UserAccount includes teams[] with BaseTeam objects, each having id like "account:team:{uuid}"
    if (Array.isArray(identity?.teams)) {
      for (const t of identity.teams) {
        if (t.id)
          scopes.push(encodeURIComponent(t.id));
      }
    }

    // Fallback: use the user's own account id
    if (scopes.length === 0 && identity?.id)
      scopes.push(encodeURIComponent(identity.id));

    if (scopes.length === 0)
      throw new Error('Could not determine KNIME Hub team from /accounts/identity');
    this.teamScopes = scopes;
    return this.teamScopes;
  }

  // ------------------------------------------------------ workflow inputs
  async getWorkflowInputs(id: string): Promise<KnimeWorkflowSpec> {
    try {
      const spec = await this.request<any>('GET', `/deployments/${id}/open-api`);
      return {
        inputs: this.parseOpenApiInputs(spec),
        hasFileOutputs: this.detectFileOutputs(spec),
      };
    }
    catch {
      return {inputs: [], hasFileOutputs: false};
    }
  }

  private parseOpenApiInputs(spec: any): KnimeInputParam[] {
    const params: KnimeInputParam[] = [];
    const inputSchema = spec?.components?.schemas?.InputParameters;

    // File inputs are in the multipart/form-data schema of the /execution endpoint, not in InputParameters
    this.parseMultipartFileInputs(spec, params);

    if (!inputSchema?.properties)
      return params;

    const required = new Set(inputSchema.required ?? []);
    for (const [name, prop] of Object.entries(inputSchema.properties) as [string, any][]) {
      // Table input: object with table-spec/table-data in example or properties
      if (this.isTableInput(prop)) {
        const tableSpec = this.extractTableSpec(prop);
        params.push({
          name, type: 'table', required: required.has(name),
          description: prop.description, tableSpec,
        });
        continue;
      }

      // File input
      if (prop?.format === 'binary' || (prop?.type === 'string' && prop?.format === 'base64')) {
        params.push({
          name, type: 'file', required: required.has(name), description: prop.description,
        });
        continue;
      }

      // Object with typed sub-properties: expand into individual inputs grouped under parent name
      if (prop?.type === 'object' && prop?.properties && Object.keys(prop.properties).length > 0) {
        const subRequired = new Set(prop.required ?? []);
        for (const [subName, subProp] of Object.entries(prop.properties) as [string, any][]) {
          params.push({
            name: subName,
            type: this.mapScalarType(subProp),
            required: subRequired.has(subName),
            description: subProp.description,
            defaultValue: subProp.default,
            group: name,
            groupDescription: prop.description || undefined,
          });
        }
        continue;
      }

      // Object with empty properties (arbitrary variable input): JSON text area
      if (prop?.type === 'object') {
        params.push({
          name, type: 'json', required: required.has(name),
          description: prop.description, defaultValue: prop.example,
        });
        continue;
      }

      // Scalar at top level (unlikely but handle)
      params.push({
        name, type: this.mapScalarType(prop), required: required.has(name),
        description: prop.description, defaultValue: prop.default,
      });
    }
    return params;
  }

  private isTableInput(prop: any): boolean {
    if (prop?.type !== 'object')
      return false;
    // Check properties for table-spec/table-data
    if (prop.properties?.['table-spec'] && prop.properties?.['table-data'])
      return true;
    // Check example for table-spec/table-data
    if (prop.example?.['table-spec'] && prop.example?.['table-data'])
      return true;
    return false;
  }

  private extractTableSpec(prop: any): {name: string; type: string}[] | undefined {
    const spec = prop.example?.['table-spec'];
    if (!Array.isArray(spec))
      return undefined;
    // Each entry is {columnName: type}
    return spec.map((s: any) => {
      const colName = Object.keys(s)[0];
      return {name: colName, type: s[colName]};
    });
  }

  private detectFileOutputs(spec: any): boolean {
    const paths = spec?.paths;
    if (!paths)
      return false;
    for (const [, pathVal] of Object.entries(paths) as [string, any][]) {
      if (pathVal?.post?.['x-knime-output-resources']?.schema?.properties)
        return true;
    }
    return false;
  }

  private parseMultipartFileInputs(spec: any, params: KnimeInputParam[]): void {
    const paths = spec?.paths;
    if (!paths)
      return;
    for (const [pathKey, pathVal] of Object.entries(paths) as [string, any][]) {
      if (!pathKey.endsWith('/execution'))
        continue;
      const multipart = pathVal?.post?.requestBody?.content?.['multipart/form-data'];
      if (!multipart?.schema?.properties)
        continue;
      for (const [name, prop] of Object.entries(multipart.schema.properties) as [string, any][]) {
        if (prop?.format === 'binary' || prop?.type === 'string')
          params.push({name, type: 'file', required: false, description: prop.description});
      }
      break;
    }
  }

  private mapScalarType(prop: any): KnimeParamType {
    if (prop?.type === 'integer' || prop?.format === 'int32' || prop?.format === 'int64')
      return 'int';
    if (prop?.type === 'number' || prop?.format === 'double' || prop?.format === 'float')
      return 'double';
    if (prop?.type === 'boolean')
      return 'boolean';
    return 'string';
  }

  // ------------------------------------------------------------ execution
  async executeSyncWorkflow(id: string, input: KnimeExecutionInput): Promise<KnimeExecutionResult> {
    // Check if any input values are File/Blob — if so, use multipart/form-data
    const fileKeys: string[] = [];
    for (const key of Object.keys(input)) {
      if (input[key] instanceof File || input[key] instanceof Blob)
        fileKeys.push(key);
    }

    let data: any;
    if (fileKeys.length > 0) {
      const formData = new FormData();
      const jsonInputs: {[key: string]: any} = {};
      for (const key of Object.keys(input)) {
        if (fileKeys.includes(key))
          formData.append(key, input[key] as Blob);
        else
          jsonInputs[key] = input[key];
      }
      for (const [key, value] of Object.entries(jsonInputs))
        formData.append(key, new Blob([JSON.stringify(value)], {type: 'application/json'}));
      data = await this.requestMultipart<any>('POST', `/deployments/${id}/execution?reset=true`, formData);
    }
    else
      data = await this.request<any>('POST', `/deployments/${id}/execution?reset=true`, input);

    // Response may come as a string (when content-type is not application/json) — parse it
    if (typeof data === 'string') {
      try { data = JSON.parse(data); }
      catch { /* leave as-is */ }
    }

    this.checkForExecutionErrors(data);
    const result = (data?.outputValues || data?.outputResources || data?.['@controls'])
      ? this.parseExecutionResult(data)
      : this.parseFlatResult(data);

    // Sync execution response includes a job ID that can be used to fetch output resources
    if (data?.id)
      result.jobId = data.id;

    return result;
  }

  async startAsyncJob(id: string, input: KnimeExecutionInput): Promise<string> {
    // Step 1: Create the job
    const createData = await this.request<any>('POST', `/deployments/${id}/jobs`);
    const jobId = createData.id ?? createData.jobId;

    // Step 2: Execute with input asynchronously
    const hasInput = Object.keys(input).length > 0;
    if (!hasInput) {
      await this.request<any>('POST', `/jobs/${jobId}?async=true`);
      return jobId;
    }

    // Check for file (Blob/File) inputs — if present, use multipart/form-data
    const fileKeys: string[] = [];
    for (const key of Object.keys(input)) {
      if (input[key] instanceof File || input[key] instanceof Blob)
        fileKeys.push(key);
    }

    if (fileKeys.length > 0) {
      const formData = new FormData();
      const jsonInputs: {[key: string]: any} = {};
      for (const key of Object.keys(input)) {
        if (fileKeys.includes(key))
          formData.append(key, input[key] as Blob);
        else
          jsonInputs[key] = input[key];
      }
      for (const [key, value] of Object.entries(jsonInputs))
        formData.append(key, new Blob([JSON.stringify(value)], {type: 'application/json'}));
      await this.requestMultipart<any>('POST', `/jobs/${jobId}?async=true`, formData);
    }
    else
      await this.request<any>('POST', `/jobs/${jobId}?async=true`, input);

    return jobId;
  }

  async getJobStatus(jobId: string): Promise<KnimeJobStatus> {
    const data = await this.request<any>('GET', `/jobs/${jobId}`);
    const state = data.workflowState ?? data.state ?? KnimeJobState.Undefined;
    return {
      id: jobId,
      state: state as KnimeJobState,
      message: data.message ?? this.extractNodeMessages(data) ?? undefined,
      _rawData: data,
    };
  }

  async getJobResult(jobId: string): Promise<KnimeExecutionResult> {
    const data = await this.request<any>('GET', `/jobs/${jobId}`);
    return this.buildJobResult(jobId, data);
  }

  /** Parse job data and immediately fetch output resources before the job gets discarded. */
  async buildJobResult(jobId: string, data: any): Promise<KnimeExecutionResult> {
    const result = this.parseExecutionResult(data);
    result.jobId = jobId;

    if (result.outputResources) {
      const fetchedResources: {name: string; resource: KnimeOutputResource}[] = [];
      await Promise.all(Object.keys(result.outputResources).map(async (resourceId) => {
        try {
          const resource = await this.fetchOutputResource(jobId, resourceId);
          fetchedResources.push({name: result.outputResources![resourceId] || resourceId, resource});
        }
        catch (e: any) {
          console.warn(`Failed to fetch output resource '${resourceId}':`, e);
        }
      }));
      if (fetchedResources.length > 0)
        result.fetchedResources = fetchedResources;
    }

    return result;
  }

  async fetchOutputResource(jobId: string, resourceId: string): Promise<KnimeOutputResource> {
    const creds = await getCredentials();
    const url = `${this.apiUrl}/jobs/${jobId}/output-resources/${resourceId}`;

    // Step 1: Request with redirect: manual to get the S3 pre-signed URL
    // (KNIME redirects to S3 which doesn't have CORS headers for our origin)
    const redirectResponse = await grok.dapi.fetchProxy(url, {
      method: 'GET',
      headers: {
        'Authorization': getBasicAuthHeader(creds.id, creds.secret),
      },
      redirect: 'manual',
    });

    // If we got a redirect, fetch the target URL through the proxy
    if (redirectResponse.type === 'opaqueredirect' || redirectResponse.status === 301 || redirectResponse.status === 302 || redirectResponse.status === 307) {
      const location = redirectResponse.headers.get('location');
      if (location) {
        const response = await grok.dapi.fetchProxy(location, {method: 'GET'});
        if (!response.ok) {
          const text = await response.text();
          throw new Error(`KNIME resource download error ${response.status}: ${text}`);
        }
        return this.parseResourceResponse(response);
      }
    }

    // No redirect — parse the response directly
    if (!redirectResponse.ok) {
      const text = await redirectResponse.text();
      throw new Error(`KNIME Hub API error ${redirectResponse.status}: ${text}`);
    }
    return this.parseResourceResponse(redirectResponse);
  }

  private async parseResourceResponse(response: Response): Promise<KnimeOutputResource> {
    const ct = response.headers.get('content-type') ?? 'application/octet-stream';
    if (ct.includes('application/json'))
      return {contentType: ct, json: await response.json()};
    if (ct.startsWith('text/'))
      return {contentType: ct, text: await response.text()};
    return {contentType: ct, blob: await response.blob()};
  }

  /** Build the viewable data-app URL from deployment info and job ID. */
  getDataAppUrl(deploymentName: string, deploymentId: string, jobId: string): string {
    const hubHost = this.apiUrl.replace('https://api.', '').replace('http://api.', '');
    const encodedName = encodeURIComponent(deploymentName);
    return `https://apps.${hubHost}/d/${encodedName}~${deploymentId}/${jobId}`;
  }

  async cancelJob(jobId: string): Promise<void> {
    await this.request<any>('DELETE', `/jobs/${jobId}/execution`);
  }

  // --------------------------------------------------------------- helpers
  private checkForExecutionErrors(data: any): void {
    const state = data?.workflowState ?? data?.state;
    if (state === KnimeJobState.ExecutionFailed || state === KnimeJobState.ExecutionFailedWithContent ||
      state === KnimeJobState.NotExecutable) {
      const messages = this.extractNodeMessages(data);
      throw new Error(messages || `Workflow execution failed (${state})`);
    }
  }

  private extractNodeMessages(data: any): string | null {
    if (Array.isArray(data?.nodeMessages) && data.nodeMessages.length > 0)
      return data.nodeMessages.map((m: any) => `${m.node}: ${m.message}`).join('\n');
    if (Array.isArray(data?.errors) && data.errors.length > 0)
      return data.errors.join('; ');
    return null;
  }

  private parseFlatResult(data: any): KnimeExecutionResult {
    const outputs: {[key: string]: any} = {};
    const outputTables: any[] = [];
    if (!data)
      return {outputs};

    // Top-level table-spec/table-data: the response itself is a single table
    if (data['table-spec'] && data['table-data']) {
      outputTables.push(data);
      return {outputs, outputTables};
    }

    for (const key of Object.keys(data)) {
      const val = data[key];
      if (Array.isArray(val) || (typeof val === 'object' && val !== null && val['table-data']))
        outputTables.push(val);
      else
        outputs[key] = val;
    }
    return {outputs, outputTables: outputTables.length > 0 ? outputTables : undefined};
  }

  private parseExecutionResult(data: any): KnimeExecutionResult {
    const outputs: {[key: string]: any} = {};
    const outputTables: any[] = [];

    // Extract outputValues (inline key-value results from REST deployments)
    if (data?.outputValues)
      for (const key of Object.keys(data.outputValues)) {
        const val = data.outputValues[key];
        if (Array.isArray(val) || (typeof val === 'object' && val !== null && val['table-data']))
          outputTables.push(val);
        else if (val !== null)
          outputs[key] = val;
      }

    // Extract outputResources (named resources fetched separately)
    const outputResources: {[id: string]: string} | undefined =
      data?.outputResources && Object.keys(data.outputResources).length > 0
        ? data.outputResources : undefined;

    // Extract session creation URL from @controls (for data-app deployments)
    let dataAppUrl: string | undefined;
    const controls = data?.['@controls'];
    if (controls?.['knime:new-session']?.href)
      dataAppUrl = controls['knime:new-session'].href;
    else if (controls?.['knime:new-view-session']?.href)
      dataAppUrl = controls['knime:new-view-session'].href;

    return {
      outputs,
      outputTables: outputTables.length > 0 ? outputTables : undefined,
      outputResources,
      dataAppUrl,
    };
  }
}

/**
 * Derive the API URL from the user-entered base URL.
 *
 * Users may enter either `https://hub.knime.com` or
 * `https://api.hub.knime.com` — we normalise to the `api.` form.
 * Trailing slashes are stripped.
 */
function toApiUrl(baseUrl: string): string {
  let url = baseUrl.replace(/\/+$/, '');
  try {
    const parsed = new URL(url);
    if (!parsed.hostname.startsWith('api.'))
      parsed.hostname = 'api.' + parsed.hostname;
    url = parsed.origin;
  }
  catch {
    // If URL parsing fails, just prefix naively
    if (!url.includes('://api.'))
      url = url.replace('://', '://api.');
  }
  return url;
}
