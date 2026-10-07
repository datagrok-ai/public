/// REST clients for the observability routes (`/alerts`, `/problems`, `/errors`, `/logging`, `/log/timeline`).
import {NodeApiClient, buildQuery} from './node-dapi';

export type Query = Record<string, string | number | boolean | undefined>;

const seg = (id: string) => encodeURIComponent(String(id));

export class NodeAlertsClient {
  constructor(private client: NodeApiClient) {}

  list(q: Query = {}): Promise<any[]> { return this.client.get(`/alerts${buildQuery(q)}`); }
  detection(): Promise<any> { return this.client.get('/alerts/detection'); }
  get(id: string): Promise<any> { return this.client.get(`/alerts/${seg(id)}`); }

  transition(id: string, action: 'ack' | 'mute' | 'unmute' | 'resolve', body: Record<string, any>): Promise<any> {
    return this.client.post(`/alerts/${seg(id)}/${action}`, body);
  }

  problems(q: Query = {}): Promise<any[]> { return this.client.get(`/problems${buildQuery(q)}`); }
  problem(id: string): Promise<any> { return this.client.get(`/problems/${seg(id)}`); }
  history(id: string, q: Query = {}): Promise<any[]> { return this.client.get(`/problems/${seg(id)}/history${buildQuery(q)}`); }
  setStatus(id: string, body: Record<string, any>): Promise<any> { return this.client.post(`/problems/${seg(id)}/status`, body); }
}

export class NodeErrorsClient {
  constructor(private client: NodeApiClient) {}

  /** Occurrences without `by`, aggregate rows with it; CSV text or a JSON array when `format` is set. */
  query(q: Query): Promise<any> { return this.client.get(`/errors${buildQuery(q)}`); }
  diff(q: Query): Promise<any> { return this.client.get(`/errors/diff${buildQuery(q)}`); }
  show(signature: string, q: Query = {}): Promise<any> { return this.client.get(`/errors/${seg(signature)}${buildQuery(q)}`); }
}

/** Logging policy, overrides and capture rules (LoggingRouter), plus the timeline (ActionLoggerRouter). */
export class NodeLoggingClient {
  constructor(private client: NodeApiClient) {}

  policy(q: Query = {}): Promise<any> { return this.client.get(`/logging/policy${buildQuery(q)}`); }
  setPolicy(body: {set: Record<string, any>; reason?: string}): Promise<any> {
    return this.client.request('PUT', '/logging/policy', body);
  }
  effective(q: Query): Promise<any> { return this.client.get(`/logging/policy/effective${buildQuery(q)}`); }
  history(q: Query = {}): Promise<any[]> { return this.client.get(`/logging/policy/history${buildQuery(q)}`); }
  revert(body: Record<string, any>): Promise<any> { return this.client.post('/logging/policy/revert', body); }
  overrides(): Promise<any[]> { return this.client.get('/logging/policy/overrides'); }
  addOverride(body: Record<string, any>): Promise<any> { return this.client.post('/logging/policy/overrides', body); }

  captureRules(q: Query = {}): Promise<any[]> { return this.client.get(`/logging/capture${buildQuery(q)}`); }
  captureRule(id: string): Promise<any> { return this.client.get(`/logging/capture/${seg(id)}`); }
  addCaptureRule(body: Record<string, any>): Promise<any> { return this.client.post('/logging/capture', body); }
  stopCaptureRule(id: string, reason?: string): Promise<any> {
    return this.client.request('DELETE', `/logging/capture/${seg(id)}`, {reason});
  }

  timeline(q: Query): Promise<any[]> { return this.client.get(`/log/timeline${buildQuery(q)}`); }
}
