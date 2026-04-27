import fs from 'fs';
import os from 'os';
import path from 'path';
import * as color from '../utils/color-utils';
import {getDevKey, getToken} from '../utils/test-utils';

const fetch = require('node-fetch');
const AdmZip = require('adm-zip');

interface ReportArgs {
  _: string[];
  help?: boolean;
  [key: string]: any;
}

export async function report(args: ReportArgs): Promise<boolean> {
  const subcommand = args._[1];

  switch (subcommand) {
  case 'fetch':
    return await handleFetch(args);
  case 'read':
    return await handleRead(args);
  case 'resolve':
    return await handleResolve(args);
  case 'ticket':
    return await handleTicket(args);
  default:
    return false;
  }
}

interface ReportBundle {
  zipBuffer: Buffer;
  metaJson: any | null;
}

async function fetchReportMeta(url: string, token: string, number: string): Promise<any | null> {
  const resp = await fetch(
    `${url}/reports?text=number%3D${encodeURIComponent(number)}`,
    {headers: {Authorization: token}},
  );
  if (!resp.ok)
    return null;
  const arr = await resp.json();
  if (!Array.isArray(arr) || arr.length === 0)
    return null;
  return arr[0];
}

async function fetchReportBundle(instance: string, number: string): Promise<ReportBundle> {
  const {url, key} = getDevKey(instance);
  const token = await getToken(url, key);

  const byNumberResp = await fetch(
    `${url}/reports/by-number/${encodeURIComponent(number)}/zip`,
    {headers: {Authorization: token}},
  );
  if (byNumberResp.ok) {
    const zipBuffer = await byNumberResp.buffer();
    const metaJson = await fetchReportMeta(url, token, number);
    return {zipBuffer, metaJson};
  }
  if (byNumberResp.status !== 404 && byNumberResp.status !== 405)
    throw new Error(`HTTP ${byNumberResp.status} from /reports/by-number/${number}/zip`);

  const metaJson = await fetchReportMeta(url, token, number);
  if (metaJson == null)
    throw new Error(`Report #${number} not found`);
  const reportId = metaJson.id || metaJson.Id;
  if (!reportId)
    throw new Error('Report found but has no id field');
  const downloadResp = await fetch(`${url}/reports/${reportId}/zip`, {
    headers: {Authorization: token},
  });
  if (!downloadResp.ok)
    throw new Error(`Report download failed (HTTP ${downloadResp.status})`);
  const zipBuffer = await downloadResp.buffer();
  return {zipBuffer, metaJson};
}

async function handleFetch(args: ReportArgs): Promise<boolean> {
  const instance = args._[2];
  const number = args._[3];
  if (!instance || !number) {
    color.error('Usage: grok report fetch <instance> <number>');
    return false;
  }

  try {
    console.log(`Fetching report #${number}...`);
    const {zipBuffer, metaJson} = await fetchReportBundle(instance, number);
    const outputPath = path.join(os.tmpdir(), `report_${instance}_${number}.zip`);
    fs.writeFileSync(outputPath, zipBuffer);

    if (metaJson != null) {
      const metaPath = outputPath.replace('.zip', '_meta.json');
      fs.writeFileSync(metaPath, JSON.stringify(metaJson, null, 2));
    }

    color.success(`Report saved to: ${outputPath}`);
    console.log(outputPath);
    return true;
  }
  catch (err: any) {
    color.error(`Error: ${err.message}`);
    return false;
  }
}

interface LoadedReport {
  data: any;
  zip: any | null;
  d42Names: string[];
}

function looksLikePath(s: string): boolean {
  return s.includes('/') || s.includes('\\') || /\.(zip|json)$/i.test(s) || fs.existsSync(s);
}

function loadFromZip(z: any): LoadedReport {
  const entries = z.getEntries();
  const reportEntry = entries.find((e: any) => e.entryName.toLowerCase().endsWith('report.json'));
  if (reportEntry == null) {
    const names = entries.map((e: any) => e.entryName).join(', ');
    throw new Error(`report.json not found in zip. Contents: ${names}`);
  }
  const data = JSON.parse(reportEntry.getData().toString('utf-8'));
  const d42Names = entries
    .map((e: any) => e.entryName)
    .filter((n: string) => n.toLowerCase().endsWith('.d42'));
  return {data, zip: z, d42Names};
}

function loadFromBuffer(buf: Buffer): LoadedReport {
  try {
    const z = new AdmZip(buf);
    return loadFromZip(z);
  }
  catch {
    const text = buf.toString('utf-8');
    return {data: JSON.parse(text), zip: null, d42Names: []};
  }
}

function loadFromPath(p: string): LoadedReport {
  if (p.toLowerCase().endsWith('.json')) {
    const text = fs.readFileSync(p, 'utf-8');
    return {data: JSON.parse(text), zip: null, d42Names: []};
  }
  return loadFromBuffer(fs.readFileSync(p));
}

function unwrapEnvelope(data: any): { meta: any, body: any } {
  if (data && typeof data === 'object' && '#type' in data && 'data' in data
      && typeof data.data === 'object' && data.data != null) {
    const meta: any = {};
    if (data.id != null) meta.id = data.id;
    if (data.createdOn != null) meta.createdOn = data.createdOn;
    if (data['#type'] != null) meta['#type'] = data['#type'];
    return {meta, body: data.data};
  }
  return {meta: {}, body: data};
}

function loadSidecar(inputPath: string): any | null {
  const stem = inputPath.replace(/\.[^./\\]+$/, '');
  const sidecar = `${stem}_meta.json`;
  if (sidecar !== inputPath && fs.existsSync(sidecar))
    return JSON.parse(fs.readFileSync(sidecar, 'utf-8'));
  return null;
}

function ensureParentDir(p: string): void {
  const dir = path.dirname(p);
  if (dir && dir !== '.' && !fs.existsSync(dir))
    fs.mkdirSync(dir, {recursive: true});
}

function extractScreenshot(zip: any, body: any, outPath: string): boolean {
  if (zip != null) {
    const entries = zip.getEntries();
    const e = entries.find((x: any) => x.entryName.toLowerCase().endsWith('screenshot.png'));
    if (e != null) {
      ensureParentDir(outPath);
      fs.writeFileSync(outPath, e.getData());
      return true;
    }
  }
  const b64 = body && body.screenshot;
  if (typeof b64 === 'string' && b64.length > 0) {
    const stripped = b64.includes(',') ? b64.slice(b64.indexOf(',') + 1) : b64;
    try {
      const buf = Buffer.from(stripped, 'base64');
      ensureParentDir(outPath);
      fs.writeFileSync(outPath, buf);
      return true;
    }
    catch {
      return false;
    }
  }
  return false;
}

function extractD42(zip: any, names: string[], outDir: string): string[] {
  if (zip == null || names.length === 0) return [];
  fs.mkdirSync(outDir, {recursive: true});
  const written: string[] = [];
  for (var name of names) {
    const entry = zip.getEntry(name);
    if (entry == null) continue;
    const out = path.join(outDir, path.basename(name));
    fs.writeFileSync(out, entry.getData());
    written.push(out);
  }
  return written;
}

async function handleRead(args: ReportArgs): Promise<boolean> {
  const positional = args._.slice(2);
  if (positional.length === 0) {
    color.error('Usage: grok report read <path> | <instance> <number>');
    return false;
  }

  const screenshotOut = args['extract-screenshot'] as string | undefined;
  const d42Dir = args['extract-d42'] as string | undefined;
  const extractActions = args['extract-actions'] === true;

  try {
    let loaded: LoadedReport;
    let sidecarMeta: any | null = null;
    let inputPath: string | null = null;
    let networkBase: string | null = null;

    if (positional.length >= 2 && !looksLikePath(positional[0])) {
      const [instance, number] = positional;
      const bundle = await fetchReportBundle(instance, number);
      loaded = loadFromBuffer(bundle.zipBuffer);
      sidecarMeta = bundle.metaJson;
      networkBase = path.join(os.tmpdir(), `report_${instance}_${number}`);
    }
    else {
      inputPath = path.resolve(positional[0]);
      if (!fs.existsSync(inputPath)) {
        color.error(`File not found: ${inputPath}`);
        return false;
      }
      loaded = loadFromPath(inputPath);
      sidecarMeta = loadSidecar(inputPath);
    }

    const {meta: envelopeMeta, body} = unwrapEnvelope(loaded.data);
    const meta = Object.assign({}, envelopeMeta, sidecarMeta || {});

    const output: any = {meta, ...body};
    const files: any = {};

    if (screenshotOut) {
      if (extractScreenshot(loaded.zip, body, screenshotOut)) {
        files.screenshot = screenshotOut;
        output.screenshot = screenshotOut;
      }
    }

    if (d42Dir) {
      const written = extractD42(loaded.zip, loaded.d42Names, d42Dir);
      if (written.length > 0)
        files.d42 = written;
    }

    if (extractActions && Array.isArray(body && body.actions)) {
      const stem = inputPath != null
        ? inputPath.replace(/\.[^./\\]+$/, '')
        : networkBase;
      if (stem != null) {
        const actionsPath = `${stem}_actions.json`;
        fs.writeFileSync(actionsPath, JSON.stringify(body.actions, null, 2));
        files.actions = actionsPath;
      }
    }

    if (Object.keys(files).length > 0)
      output.files = files;

    process.stdout.write(JSON.stringify(output));
    process.stdout.write('\n');
    return true;
  }
  catch (err: any) {
    color.error(`Error: ${err.message}`);
    return false;
  }
}

async function handleResolve(args: ReportArgs): Promise<boolean> {
  const instance = args._[2];
  const number = args._[3];
  if (!instance || !number) {
    color.error('Usage: grok report resolve <instance> <number>');
    return false;
  }

  try {
    const metaPath = path.join(os.tmpdir(), `report_${instance}_${number}_meta.json`);
    if (!fs.existsSync(metaPath)) {
      color.error(`Meta file not found: ${metaPath}`);
      color.warn('Hint: was the report fetched via `grok report fetch` first?');
      return false;
    }

    const meta = JSON.parse(fs.readFileSync(metaPath, 'utf-8'));
    const reportId = meta.id || meta.Id;
    if (!reportId) {
      color.error(`No report id in meta file: ${metaPath}`);
      return false;
    }

    const {url, key} = getDevKey(instance);
    const token = await getToken(url, key);

    console.log(`Resolving report #${number} (id: ${reportId})...`);
    const resp = await fetch(`${url}/reports/${reportId}/resolve`, {
      method: 'POST',
      headers: {Authorization: token},
    });

    if (resp.ok) {
      color.success(`Report #${number} resolved on ${instance}`);
      return true;
    }

    const body = await resp.text();
    color.error(`Resolve failed (HTTP ${resp.status}): ${body}`);
    return false;
  }
  catch (err: any) {
    color.error(`Error: ${err.message}`);
    return false;
  }
}

async function handleTicket(args: ReportArgs): Promise<boolean> {
  const instance = args._[2];
  const reportId = args._[3];
  if (!instance || !reportId) {
    color.error('Usage: grok report ticket <instance> <report-id>');
    return false;
  }

  try {
    const {url, key} = getDevKey(instance);
    const token = await getToken(url, key);

    console.log('Getting current user...');
    const userResp = await fetch(`${url}/users/current`, {
      headers: {Authorization: token},
    });
    if (!userResp.ok) {
      color.error(`Failed to get current user (HTTP ${userResp.status})`);
      return false;
    }
    const user = await userResp.json();
    const userId = user.id || user.Id;
    if (!userId) {
      color.error('No user id in response');
      return false;
    }

    console.log(`Creating JIRA ticket for report ${reportId}...`);
    const ticketResp = await fetch(
      `${url}/reports/${reportId}/jira?assigneeId=${userId}`,
      {
        method: 'POST',
        headers: {Authorization: token, 'Content-Type': 'application/json'},
      },
    );

    if (ticketResp.status !== 200 && ticketResp.status !== 201) {
      const body = await ticketResp.text();
      color.error(`JIRA ticket creation failed (HTTP ${ticketResp.status}): ${body}`);
      return false;
    }

    const result = await ticketResp.json();
    const ticketKey = result.key;
    if (!ticketKey) {
      color.error(`No ticket key in response: ${JSON.stringify(result).slice(0, 200)}`);
      return false;
    }

    color.success(`Created ticket: ${ticketKey}`);
    console.log(ticketKey);
    return true;
  }
  catch (err: any) {
    color.error(`Error: ${err.message}`);
    return false;
  }
}
