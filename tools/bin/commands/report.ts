import fs from 'fs';
import os from 'os';
import path from 'path';
import * as color from '../utils/color-utils';
import {getDevKey, getToken} from '../utils/test-utils';

const fetch = require('node-fetch');

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
  case 'resolve':
    return await handleResolve(args);
  case 'ticket':
    return await handleTicket(args);
  default:
    return false;
  }
}

async function handleFetch(args: ReportArgs): Promise<boolean> {
  const instance = args._[2];
  const number = args._[3];
  if (!instance || !number) {
    color.error('Usage: grok report fetch <instance> <number>');
    return false;
  }

  try {
    const {url, key} = getDevKey(instance);
    const token = await getToken(url, key);

    console.log(`Searching for report #${number}...`);
    const searchResp = await fetch(
      `${url}/reports?text=number%3D${encodeURIComponent(number)}`,
      {headers: {Authorization: token}},
    );
    if (!searchResp.ok) {
      color.error(`Report search failed (HTTP ${searchResp.status})`);
      return false;
    }
    const results = await searchResp.json();
    if (!Array.isArray(results) || results.length === 0) {
      color.error(`Report #${number} not found`);
      return false;
    }

    const reportData = results[0];
    const reportId = reportData.id || reportData.Id;
    if (!reportId) {
      color.error('Report found but has no id field');
      return false;
    }

    console.log(`Downloading report ${reportId}...`);
    const downloadResp = await fetch(`${url}/reports/${reportId}/zip`, {
      headers: {Authorization: token},
    });
    if (!downloadResp.ok) {
      color.error(`Report download failed (HTTP ${downloadResp.status})`);
      return false;
    }

    const buffer = await downloadResp.buffer();
    const outputPath = path.join(os.tmpdir(), `report_${instance}_${number}.zip`);
    fs.writeFileSync(outputPath, buffer);

    const metaPath = outputPath.replace('.zip', '_meta.json');
    fs.writeFileSync(metaPath, JSON.stringify(reportData, null, 2));

    color.success(`Report saved to: ${outputPath}`);
    console.log(outputPath);
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
