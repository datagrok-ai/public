import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {imageRenderer} from '@datagrok-libraries/db-explorer/src/renderer';
import {matchAndParseQuery} from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';
import {powerSearchQueryTable} from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';

export async function registerDGUserHandler(_package: DG.Package) {
  // chech if datagrok connection is available
  const dgConn = await grok.dapi.connections.filter('name="Datagrok"').first();
  if (!dgConn) {
    console.warn('Datagrok connection not found. Datagrok user object handlers not registered');
    return;
  }
  const exp = await DBExplorer.initFromConfigPath(_package);
  if (!exp) {
    grok.shell.error('Failed to load db-explorer config');
    return;
  }
  exp.addEntryPointValueConverter((v) => {
    v ??= '';
    v = v.toString();
    if (v.startsWith('DGUSER-'))
      return v.substring(7);
    return v;
  });
  exp.addCustomRenderer((_, colName, value) => {
    const lc = colName?.toLowerCase() || '';
    return (lc === 'picture') && typeof value === 'string' && !!value && value.startsWith('http');
  }, (value) => imageRenderer(value as string));
  exp.addDefaultHeaderReplacerColumns(['type']); // name is already default

  exp.addCustomRenderer((_, colName, value) => {
    return colName?.toLowerCase() === 'jira_ticket' && typeof value === 'string' && !!value;
  },
  (value) => {
    const jiraUrl = `https://reddata.atlassian.net/browse/${value}`;
    return ui.link(value.toString(), jiraUrl);
  });


  console.log('Datagrok user object handlers registered');
}

export async function newUsersSearch(s: string) {
  const matches = matchAndParseQuery('new users ${}', s);
  if (!matches || !matches[0])
    return null;


  let firstDate: string | null = null;
  let secondDate: string | null = null;

  const curDate = new Date();
  const matchVal = matches[0].toLowerCase();
  if (matchVal === 'today') {
    firstDate = curDate.toISOString().split('T')[0];
    secondDate = curDate.toISOString().split('T')[0];
  } else if (matchVal === 'this month') {
    firstDate = `${curDate.getFullYear()}-${curDate.getMonth() + 1}-01`;
    secondDate = curDate.toISOString().split('T')[0];
  } else if (matchVal === 'this year') {
    firstDate = `${curDate.getFullYear()}-01-01`;
    secondDate = curDate.toISOString().split('T')[0];
  } else if (matchVal === 'last month') {
    const lastMonth = new Date(curDate.getFullYear(), curDate.getMonth() - 1, 3);
    firstDate = `${lastMonth.getFullYear()}-${lastMonth.getMonth() + 1}-01`;
    secondDate = `${curDate.getFullYear()}-${curDate.getMonth() + 1}-01`;
  } else if (matchVal === 'last year') {
    firstDate = `${curDate.getFullYear() - 1}-01-01`;
    secondDate = `${curDate.getFullYear()}-01-01`;
  } else if (matchVal === 'yesterday') {
    const yesterday = new Date(curDate);
    yesterday.setDate(yesterday.getDate() - 1);
    firstDate = yesterday.toISOString().split('T')[0];
    secondDate = yesterday.toISOString().split('T')[0];
  } else {
    const extraMatchMonths = matchAndParseQuery('last ${} months', matchVal);
    const extraMatchYears = matchAndParseQuery('last ${} years', matchVal);
    const extraMatchDays = matchAndParseQuery('last ${} days', matchVal);
    if (extraMatchMonths && extraMatchMonths[0]) {
      const months = parseInt(extraMatchMonths[0]);
      const lastMonth = new Date(curDate.getFullYear(), curDate.getMonth() - months, 3);
      firstDate = `${lastMonth.getFullYear()}-${lastMonth.getMonth() + 1}-01`;
      secondDate = curDate.toISOString().split('T')[0];
    } else if (extraMatchDays && extraMatchDays[0]) {
      const days = parseInt(extraMatchDays[0]);
      const yesterday = new Date(curDate);
      yesterday.setDate(yesterday.getDate() - days);
      firstDate = yesterday.toISOString().split('T')[0];
      secondDate = yesterday.toISOString().split('T')[0];
    } else if (extraMatchYears && extraMatchYears[0]) {
      const years = parseInt(extraMatchYears[0]);
      firstDate = `${curDate.getFullYear() - years}-01-01`;
      secondDate = curDate.toISOString().split('T')[0];
    } else
      return null;
  }

  if (!firstDate || !secondDate)
    return null;

  const connection = await grok.dapi.connections.filter('name="Datagrok"').first();
  if (!connection)
    return null;

  const q = connection.query('newUsers', `
    --name: newUsers
    --output: dataframe result
    select * from users
    where joined between TO_DATE('${firstDate}', 'YYYY-MM-DD') and TO_DATE('${secondDate}', 'YYYY-MM-DD')
    `);

  try {
    const fc = q.prepare();
    return powerSearchQueryTable(fc, {
      postProcess: (tv: DG.TableView) => {
        const col = tv.dataFrame.col('login');
        if (col)
          col.semType = 'DG_USER_LOGIN';
      },
    });
  } catch (e) {
    return null;
  }
}
