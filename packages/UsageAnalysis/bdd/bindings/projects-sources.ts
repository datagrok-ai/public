/* What only the Projects features of projects built from a database, a query or a script need: the
   Save project dialog's creation-script text and table rows, and queries on the platform's own
   database. */
import type {Page} from '@playwright/test';
import {element, Given, When} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The text a table's Creation script button of the Save project dialog unfolds: a span next to the
   button, hidden until the button is clicked. */
element('creation script text', {selector: '.ui-creation-script > span'});

/* The table rows of the Save project dialog's entity list, without the rows of its views and
   connections (the "project table" kind takes every row). */
element('saved table row', {selector: '.d4-dialog [name^="project-table-"]'});

const listOf = (text: string): string[] => text.split(',').map((s) => s.trim()).filter(Boolean);

/* A query saved on the platform's own database connection (System:Datagrok, which every user may
   read and query; the entity type names are the platform's, the same on every server), for a
   feature whose claim is what a project does with it, not the query editor. Its --input lines become
   its parameters. The feature removes it through "no query named … is on the server", which runs
   before this and again at feature end. */
export const datagrokQuery = Given('a query {string} on the Datagrok connection is:', async (page: Page, name: string, sql: string) => {
  await page.evaluate(async ([n, s]) => {
    const connection = await grok.dapi.connections.filter('name = "Datagrok"').first();
    if (connection?.nqName !== 'System:Datagrok')
      throw new Error('no System:Datagrok connection on the server');
    await grok.dapi.queries.save(connection.query(n, s));
  }, [name, sql]);
}, {tier: 'api', description: 'saved through the JS API on System:Datagrok, its --input lines becoming its parameters'});

export const changeQuerySql = When('the query {string} on the server is changed to:', async (page: Page, name: string, sql: string) => {
  const inputs = await page.evaluate(async ([n, s]) => {
    const q = await grok.dapi.queries.filter(`name = "${n}"`).first();
    if (!q)
      throw new Error(`no query "${n}" on the server`);
    // the parameters are parsed from the SQL when a query is made, not when its text is set: the
    // edited query is made anew on the connection and saved over the old one
    const connection = await grok.dapi.connections.filter('name = "Datagrok"').first();
    const edited = connection.query(n, s);
    edited.id = q.id;
    await grok.dapi.queries.save(edited);
    return (await grok.dapi.queries.find(q.id)).inputs.map((p: any) => p.name).join(', ');
  }, [name, sql]);
  expect(inputs, `the parameters of the query "${name}" after the edit`).toBe(listOf(
    sql.split('\n').filter((l) => l.startsWith('--input:')).map((l) => l.split(/\s+/)[2]).join(', ')).join(', '));
}, {tier: 'api', description: 'the query\'s SQL (and so its parameters) rewritten on the server, as its editor would save it; the parameters are read back'});
