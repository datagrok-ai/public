/* The names the portal's features use. Everything the app draws carries a data-u2-name derived
   from a code (src/views/portal-view.ts: `slug`), so each name below is one attribute selector. */
import {context} from '@datagrok-libraries/bdd';

const u2 = (name: string) => `[data-u2-name="${name}"]`;

export const app = context('Project Portal app', {selector: u2('Project Portal'), aliases: ['portal']});

app.element('portfolio list', {selector: u2('portfolio-list')});
app.element('small molecule discovery portfolio', {selector: u2('portfolio-ex-smd'), aliases: ['EX-SMD portfolio']});
app.element('biologics discovery portfolio', {selector: u2('portfolio-ex-biod'), aliases: ['EX-BIOD portfolio']});
app.element('biologics development portfolio', {selector: u2('portfolio-ex-bioc'), aliases: ['EX-BIOC portfolio']});

app.element('stage board', {selector: u2('stage-board'), aliases: ['board']});
app.element('board header', {selector: u2('board-header')});
app.element('unstaged column', {selector: u2('stage-unstaged')});
app.element('gate D4 column', {selector: u2('stage-d4')});
app.element('target validation column', {selector: u2('stage-tv'), aliases: ['TV column']});

app.element('EXS-KIN1 card', {selector: u2('program-exs-kin1')});
app.element('EXB-101 card', {selector: u2('program-exb-101')});
app.element('EXB-201 card', {selector: u2('program-exb-201')});
app.element('scratch project card', {selector: u2('program-bdd-scratch')});
app.element('KINASE-X card', {selector: u2('program-kinase-x')});

app.element('program detail', {selector: u2('program-detail')});
app.element('stage track', {selector: u2('stage-track')});
app.element('gate criteria', {selector: u2('gate-criteria')});
app.element('gate reasons', {selector: u2('gate-reasons')});
app.element('transition history', {selector: u2('transition-history')});
app.element('advance button', {selector: u2('advance-button')});
app.element('open workspace button', {selector: u2('open-workspace-button')});
app.element('assay criterion met button', {selector: u2('mark-met-assay-ready')});
app.element('gate heading', {selector: u2('gate-heading')});
app.element('readiness score', {selector: u2('readiness-score')});
app.element('H2L step', {selector: u2('step-h2l')});
app.element('LO step', {selector: u2('step-lo')});
app.element('potency criterion', {selector: u2('criterion-potency')});
app.element('ADME panel criterion detail', {selector: u2('criterion-detail-adme-panel')});
