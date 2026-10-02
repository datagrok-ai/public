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
app.element('work units', {selector: u2('work-units')});
app.element('KIN-A series', {selector: u2('unit-series-kin-a')});
app.element('TOX-201 study', {selector: u2('unit-study-tox-201')});
app.element('compound roles', {selector: u2('compound-roles')});
app.element('EXC-3001 compound', {selector: u2('compound-exc-3001')});

app.element('table button', {selector: u2('table-button')});
app.element('all portfolios button', {selector: u2('all-portfolios-button')});

app.element('scratch runs project card', {selector: u2('program-bdd-runs')});
app.element('saved runs', {selector: u2('saved-runs')});
app.element('example model link', {selector: u2('model-link-projectportalcore-fedbatchtiter')});
app.element('model card', {selector: u2('model-card')});
app.element('model card question', {selector: u2('model-card-question')});
app.element('model card inputs', {selector: u2('model-card-inputs')});

app.element('scratch save project card', {selector: u2('program-bdd-save')});
app.element('save run link', {selector: u2('save-run-link')});
app.element('save run dialog', {selector: u2('save-run-dialog')});
app.element('run portfolio input', {selector: u2('portfolio-input')});
app.element('run program input', {selector: u2('program-input')});
app.element('run unit operation input', {selector: u2('unit-operation-input')});
app.element('run type input', {selector: u2('run-type-input')});
app.element('run purpose input', {selector: u2('purpose-input')});
app.element('save run button', {selector: u2('save-run-button')});
app.element('save a run button', {selector: u2('save-a-run-button')});

app.element('scratch pick project card', {selector: u2('program-bdd-pick')});
app.element('process description criterion', {selector: u2('criterion-process-description')});
app.element('process description run', {selector: u2('selected-run-process-description')});
app.element('process description run history', {selector: u2('run-history-process-description')});
app.element('choose process description run button', {selector: u2('choose-run-process-description')});
app.element('choose run dialog', {selector: u2('choose-run-dialog')});
app.element('run filter', {selector: u2('run-filter')});
app.element('feed 2 candidate', {selector: u2('candidate-feed-2')});
app.element('feed 4 candidate', {selector: u2('candidate-feed-4')});
app.element('run comparison', {selector: u2('run-comparison')});
app.element('selected run input', {selector: u2('selected-run-input')});
app.element('selection rationale input', {selector: u2('selection-rationale-input')});
app.element('select run button', {selector: u2('select-run-button')});
