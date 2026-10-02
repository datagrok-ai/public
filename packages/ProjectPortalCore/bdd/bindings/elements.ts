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

app.element('open data button', {selector: u2('open-data-button')});
app.element('open assay results button', {selector: u2('open-dataset-assays')});

app.element('needs attention', {selector: u2('attention')});
app.element('needs attention heading', {selector: u2('attention-heading')});
app.element('hERG not measured count', {selector: u2('attention-rule-herg-unmeasured')});
app.element('lead without solubility count', {selector: u2('attention-rule-lead-no-solubility')});
app.element('EXC-1006 attention item', {selector: u2('attention-item-exc-1006')});
app.element('EXC-1006 attention link', {selector: u2('attention-open-exc-1006')});

app.element('sync history button', {selector: u2('sync-history-button')});
app.element('sync history dialog', {selector: u2('sync-history-dialog')});
app.element('sync schedule', {selector: u2('sync-schedule')});
app.element('sync owners', {selector: u2('sync-owners')});
app.element('latest sync run', {selector: u2('sync-run-0')});
app.element('sync auto-advance', {selector: u2('sync-auto-advance')});
app.element('preview advance by data button', {selector: u2('auto-advance-preview-button')});
app.element('advance by data preview', {selector: u2('auto-advance-preview')});

app.element('export button', {selector: u2('export-button')});
app.element('export dialog', {selector: u2('export-dialog')});
app.element('exported portfolio file', {selector: u2('export-json')});

app.element('missing assays heading', {selector: u2('missing-assays-heading')});
app.element('missing solubility', {selector: u2('missing-assay-adme-panel-solubility-um')});
app.element('suggested models heading', {selector: u2('suggested-models-heading')});
app.element('suggested solubility estimate', {selector: u2('suggested-model-projectportalcore-solubilityestimate')});
app.element('solubility estimate card link', {selector: u2('recommended-card-projectportalcore-solubilityestimate')});

app.element('follow button', {selector: u2('follow-button')});
app.element('what\'s new', {selector: u2('whats-new')});
app.element('KIN-A follow link', {selector: u2('follow-series-kin-a')});

app.element('mode input', {selector: u2('mode-input')});
app.element('persona input', {selector: u2('persona-input')});
app.element('quick actions', {selector: u2('quick-actions')});
app.element('quick actions heading', {selector: u2('quick-actions-heading')});
app.element('more button', {selector: u2('more-button')});
app.element('saved runs heading', {selector: u2('saved-runs-heading')});

app.element('program roles', {selector: u2('program-roles')});
app.element('role principal input', {selector: u2('role-principal-input')});
app.element('program role input', {selector: u2('role-input')});
app.element('add role button', {selector: u2('add-role-button')});
app.element('role result', {selector: u2('role-result')});
app.element('remove admin role link', {selector: u2('remove-role-admin')});

app.element('evidence project card', {selector: u2('program-bdd-evidence')});
app.element('decision evidence', {selector: u2('decision-evidence')});
app.element('evidence readiness checkbox', {selector: u2('evidence-readiness')});
app.element('evidence link input', {selector: u2('evidence-link-input')});
app.element('attach evidence button', {selector: u2('add-evidence-link-button')});
app.element('rationale input', {selector: u2('rationale-input')});
app.element('latest decision evidence', {selector: u2('decision-evidence-0')});
app.element('evidence snapshot checkbox', {selector: u2('evidence-snapshot')});
app.element('program snapshots', {selector: u2('program-snapshots')});
app.element('latest decision snapshot link', {selector: '[data-u2-name^="evidence-0-snapshot-"]'});
app.element('snapshot dialog', {selector: u2('snapshot-dialog')});
app.element('snapshot stage', {selector: u2('snapshot-stage')});
app.element('snapshot criteria', {selector: u2('snapshot-criteria')});
app.element('snapshot table', {selector: u2('snapshot-table')});

app.element('scratch docs project card', {selector: u2('program-bdd-docs')});
app.element('documents heading', {selector: u2('documents-heading')});
app.element('process description document card', {selector: u2('document-card-process-description')});
app.element('link process description document button', {selector: u2('link-document-process-description')});
app.element('link document dialog', {selector: u2('link-document-dialog')});
app.element('document file input', {selector: u2('document-file-input')});
app.element('document title input', {selector: u2('document-title-input')});
app.element('document version input', {selector: u2('document-version-input')});
app.element('link document button', {selector: u2('link-document-button')});
app.element('approve process description document link', {selector: u2('approve-document-process-description-0')});

app.element('provenance project card', {selector: u2('program-bdd-prov')});
app.element('program title', {selector: u2('program-title')});
app.element('latest decision provenance link', {selector: u2('provenance-decision-0')});
app.element('provenance dialog', {selector: u2('provenance-dialog')});
app.element('provenance start', {selector: u2('provenance-start')});
app.element('provenance rests on', {selector: u2('provenance-up')});
app.element('provenance used by', {selector: u2('provenance-down')});
app.element('first used-by record', {selector: u2('provenance-down-0-open')});

app.element('feed project card', {selector: u2('program-bdd-feed')});
app.element('activity heading', {selector: u2('feed-heading')});
app.element('activity entries', {selector: u2('feed-entries')});
app.element('activity kind input', {selector: u2('feed-kind-input')});
app.element('my feed button', {selector: u2('my-feed-button')});
app.element('my feed dialog', {selector: u2('my-feed-dialog')});
app.element('my feed entries', {selector: u2('my-feed-entries')});
app.element('first my feed program link', {selector: u2('my-feed-program-0')});

app.element('score solubility estimate button', {selector: u2('score-model-projectportalcore-solubilityestimate')});
app.element('score result', {selector: u2('score-result')});
app.element('predictions heading', {selector: u2('predictions-heading')});
app.element('solubility predictions fit', {selector: u2('prediction-fit-projectportalcore-solubilityestimate-solubility-um')});
app.element('missing solubility prediction', {selector: u2('missing-assay-adme-panel-solubility-um-predicted')});
app.element('solubility predicted vs measured button', {selector: u2('predicted-vs-measured-projectportalcore-solubilityestimate-solubility-um')});

app.element('approval project card', {selector: u2('program-bdd-approval')});
app.element('binding criterion met button', {selector: u2('mark-met-binding')});
app.element('approval setting', {selector: u2('approval-setting')});
app.element('approval result', {selector: u2('approval-result')});
app.element('pending approvals', {selector: u2('pending-approvals')});
app.element('no pending approvals', {selector: u2('no-pending-approvals')});
app.element('closed approvals', {selector: u2('closed-approvals')});
app.element('first approval approve button', {selector: u2('approval-0-approve')});
app.element('first approval withdraw button', {selector: u2('approval-0-withdraw')});
app.element('waiting for me button', {selector: u2('waiting-button')});
app.element('waiting dialog', {selector: u2('waiting-dialog')});
app.element('waiting list', {selector: u2('waiting-list')});

app.element('workspace button', {selector: u2('workspace-button')});
app.element('workspace dialog', {selector: u2('workspace-dialog')});
app.element('orphans', {selector: u2('orphans')});
app.element('park preview button', {selector: u2('park-preview-button')});
app.element('park result', {selector: u2('park-result')});
app.element('park scratch folder link', {selector: u2('park-orphan-bddparkscratch')});
app.element('parked folders', {selector: u2('parked-folders')});
app.element('no parked folders', {selector: u2('no-parked')});
app.element('unpark scratch folder button', {selector: u2('unpark-bddparkscratch')});

app.element('scratch evidence project card', {selector: u2('program-bdd-et')});
app.element('generate process description evidence button', {selector: u2('generate-evidence-process-description')});
app.element('scratch stale evidence project card', {selector: u2('program-bdd-et2')});
app.element('regenerate process description evidence link', {selector: u2('regenerate-document-process-description-0')});
