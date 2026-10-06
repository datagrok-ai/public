/* Chem's own screen parts: the calculators listed on the left of the Chemical Properties dialog
   (a checkbox and a name per registered calculator function), the MPO Profiles app — the
   editable title and description of a profile tab, its property rows and the rows of the list —
   and the 3D view of a molecule. */
import {element, kind} from '@datagrok-libraries/bdd';

kind('calculator', {
  selector: '.biochem-calc-nav-item',
  match: ['text'],
  editorSelector: 'input[type="checkbox"]',
  description: 'a calculator row of the Chemical Properties dialog, by its name ("Chemical Properties (OCL)", "logP")',
});

kind('reaction', {
  aliases: ['reaction card'],
  selector: '.d4-dialog-contents .d4-flex-col:has(> canvas + label)',
  match: ['label'],
  labelSelector: ':scope > label:first-of-type',
  description: 'a reaction card of the Run Reaction and Two-Component Reaction dialogs, by its name ("Amide Coupling")',
});

element('MPO profile title', {selector: '.chem-profile-header',
  description: 'the editable name at the top of a profile tab; "Untitled Profile" until one is typed'});
element('MPO profile description', {selector: '.chem-profile-description',
  description: 'the editable description under the profile title'});

kind('MPO property', {
  selector: '.statistics-mpo-row',
  match: ['text'],
  editorSelector: '.statistics-mpo-property-cell input',
  description: 'a property row of the profile editor; its editor is the property-name field ("first MPO property")',
});

kind('MPO profile', {
  selector: '.chem-mpo-profiles-table tr:has(.chem-mpo-actions-button)',
  match: ['label'],
  labelSelector: ':scope > td:nth-child(2)',
  parts: {actions: '.chem-mpo-actions-button'},
  description: 'a row of the Manage Profiles list, by the profile name; "actions of X MPO profile" is its ⋮ button',
});

element('R-Groups settings icon', {selector: '.chem-rgroup-settings-icon',
  description: 'the gear of the R-Groups Analysis dialog: shows Matching strategy and Only match at R groups, which the dialog remembers per account'});

/* The SAR Matrix Summary tab's own parts. A Summary row is a plain div carrying its heading in its
   own child, and the segments and tabs are plain divs rather than tab controls, so the platform's
   `tab` kind does not see them. */
kind('summary row', {selector: '.chem-sar-sum-row', labelSelector: '.chem-sar-card-name',
  match: ['label', 'text'], description: 'a clickable finding on the SAR Matrix Summary tab'});

kind('summary segment', {selector: '.chem-sar-sum-seg', match: ['text'],
  description: 'one of the Summary tab\'s segments (Overview, Effects, Worth making, Method)'});

kind('tier chip', {selector: '.chem-sar-sum-seg-tiers .chem-sar-sum-role',
  labelSelector: '.chem-sar-sum-tier-name', match: ['label'],
  description: 'a fold tier (All, L1, L2, ...) the Summary tab reads its rankings at'});

kind('effects tab', {selector: '.chem-sar-sum-sub', labelSelector: '.chem-sar-sum-sub-name',
  match: ['label'],
  description: 'one component column\'s tab on the Summary tab\'s Effects segment'});
