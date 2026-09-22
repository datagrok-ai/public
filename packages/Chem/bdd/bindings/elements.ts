/* Chem's own screen parts: the calculators listed on the left of the Chemical Properties dialog
   (a checkbox and a name per registered calculator function). */
import {kind} from '@datagrok-libraries/bdd';

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
