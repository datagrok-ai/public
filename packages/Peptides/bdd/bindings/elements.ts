import {dataset, element} from '@datagrok-libraries/bdd';

dataset('peptides', {path: 'System:DemoFiles/bio/peptides.csv',
  description: '647 aligned separator peptides, ID and IC50; 17 positions including NH2 and COOH'});

element('Peptides landing view', {selector: '.grok-view-container > .grok-view'});
