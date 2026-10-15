/* What FileEditors' PDF preview draws: one canvas per page, up to ten, filled as the pages render. */
import {element} from '@datagrok-libraries/bdd';

element('pdf preview page', {selector: '.fe-pdf-preview canvas',
  description: 'a page of a PDF file\'s preview (FileEditors), a canvas the page is drawn on'});
