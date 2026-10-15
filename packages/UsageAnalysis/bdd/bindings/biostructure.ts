/* The Mol* viewport overlay of the Biostructure viewer (BiostructureViewer): its buttons carry no text,
   only a title, and an open panel can show a button whose text is the same title (Settings / Controls
   Info), so the generic button kind, which matches text before titles, can land on the panel's. */
import {kind} from '@datagrok-libraries/bdd';

kind('overlay button', {
  selector: '.msp-viewport-controls-buttons button[title], button.bsv-bs-icon-btn',
  match: ['aria'],
  description: 'a button of the Mol* viewport overlay, by its title ("Toggle Controls Panel" overlay button); selected while Mol* shows it toggled on',
});
