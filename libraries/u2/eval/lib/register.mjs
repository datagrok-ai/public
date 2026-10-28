import {register} from 'node:module';
import '../../tests/dom-shim.js';

// the platform's stylesheet link first, so u2 installs its sheet before package CSS, as on a real page
const link = document.createElement('link');
link.setAttribute('rel', 'stylesheet');
link.setAttribute('href', 'datagrok.css');
document.head.append(link);

register('../../tests/dg-stub.mjs', import.meta.url);
register('./hooks.mjs', import.meta.url);
