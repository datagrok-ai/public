import * as grokNamespace from './js-api/grok';
import * as uiNamespace from './js-api/ui';
import * as DGNamespace from './js-api/dg';

declare global {
    const grok: typeof grokNamespace;
    const ui: typeof uiNamespace;
    const DG: typeof DGNamespace;
}