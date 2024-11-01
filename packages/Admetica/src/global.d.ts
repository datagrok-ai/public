import * as grokNamespace from 'datagrok-api/grok';
import * as uiNamespace from 'datagrok-api/ui';
import * as DGNamespace from 'datagrok-api/dg';
import * as rxjsNamespace from 'rxjs';
import $Namespace from 'cash-dom';

declare global {
    const grok: typeof grokNamespace;
    const ui: typeof uiNamespace;
    const DG: typeof DGNamespace;
    const rjxs: typeof rxjsNamespace;
    const $: typeof $Namespace;
}
