import * as DG from "datagrok-api/dg";
import {Subscription} from "rxjs";
import * as grok from "datagrok-api/grok";
import {first} from "rxjs/operators";

export function applyCodeMirror(view: DG.View, apply: (codeMirror: any) => void, timeout: number = 1500): void {
    const startTime = new Date().getTime();
    let timeoutId: number | undefined;

    const sub: Subscription = grok.events.onViewChanged.pipe(first()).subscribe((_) => {
        if (timeoutId)
            clearTimeout(timeoutId);
    });

    const check = () => {
        if (new Date().getTime() - startTime >= timeout) {
            sub.unsubscribe();
            return;
        }
        const codeMirror = view.root.querySelector(".CodeMirror");
        if (codeMirror) {
            // @ts-ignore
            apply(codeMirror.CodeMirror);
            sub.unsubscribe();
            return;
        }

        timeoutId = window.setTimeout(check, 100);
    };

    check();
}