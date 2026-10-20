import * as DG from "datagrok-api/dg";

/** Applies [apply] to the view's CodeMirror once the view has built it. The view that is added
 * becomes current right after, so a view change must not end the wait. */
export function applyCodeMirror(view: DG.View, apply: (codeMirror: any) => void, timeout: number = 10000): void {
    const startTime = new Date().getTime();

    const check = () => {
        const codeMirror = view.root.querySelector(".CodeMirror");
        if (codeMirror)
            apply((codeMirror as any).CodeMirror);
        else if (new Date().getTime() - startTime < timeout)
            window.setTimeout(check, 100);
    };

    check();
}
