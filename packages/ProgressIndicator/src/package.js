/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name test
//input: string s
//output: string out
//meta.progress: true
export async function countToTen(s, progress) {
    for(let i = 0; i < 10; i++) {
        if (progress != null)
            progress.update(i * 10);
        await new Promise(resolve => setTimeout(resolve, 1000));
    }
    return `result: ${s}`
}