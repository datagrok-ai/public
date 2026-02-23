//@ts-ignore
import initScaffoldNetworkModule from '../scaffold_network.js';

let module: any = null;

async function ensureInit(webRoot: string): Promise<void> {
  if (!module) {
    module = await initScaffoldNetworkModule({
      locateFile: () => `${webRoot}/dist/scaffold_network.wasm`,
    });
  }
}

onmessage = async (event: MessageEvent) => {
  const {smilesJson, ringCutoff, dischargeAndDeradicalize, webRoot}: {
    smilesJson: string,
    ringCutoff: number,
    dischargeAndDeradicalize: boolean,
    webRoot: string,
  } = event.data;

  try {
    await ensureInit(webRoot);
    const result: string = module.generateScaffoldTreeJson(smilesJson, ringCutoff, dischargeAndDeradicalize);
    postMessage({result});
  } catch (e) {
    postMessage({error: e instanceof Error ? e.message : String(e)});
  }
};
