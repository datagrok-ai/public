export interface ScaffoldNetworkModule {
  generateScaffoldTreeJson(smilesArrayJson: string, ringCutoff: number, dischargeAndDeradicalize: boolean): string;
}

export default function initScaffoldNetworkModule(
  overrides?: {locateFile?: (path: string) => string}
): Promise<ScaffoldNetworkModule>;
