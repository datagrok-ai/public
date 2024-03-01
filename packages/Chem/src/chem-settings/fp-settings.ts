import { Fingerprint } from '../utils/chem-common';

export let _fpParamsMap: Map<Fingerprint, FpParams> = new Map<Fingerprint, FpParams>();

export type FpParam = {
    uiName: string,
    value: number,
    tooltip: string,
}

export type FpParams = MorganFpParams | AtomPairFpParams | RDKitFpParams | TopologicalTorsionFpParams;

export class FpParamsBase {
    constructor() {};

    get params(): { [key: string]: number } {
        const obj: { [key: string]: number } = {};
        Object.keys(this).forEach(key => obj[key] = (this as any)[key].value);
        return obj;
    }

    get paramsAsString(): string {
        return JSON.stringify(this.params);
    }

    applyParams(params?: { [key: string]: number }) {
        if (params) {
            Object.keys(params).forEach((key) => {
                if (this.hasOwnProperty(key))
                    (this as any)[key].value = params[key];
            });
        }
    };
}

export class MorganFpParams extends FpParamsBase {
    radius: FpParam = { uiName: 'Radius', value: 2, tooltip: '' };
    numBits: FpParam = { uiName: 'Bits', value: 2048, tooltip: '' };
    constructor(params: { [key: string]: number }) {
        super();
        this.applyParams(params);
    };
}

export class AtomPairFpParams extends FpParamsBase {
    minLength: FpParam = { uiName: 'Min length', value: 1, tooltip: '' };
    maxLength: FpParam = { uiName: 'Max length', value: 30, tooltip: '' };
    numBits: FpParam = { uiName: 'Bits', value: 2048, tooltip: '' };
    constructor(params: { [key: string]: number }) {
        super();
        this.applyParams(params);
    };
}

export class RDKitFpParams extends FpParamsBase {
    maxPath: FpParam = { uiName: 'Max path', value: 7, tooltip: '' };
    numBits: FpParam = { uiName: 'Bits', value: 2048, tooltip: '' };
    constructor(params: { [key: string]: number }) {
        super();
        this.applyParams(params);
    };
}

export class TopologicalTorsionFpParams extends FpParamsBase {
    numBits: FpParam = { uiName: 'Bits', value: 2048, tooltip: '' };
    constructor(params: { [key: string]: number }) {
        super();
        this.applyParams(params);
    };
}

export const fpBuilderMap: { [key: string]: (params: { [key: string]: number }) => FpParams } = {
    [Fingerprint.Morgan]: (params) => new MorganFpParams(params),
    [Fingerprint.AtomPair]: (params) => new AtomPairFpParams(params),
    [Fingerprint.RDKit]: (params) => new RDKitFpParams(params),
    [Fingerprint.TopologicalTorsion]: (params) => new TopologicalTorsionFpParams(params),
    //no params can be passed for MACCS fp according to rdkit minilib api
}