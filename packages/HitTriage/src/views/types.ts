import * as DG from 'datagrok-api/dg';

export type IDescriptorTree = {
    [key: string]: {
      descriptors: Array<Descriptor>,
    } & Descriptor;
}

export type Descriptor = {
    name: string,
    description: string,
};

export type IFunctionInput = {
    name: string,
    displayName: string,
    description?: string | null,
    type: DG.TYPE,
    defaultValue?: any,
    choices?: string[] | null,
    value?: any
}

export type IFunctionProps = {
    name: string,
    package: string,
    displayName: string,
    description?: string | null,
    inputs: IFunctionInput[],
    calculate: boolean,
    onlyCheckboxes: boolean,
}

export type IComputeDialogResult = {
    descriptors: string[],
    externals: {
        [_: string]: {[key: string]: any}
    }
}
