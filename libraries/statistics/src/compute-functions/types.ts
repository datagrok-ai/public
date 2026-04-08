import * as DG from 'datagrok-api/dg';

export type IFunctionArgs<T = any> = {
    [key: string]: T,
}

export type TemplateFunction = {
    package: string,
    name: string,
    args: IFunctionArgs,
}

export type TemplateScript = {
    name: string,
    args: IFunctionArgs,
    id: string,
}

export type ComputeQuery = TemplateScript & {
    inputName: string,
}

export type TemplateCompute = {
    descriptors: {
        enabled: boolean,
        args: string[],
    }
    functions: TemplateFunction[],
    scripts?: TemplateScript[],
    queries?: TemplateScript[],
};

export type ComputeFunctions = {
    functions: DG.Func[],
    scripts: DG.Script[],
    queries: DG.DataQuery[],
};

export type IComputeDialogResult = {
    descriptors: string[],
    externals: {
        [_: string]: IFunctionArgs
    },
    scripts?: {
        [_: string]: IFunctionArgs
    }
    queries?: {
        [_: string]: IFunctionArgs
    }
}

export type IDescriptorTree = {
    [key: string]: {
      descriptors: Array<Descriptor>,
    } & Descriptor;
}

export type Descriptor = {
    name: string,
    description: string,
};

export type IChemFunctionsDialogResult = {
    okProxy: () => void,
    root: HTMLElement,
};
