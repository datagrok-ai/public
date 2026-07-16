import * as DG from 'datagrok-api/dg';

export type IFunctionArgs<T = any> = {
    [key: string]: T,
}

export type TemplateFunction = {
    package: string,
    name: string,
    args: IFunctionArgs,
    // When true, the function is recomputed over the whole molecule column every time a
    // campaign that uses this compute config is opened (e.g. database lookups that resolve later).
    rerunOnOpen?: boolean,
}

export type TemplateScript = {
    name: string,
    args: IFunctionArgs,
    id: string,
    // See TemplateFunction.rerunOnOpen.
    rerunOnOpen?: boolean,
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
    // Keyed by the same function key used in externals/scripts/queries; only enabled functions
    // that are flagged to re-run on campaign open are present. Absent when the caller did not
    // enable the re-run-on-open option in the dialog.
    rerunOnOpen?: {
        [_: string]: boolean
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
