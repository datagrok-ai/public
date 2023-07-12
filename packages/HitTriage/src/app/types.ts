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

export type IFunctionArgs = {
    [key: string]: any,
}

export type ITemplateFunction = {
    package: string,
    name: string,
    // args: IFunctionArgs,
}

export type IComputeDialogResult = {
    descriptors: string[],
    externals: {
        [_: string]: IFunctionArgs
    }
}

export type ITemplate = {
    name: string,
    compute: ITemplateCompute,
    submit?: ITemplateSubmit,
}

export type IngestType = 'File' | 'Function';

export type ITemplateIngest = {
    type: IngestType,
    query: string,
    molColName: string,
};

export type ITemplateCompute = {
    descriptors: {
        enabled: boolean,
    }
    functions: ITemplateFunction[],
};

export type ITemplateSubmit = {
    fName: string,
    package: string
};

export type ICampaign = {
    name: string,
    templateName: string,
    filters: {[key: string]: any}[],
    ingest: ITemplateIngest,
};
