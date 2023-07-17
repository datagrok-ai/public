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
    args: IFunctionArgs,
}

export type IComputeDialogResult = {
    descriptors: string[],
    externals: {
        [_: string]: IFunctionArgs
    }
}

export type ITemplate = {
    name: string,
    key: string,
    compute: ITemplateCompute,
    submit?: ITemplateSubmit,
    campaignFields: ICampaignField[],
    dataSourceType: IngestType,
}

export const CampaignFieldTypes = {
  String: DG.TYPE.STRING,
  Number: DG.TYPE.FLOAT,
  Boolean: DG.TYPE.BOOL,
  Date: DG.TYPE.DATE_TIME,
} as const;

export type ICampaignFieldType = keyof typeof CampaignFieldTypes;

export type ICampaignField = {name: string, type: ICampaignFieldType, required?: boolean};

export type IngestType = 'File' | 'Query';

export type ITemplateIngest = {
    type: IngestType,
    query: string,
    molColName: string,
};

export type ITemplateCompute = {
    descriptors: {
        enabled: boolean,
        args: string[],
    }
    functions: ITemplateFunction[],
};

export type ITemplateSubmit = {
    fName: string,
    package: string
};

export type ICampaignStatus = 'In Progress' | 'Submitted';

export type ICampaign = {
    name: string,
    templateName: string,
    status: ICampaignStatus,
    createDate: string,
    campaignFields: {[key: string]: any},
    filters: {[key: string]: any}[],
    ingest: ITemplateIngest,
};

export type IChemFunctionsDialogResult = {
    okProxy: () => void,
    root: HTMLElement,
};
