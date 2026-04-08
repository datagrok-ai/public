import * as DG from 'datagrok-api/dg';

// Re-export shared compute-function types from the statistics library
export type {TemplateFunction, TemplateScript, ComputeQuery, TemplateCompute,
  ComputeFunctions, IComputeDialogResult, IFunctionArgs,
  IDescriptorTree, Descriptor, IChemFunctionsDialogResult} from '@datagrok-libraries/statistics/src/compute-functions/types';
import type {TemplateFunction, TemplateScript, TemplateCompute} from '@datagrok-libraries/statistics/src/compute-functions/types';

// Backward-compatible aliases
export type HitTriageTemplateFunction = TemplateFunction;
export type HitTriageTemplateScript = TemplateScript;
export type HitTriageComputeQuery = TemplateScript & { inputName: string };
export type HitTriageTemplateCompute = TemplateCompute;

export type AppName = keyof CampaignsType;

export type CampaignsType = {
    'Hit Triage': HitTriageCampaign,
    'Hit Design': HitDesignCampaign,
    'PeptiHit': HitDesignCampaign,
}

export type HitTriageTemplate = {
    name: string,
    key: string,
    compute: TemplateCompute,
    submit?: HitTriageTemplateSubmit,
    campaignFields: HitTriageCampaignField[],
    dataSourceType: IngestType,
    queryFunctionName?: string,
    isDataSourceQuery?: boolean,
    layoutViewState?: string,
    localLayoutPath?: string,
}

export const CampaignFieldTypes = {
  String: DG.TYPE.STRING,
  Number: DG.TYPE.FLOAT,
  Boolean: DG.TYPE.BOOL,
  Date: DG.TYPE.DATE_TIME,
  Molecule: DG.SEMTYPE.MOLECULE,
} as const;

export type HitTriageCampaignFieldType = keyof typeof CampaignFieldTypes;

export type HitTriageCampaignField =
    {name: string, type: HitTriageCampaignFieldType, semtype?: string, required?: boolean};

export type IngestType = 'File' | 'Query';

export type HitTriageTemplateIngest = {
    type: IngestType,
    query: string,
    molColName: string,
};

export type HitTriageTemplateSubmit = {
    fName: string,
    package?: string | undefined
};

export type HitTriageCampaignStatus = string;

export type HitTriageCampaign = {
    name: string,
    friendlyName?: string,
    templateName: string,
    template?: HitDesignTemplate,
    savePath?: string,
    status: HitTriageCampaignStatus,
    createDate: string,
    campaignFields: {[key: string]: any},
    filters: {[key: string]: any}[],
    ingest: HitTriageTemplateIngest,
    columnSemTypes?: {[key: string]: string},
    rowCount?: number,
    filteredRowCount?: number,
    layout?: string,
    columnTypes?: {[key: string]: string},
    version?: number,
    permissions?: TriagePermissions,
    authorUserId?: string,
    authorUserFriendlyName?: string,
    lastModifiedUserName?: string,
};

export type TriagePermissions = {
    // will store the ids of UserGroups that have access to the campaign
    edit: string[],
    view: string[],
};

export type INewTemplateResult<T> = {
    template: Promise<T>,
    root: HTMLElement,
    cancelPromise: Promise<void>,
}

// ##################### HIT DESIGN TYPES #####################

export type HitDesignTemplate = Omit<HitTriageTemplate, 'dataSourceType' | 'queryFunctionName'> &
    {stages: string[]};

export type HitDesignCampaign = Omit<HitTriageCampaign, 'filters' | 'ingest'> & {tilesViewerFormSketch?: string};

// todo: probably add some more stuff
export type PeptiHitTemplate = HitDesignTemplate & {toAtomiLevelProps?: {[key: string]: any}}
