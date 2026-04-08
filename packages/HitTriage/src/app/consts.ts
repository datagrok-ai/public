import type {HitDesignCampaign} from './types';

// Re-export shared compute-function constants from the statistics library
export {funcTypeNames, HTScriptPrefix, HTQueryPrefix, ComputeQueryMolColName} from '@datagrok-libraries/statistics/src/compute-functions/consts';

export const CampaignIdKey = 'campaignId';
export const HitDesignCampaignIdKey = 'campaignId';
export const CampaignTableName = 'enriched_table.csv';
export const CampaignJsonName = 'campaign.json';
export const HitTriageComputeFunctionTag = 'HitTriageFunction';
export const HitDesignerFunctionTag = 'HitDesignerFunction';
export const HitTriageDataSourceTag = 'HitTriageDataSource';
export const HitTriageSubmitTag = 'HitTriageSubmitFunction';
export const HitSelectionColName = 'Selected hits';
export const TileCategoriesColName = 'Stage';
export const HitDesignMolColName = 'Molecule';
export const PeptiHitHelmColName = 'Helm';
export const ViDColName = 'V-iD';
export const ViDSemType = 'HIT_DESIGN_VID';
export const HTcampaignName = 'HTcampaignName';
export const HDcampaignName = 'HDcampaignName';
export const i18n = {
  startNewCampaign: 'New Campaign',
  createNewCampaign: 'New Campaign',
  dataSourceFunction: 'Source',
  createNewTemplate: 'New Template',
  StartCampaign: 'Start',
  createTemplate: 'Create',
  createCampaign: 'Create',
  download: 'Download',
  cancel: 'Cancel',
  continueCampaigns: 'Continue Campaign',
  createNewCampaignHeader: 'New Campaign',
  selectTemplate: 'Template',
  noInformation: 'No Information',
} as const;

export const HDCampaignsGroupingLSKey = 'HDCampaignsGrouping';
export const HDCampaignTableColumnsLSKey = 'HDCampaignTableColumns';
export const HDCampaignsTableSortingLSKey = 'HDCampaignsTableSorting';

export enum CampaignGrouping {
  None = 'None',
  Template = 'Template',
  Status = 'Status',
  Author = 'Author',
  LastModifiedUser = 'Last Modified User',
}

export type CampaignGroupingType = CampaignGrouping | `campaignFields.${string}`;

export const DefaultCampaignTableInfoGetters = {
  'Name': (info: HitDesignCampaign) => info.friendlyName ?? info.name,
  'Code': (info: HitDesignCampaign) => info.name,
  'Created': (info: HitDesignCampaign) => info.createDate,
  'Author': (info: HitDesignCampaign) => info.authorUserFriendlyName ?? '',
  'Last Modified by': (info: HitDesignCampaign) => info.lastModifiedUserName ?? '',
  'Molecules': (info: HitDesignCampaign) => (info.rowCount ?? 0).toString(),
  'Status': (info: HitDesignCampaign) => info.status,
} as const;

export type CampaignTableColumns = keyof typeof DefaultCampaignTableInfoGetters | `campaignFields.${string}`;

export const HTFunctionOrderingLSKey = 'HTFunctionOrderingLS';
