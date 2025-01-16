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
export const ViDColFormat = 'V######';
export const HTcampaignName = 'HTcampaignName';
export const HDcampaignName = 'HDcampaignName';
export const HTScriptPrefix = 'HTScript';
export const HTQueryPrefix = 'HTQuery';
export const ComputeQueryMolColName = 'molecules';
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

export const funcTypeNames = {
  script: 'script',
  function: 'function-package',
  query: 'data-query',
} as const;

export const HDCampaignsGroupingLSKey = 'HDCampaignsGrouping';

export enum CampaignGroupingType {
  None = 'None',
  Template = 'Template',
  Status = 'Status',
  Author = 'Author',
  LastModifiedUser = 'Last Modified User',
}
