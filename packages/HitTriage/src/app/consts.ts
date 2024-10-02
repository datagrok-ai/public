export const CampaignIdKey = 'campaignId';
export const HitDesignCampaignIdKey = 'campaignId';
export const CampaignTableName = 'enriched_table.csv';
export const CampaignJsonName = 'campaign.json';
export const HitTriageComputeFunctionTag = 'HitTriageFunction';
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
  startNewCampaign: 'New campaign',
  createNewCampaign: 'New campaign',
  dataSourceFunction: 'Source',
  createNewTemplate: 'New template',
  StartCampaign: 'Start',
  createTemplate: 'Create',
  createCampaign: 'Create',
  download: 'Download',
  cancel: 'Cancel',
  continueCampaigns: 'Continue a campaign',
  createNewCampaignHeader: 'New campaign',
  selectTemplate: 'Template',
} as const;

export const funcTypeNames = {
  script: 'script',
  function: 'function-package',
  query: 'data-query',
} as const;
