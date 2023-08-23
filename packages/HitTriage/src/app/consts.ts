export const CampaignIdKey = 'campaignId';
export const HitDesignCampaignIdKey = 'hitDesignCampaignId';
export const CampaignTableName = 'enriched_table.csv';
export const CampaignJsonName = 'campaign.json';
export const HitTriageComputeFunctionTag = 'HitTriageFunction';
export const HitTriageDataSourceTag = 'HitTriageDataSource';
export const HitTriageSubmitTag = 'HitTriageSubmitFunction';
export const HitSelectionColName = 'Selected hits';
export const TileCategoriesColName = 'Stage';
export const HitDesignMolColName = 'Molecule';
export const EmptyStageCellValue = 'StageDefiningRow'; // Kostilj, for a while
export const ViDColName = 'V-iD';
export const ViDColFormat = 'V######';

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
  continueCampaigns: 'Continue Campaigns',
  createNewCampaignHeader: 'Or create new one',
  selectTemplate: 'Template',
} as const;
