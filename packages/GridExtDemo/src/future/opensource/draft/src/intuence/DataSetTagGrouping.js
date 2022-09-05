export class DataSetTagGrouping
{
 constructor(strName)
 {
     this.m_strName = strName;
     DataSetTagGrouping.MAP_ENTRIES.set(strName, this);
 }
}

DataSetTagGrouping.MAP_ENTRIES = new Map();
DataSetTagGrouping.get = function(strName)
{
    this.m_strName = strName;
}


DataSetTagGrouping.KEY = "DATASET_GROUPING";
DataSetTagGrouping.COMPOUND_SELECTORS = new DataSetTagGrouping("COMPOUND_SELECTORS");
DataSetTagGrouping.ANALYSIS_ASSAY = new DataSetTagGrouping("ANALYSIS_ASSAY");
DataSetTagGrouping.CHEMISTRY_INFO = new DataSetTagGrouping("CHEMISTRY_INFO");
DataSetTagGrouping.CHEM_ENRICH = new DataSetTagGrouping("CHEM_ENRICH");
DataSetTagGrouping.CHEM_ADMET = new DataSetTagGrouping("ADMET");