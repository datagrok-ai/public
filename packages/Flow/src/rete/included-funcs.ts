/** The Flow toolbox is **allowlist-based**: a function appears as a pipeline
 *  node only when its `Func.nqName` (namespace-qualified, e.g. `Chem:foo`) is
 *  in this set, OR it declares **`meta.includeInFlow: true`** on itself, OR it
 *  is one of the always-included kinds — saved flows, **queries**
 *  (`DG.DataQuery`, except from dev/test packages), SetVar / GetVar — see
 *  `shouldIncludeFunc` in `node-factory.ts`. `meta.includeInFlow: false` still
 *  opts out even an allowlisted function.
 *
 *  This list was produced empirically on 2026-07-20 by dumping the catalog the
 *  previous filter pipeline (structural rules + role/tag rules + a curated
 *  denylist) produced on a live server — i.e. it IS that catalog, frozen as an
 *  explicit include list (minus queries, which are included by kind). Edit
 *  freely: add an nqName to surface a function, remove one to hide it. For a
 *  function you own, prefer declaring `meta.includeInFlow: true` in its own
 *  metadata instead of editing this. */
export const INCLUDED_FUNC_NQNAMES: ReadonlySet<string> = new Set<string>([
  // core
  'core:AddNewColumn',
  'core:AddNewColumnList',
  // 'core:AddStickyPropertiesColumns',
  // 'core:Aggregate',
  'core:AppendTables',
  // 'core:ApplyModel',
  // 'core:BinByDateTime',
  // 'core:BinBySpecificLimits',
  'core:ChangeColumnsType',
  'core:ChangeType',
  'core:CloneColumn',
  'core:CloneTable',
  // 'core:Cluster',
  // 'core:CompareTables',
  // 'core:Date',
  // 'core:DateAdd',
  // 'core:DateDiff',
  // 'core:DateDifference',
  // 'core:DateDifferenceWithColumn',
  // 'core:DateFromUnixTimestamp',
  // 'core:DateNow',
  'core:DateParse',
  // 'core:DateTime',
  'core:DayOfMonth',
  'core:DayOfWeek',
  'core:DayOfYear',
  // 'core:DeleteColumns',
  // 'core:DeleteRows',
  'core:ExtractColumns',
  'core:ExtractNumbers',
  // 'core:ExtractRows',
  'core:ExtractValue',
  // 'core:FilterRows',
  // 'core:FormatColumns',
  // 'core:GetColumn',
  'core:GetTable',
  // 'core:Hour',
  'core:InsertRows',
  'core:InvertSelection',
  'core:JoinTables',
  // 'core:KeepColumns',
  // 'core:KeepRows',
  // 'core:KSTest',
  'core:LinkTables',
  // 'core:Millisecond',
  // 'core:Minute',
  // 'core:MissingValuesImputation',
  'core:Month',
  'core:Normalize',
  'core:OpenFile',
  'core:OpenTable',
  'core:ParseQnum',
  // 'core:PCA',
  'core:Percentile',
  'core:Qnum',
  'core:Quarter',
  // 'core:RandomData',
  'core:RenameColumn',
  'core:ReplaceValues',
  // 'core:RoundNumbers',
  // 'core:SaveAsColumn',
  // 'core:Second',
  'core:SelectAll',
  // 'core:SelectCols',
  'core:SelectionToColumn',
  'core:SelectNone',
  // 'core:SelectRows',
  'core:SplitByRegExp',
  'core:SplitColumn',
  // 'core:Subset',
  // 'core:TagColumns',
  'core:Time',
  'core:TimeParse',
  'core:Today',
  // 'core:TransformString',
  'core:TTest',
  'core:UnixTimestamp',
  'core:Unpivot',
  'core:Weeknum',
  'core:Year',
  // Admetica
  // 'Admetica:run_admetica',
  // Source-only as of 2026-07-30: the deployed Admetica predates it, so the
  // entry is inert until the package is republished.
  'Admetica:getAdmeProperties',
  // Bio
  'Bio:getRegionTopMenu',
  'Bio:moleculesToHelmTopMenu',
  'Bio:pepseaMsa',
  'Bio:sequenceIdentityScoring',
  'Bio:sequenceSimilarityScoring',
  'Bio:splitToMonomersTopMenu',
  'Bio:toAtomicLevel',
  // Bionemo
  // 'Bionemo:DiffdockPython',
  // BiostructureViewer
  'BiostructureViewer:fetchSequencesFromPdb',
  'BiostructureViewer:importPdb',
  'BiostructureViewer:importPdbqt',
  'BiostructureViewer:importWithNgl',
  'BiostructureViewer:importXYZ',
  'BiostructureViewer:ProteinLigandInteractionDiagram',
  // Chem
  'Chem:addChemPropertiesColumns',
  'Chem:addChemRisksColumns',
  // Adds the InChI / InChI Key column to the table; the `getInchis` /
  // `getInchiKeys` twins return a detached column instead, which is less useful
  // on a canvas (nothing downstream can bind it back to its rows).
  'Chem:addInchisKeysTopMenu',
  'Chem:addInchisTopMenu',
  'Chem:applyReaction',
  'Chem:bitbirchClusteringTopMenu',
  'Chem:ButinaMoleculesClustering',
  'Chem:CalculateLogD',
  'Chem:CalculateLogP',
  'Chem:CalculateLogS',
  'Chem:CalculatePI',
  'Chem:CalculatePKa',
  // Parameterized metric/fingerprint, and it takes the table, so the query
  // molecule gets a sketcher — unlike `findSimilar`, which hard-codes Morgan.
  'Chem:callChemSimilaritySearch',
  // 'Chem:ChemicalSpaceUsingTSNE',
  // 'Chem:ChemicalSpaceUsingUMAP',
  'Chem:chemSpaceColumns',
  'Chem:clusterMCSTopMenu',
  'Chem:convertNotation',
  'Chem:Curate',
  'Chem:deprotect',
  'Chem:descriptorsDocker',
  'Chem:FilterByCatalogs',
  // 'Chem:findSimilar',  → superseded by Chem:callChemSimilaritySearch
  'Chem:GenerateConformers',
  'Chem:getDiversities',
  // 'Chem:getInchiKeys',  → superseded by Chem:addInchisKeysTopMenu
  // 'Chem:getInchis',     → superseded by Chem:addInchisTopMenu
  // 'Chem:getMorganFingerprints',  bit-string column nothing downstream reads
  'Chem:getSimilarities',
  // The four file importers take file *content*, which nothing on a canvas can
  // produce until `core:ReadFileBytes` exists (Phase 6).
  // 'Chem:importMol',
  // 'Chem:importMol2',
  // 'Chem:importSdf',
  // 'Chem:importSmi',
  'Chem:mpoScoreByProfile',
  'Chem:MurckoScaffolds',
  'Chem:Mutate',
  'Chem:namesToSmiles',
  'Chem:pharmacophoreFeaturesTopMenu',
  'Chem:recalculateCoords',
  'Chem:removeWaterAndSaltsTopMenu',
  'Chem:runElementalAnalysis',
  'Chem:searchSubstructure',
  'Chem:similarityMatrixTopMenu',
  'Chem:structuralAlertsTopMenu',
  // Both need a synthon library uploaded first, and offer no readiness signal —
  // the node shows an empty combo with no explanation.
  // 'Chem:SynthonSearch',
  // 'Chem:synthonSearchFunc',
  'Chem:toSdf',
  'Chem:TwoComponentReaction',
  // 'Chem:USRCAT',
  // Chembl
  'Chembl:getChemblCompounds',
  'Chembl:getChemblCompoundsByOrganism',
  // ChemblApi
  // 'ChemblApi:chemblSearchWidget',
  'ChemblApi:getById',
  // Chemspace
  'Chemspace:getChemspaceIds',
  'Chemspace:getChemspacePrices',
  // Curves
  'Curves:CalculateMSR',
  'Curves:pzfxFileHandler',
  // DiffStudio
  'DiffStudio:ballFlight',
  'DiffStudio:ivpModel_Acid_Production',
  'DiffStudio:ivpModel_Bioreactor',
  'DiffStudio:ivpModel_PK_PD',
  'DiffStudio:ivpModel_Pollution',
  'DiffStudio:solve',
  'DiffStudio:solveEquations',
  'DiffStudio:solveODE',
  // Eda
  // 'Eda:applyLinearKernelSVM',
  // 'Eda:applyLinearRegression',
  // 'Eda:applyPLSRegression',
  // 'Eda:applyPolynomialKernelSVM',
  // 'Eda:applyRBFkernelSVM',
  // 'Eda:applySigmoidKernelSVM',
  // 'Eda:applySoftmax',
  // 'Eda:applyXGBooster',
  'Eda:dbScan',
  'Eda:generatePmpoDataset',
  'Eda:kNNImputationForTable',
  'Eda:MCLClustering',
  'Eda:PCA',
  'Eda:PLS',
  // FileEditors
  // 'FileEditors:viewDocx',
  // 'FileEditors:viewPdf',
  // Proteomics
  'Proteomics:DeqmsDE',
  'Proteomics:LimmaDE',
  'Proteomics:VsnNormalize',
  // SequenceTranslator
  'SequenceTranslator:combineSenseAntisenseToOligoNucleotide',
  'SequenceTranslator:convertHelmToOligoNucleotide',
  'SequenceTranslator:enumerateSingleHelmSequence',
  'SequenceTranslator:enumerateSingleHelmSequenceWithNaturalAAs',
  // Sequenceutils
  'Sequenceutils:helmMsa',
  // Surechembl
  'Surechembl:sureChemblSimilaritySearch',
  'Surechembl:sureChemblSubstructureSearch',]);
