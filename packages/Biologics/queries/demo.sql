
--name: assaysByOrganism
--connection: Biologics:biologics
--input: string organism {choices: query("SELECT distinct name FROM biologics.target_organisms")}
--description: "Find biologics assays for a specified organism."
--meta.searchPattern: "biologics assays for ${organism}"
SELECT org.name AS organism, org.identifier as organism_identifier, at.name AS assay_type,
    ar.result_value, ar.units, seq.sequence AS antibody_sequence, seq.identifier as sequence_identifier,
    seq.name AS antibody_name, adc.name as adc_name, adc.identifier as adc_identifier, adc.glyph AS glyph,
    d.name AS drug_name, d.identifier as drug_identifier, d.smiles AS drug_smiles,
    l.identifier as linker_identifier, l.linker_type, l.linker_molecule_smiles, l.linker_sequence,
    org.id as orgID, at.id as assay_type_id, seq.id as seqId, ar.id as assay_res_id, adc.id as adc_id, d.id as drug_id
FROM biologics.assay_results ar
LEFT JOIN biologics.assay_types at ON ar.assay_id = at.id
LEFT JOIN biologics.target_organisms org ON ar.target_organism_id = org.id
LEFT JOIN biologics.adc adc ON ar.adc_id = adc.id
LEFT JOIN biologics.sequences seq ON adc.antibody_id = seq.id
LEFT JOIN biologics.drugs d ON adc.drug_id = d.id
LEFT JOIN biologics.linkers l ON adc.linker_id = l.id
WHERE org.name ilike @organism

-- end

--name: ADCsWithCapsazeActivityHigherThan
--connection: Biologics:biologics
--description: "Find biologics ADCs with caspase activity higher than a specified value."
--input: double minActivity
--meta.searchPattern: "biologics ADCs with caspase activity higher than ${minActivity}"
SELECT org.name AS organism, org.identifier as organism_identifier, at.name AS assay_type,
    ar.result_value, ar.units, seq.sequence AS antibody_sequence, seq.identifier as sequence_identifier,
    seq.name AS antibody_name, adc.name as adc_name, adc.identifier as adc_identifier, adc.glyph AS glyph,
    d.name AS drug_name, d.identifier as drug_identifier, d.smiles AS drug_smiles,
    l.identifier as linker_identifier, l.linker_type, l.linker_molecule_smiles, l.linker_sequence, l.linker_type
FROM biologics.assay_results ar
LEFT JOIN biologics.assay_types at ON ar.assay_id = at.id
LEFT JOIN biologics.target_organisms org ON ar.target_organism_id = org.id
LEFT JOIN biologics.adc adc ON ar.adc_id = adc.id
LEFT JOIN biologics.sequences seq ON adc.antibody_id = seq.id
LEFT JOIN biologics.drugs d ON adc.drug_id = d.id
LEFT JOIN biologics.linkers l ON adc.linker_id = l.id
WHERE at.name ILIKE 'caspase activity' AND ar.result_value > @minActivity;

-- end

--name: ADCsWithIC50HLThan
--connection: Biologics:biologics
--description: "Find biologics ADCs with IC50 higher or lower than a specified value. Use 'higher' or 'lower' for valueTarget to indicate the comparison direction."
--input: double value
--input: string valueTarget {choices: ['higher', 'lower']}
--meta.searchPattern: "biologics ADCs with IC50 ${valueTarget} than ${value}"
SELECT org.name AS organism, org.identifier as organism_identifier, at.name AS assay_type,
    ar.result_value, ar.units, seq.sequence AS antibody_sequence, seq.identifier as sequence_identifier,
    seq.name AS antibody_name, adc.name as adc_name, adc.identifier as adc_identifier, adc.glyph AS glyph,
    d.name AS drug_name, d.identifier as drug_identifier, d.smiles AS drug_smiles,
    l.identifier as linker_identifier, l.linker_type, l.linker_molecule_smiles, l.linker_sequence, l.linker_type
FROM biologics.assay_results ar
LEFT JOIN biologics.assay_types at ON ar.assay_id = at.id
LEFT JOIN biologics.target_organisms org ON ar.target_organism_id = org.id
LEFT JOIN biologics.adc adc ON ar.adc_id = adc.id
LEFT JOIN biologics.sequences seq ON adc.antibody_id = seq.id
LEFT JOIN biologics.drugs d ON adc.drug_id = d.id
LEFT JOIN biologics.linkers l ON adc.linker_id = l.id
WHERE at.name ILIKE 'IC50' AND (
    ((@valueTarget ILIKE 'lower' or @valueTarget ILIKE 'less') AND ar.result_value < @value) OR
    ((@valueTarget ILIKE 'higher' or @valueTarget ILIKE 'more') AND ar.result_value > @value)
);

-- end


--name: adcsLinkedToDrug
--connection: Biologics:biologics
--input: string drugID {choices: query("SELECT distinct identifier FROM biologics.drugs")}
--description: "Find ADCs in biologics database linked to a specified drug identifier."
--meta.searchPattern: "adcs linked to ${drugID}"
SELECT adc.identifier AS adc_identifier, adc.name AS adc_name, adc.glyph AS glyph,
    seq.identifier as sequence_identifier, seq.name AS antibody_name, seq.sequence AS antibody_sequence,
    l.identifier as linker_identifier, l.linker_type, l.linker_molecule_smiles, l.linker_sequence,
    d.identifier as drug_identifier, d.name AS drug_name, d.smiles AS drug_smiles
FROM biologics.adc adc
JOIN biologics.sequences seq ON adc.antibody_id = seq.id
JOIN biologics.drugs d ON adc.drug_id = d.id
JOIN biologics.linkers l ON adc.linker_id = l.id
WHERE d.identifier = @drugID;

-- end

--name: getBiologicsPeptideHelmByIdentifier
--connection: Biologics:biologics
--description: "Retrieve the HELM sequence of a biologics peptide given its identifier that should follow the pattern GROKPEP-######."
--input: string peptideIdentifier
--output: string result {semType: Macromolecule; units: helm}
SELECT helm
FROM biologics.peptides
WHERE identifier = @peptideIdentifier
limit 1;
-- end
