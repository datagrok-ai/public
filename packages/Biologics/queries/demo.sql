
--name: assaysByOrganism
--connection: Biologics:biologics
--input: string organism {choices: query("SELECT distinct name FROM biologics.target_organisms")}
--meta.searchPattern: "biologics assays for ${organism}"
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
WHERE org.name = @organism

-- end

--name: adcsLinkedToDrug
--connection: Biologics:biologics
--input: string drugID {choices: query("SELECT distinct identifier FROM biologics.drugs")}
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