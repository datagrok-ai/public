--name: ADCs Linked to Drug
--friendlyName: ADCs Linked to Drug
--connection: Biologics:biologics
--input: string drugID {choices: query("SELECT distinct identifier FROM biologics.drugs"); semType: DG_BIOLOGICS_DRUG_ID} [Drug identifier]
--description: "Find ADCs in biologics database linked to a specified drug identifier."
--meta.searchPattern: "adcs linked to ${drugID}"
--meta.role: panel
--output: dataframe result
SELECT adc.identifier AS adc_identifier, adc.name AS adc_name, adc.glyph AS glyph,
    seq.identifier as sequence_identifier, seq.name AS antibody_name,
    seq.heavy_chain AS antibody_heavy_chain, seq.light_chain AS antibody_light_chain,
    l.identifier as linker_identifier, l.linker_type, l.linker_molecule_smiles, l.linker_sequence,
    d.identifier as drug_identifier, d.name AS drug_name, d.smiles AS drug_smiles
FROM biologics.adc adc
JOIN biologics.sequences seq ON adc.antibody_id = seq.id
JOIN biologics.drugs d ON adc.drug_id = d.id
JOIN biologics.linkers l ON adc.linker_id = l.id
WHERE d.identifier = @drugID
AND (seq.heavy_chain IS NOT NULL OR seq.light_chain IS NOT NULL);

-- end