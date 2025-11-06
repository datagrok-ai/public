--name: getMonomerByGrokID
--connection: MonomerDB:mdb1
--input: string grokID
--meta.searchPattern: "monomer with ID ${grokID}"
SELECT 
    m.id,
    m.symbol,
    m.name,
    m.monomer_type AS "monomerType",
    m.polymer_type AS "polymerType",
    m.smiles,
    m.molblock AS "molfile",
    m.capped_smiles AS "cappedSmiles",
    m.natural_analog AS "naturalAnalog",
    m.created_by AS "author",
    m.created_at AS "createDate",
    -- R-groups array
    (
        SELECT COALESCE(json_agg(
            json_build_object(
                'label', rg.label,
                'capGroupSMILES', rg.cap_group_smiles,
                'capGroupName', rg.cap_group_name,
                'alternateId', rg.label || '-' || COALESCE(rg.cap_group_name, ''),
                'description', rg.description
            )
            ORDER BY rg.label
        )::text, '[]')
        FROM (
            SELECT r1.id, r1.cap_group_smiles, r1.cap_group_name, r1.description, 'R1' as label 
            FROM mdb1.r_groups r1 WHERE r1.id = m.r1_id
            UNION ALL
            SELECT r2.id, r2.cap_group_smiles, r2.cap_group_name, r2.description, 'R2' as label 
            FROM mdb1.r_groups r2 WHERE r2.id = m.r2_id
            UNION ALL
            SELECT r3.id, r3.cap_group_smiles, r3.cap_group_name, r3.description, 'R3' as label 
            FROM mdb1.r_groups r3 WHERE r3.id = m.r3_id
            UNION ALL
            SELECT r4.id, r4.cap_group_smiles, r4.cap_group_name, r4.description, 'R4' as label 
            FROM mdb1.r_groups r4 WHERE r4.id = m.r4_id
            UNION ALL
            SELECT r5.id, r5.cap_group_smiles, r5.cap_group_name, r5.description, 'R5' as label 
            FROM mdb1.r_groups r5 WHERE r5.id = m.r5_id
            UNION ALL
            SELECT r6.id, r6.cap_group_smiles, r6.cap_group_name, r6.description, 'R6' as label 
            FROM mdb1.r_groups r6 WHERE r6.id = m.r6_id
            UNION ALL
            SELECT r7.id, r7.cap_group_smiles, r7.cap_group_name, r7.description, 'R7' as label 
            FROM mdb1.r_groups r7 WHERE r7.id = m.r7_id
            UNION ALL
            SELECT r8.id, r8.cap_group_smiles, r8.cap_group_name, r8.description, 'R8' as label 
            FROM mdb1.r_groups r8 WHERE r8.id = m.r8_id
        ) rg
    ) AS "rgroups",
    -- Meta object containing all properties
    (
        SELECT COALESCE(json_build_object(
            'properties', (
                SELECT COALESCE(json_object_agg(
                    mp.name,
                    CASE 
                        WHEN mp.property_type = 'number' THEN to_jsonb(mp.value_number)
                        WHEN mp.property_type = 'boolean' THEN to_jsonb(mp.value_boolean)
                        ELSE to_jsonb(mp.value_string)
                    END
                ), '{}'::json)
                FROM mdb1.monomer_properties mp
                WHERE mp.monomer_id = m.id
            )
        )::text, '{}')
    ) AS "meta"
FROM mdb1.monomers m
WHERE m.grok_identifier = @grokID;
--end

--name: upsertMonomer
--connection: MonomerDB:mdb1
--input: string symbol
--input: string name
--input: string monomerType
--input: string smiles
--input: string molblock
--input: string libraryName
--input: string username
--input: string polymerType {nullable: true}
--input: string cappedSmiles {nullable: true}
--input: string naturalAnalog {nullable: true}
--input: int r1Id {nullable: true}
--input: int r2Id {nullable: true}
--input: int r3Id {nullable: true}
--input: int r4Id {nullable: true}
--input: int r5Id {nullable: true}
--input: int r6Id {nullable: true}
--input: int r7Id {nullable: true}
--input: int r8Id {nullable: true}
--output: int id
INSERT INTO mdb1.monomers (
    symbol,
    name,
    monomer_type,
    polymer_type,
    smiles,
    molblock,
    capped_smiles,
    natural_analog,
    library_id,
    r1_id,
    r2_id,
    r3_id,
    r4_id,
    r5_id,
    r6_id,
    r7_id,
    r8_id,
    created_by,
    updated_by
)
VALUES (
    @symbol,
    @name,
    @monomerType,
    @polymerType,
    @smiles,
    @molblock,
    @cappedSmiles,
    @naturalAnalog,
    (SELECT id FROM mdb1.monomer_libraries WHERE name = @libraryName),
    @r1Id,
    @r2Id,
    @r3Id,
    @r4Id,
    @r5Id,
    @r6Id,
    @r7Id,
    @r8Id,
    @username,
    @username
)
ON CONFLICT (symbol, library_id) 
DO UPDATE SET
    name = EXCLUDED.name,
    monomer_type = EXCLUDED.monomer_type,
    polymer_type = EXCLUDED.polymer_type,
    smiles = EXCLUDED.smiles,
    molblock = EXCLUDED.molblock,
    capped_smiles = EXCLUDED.capped_smiles,
    natural_analog = EXCLUDED.natural_analog,
    r1_id = EXCLUDED.r1_id,
    r2_id = EXCLUDED.r2_id,
    r3_id = EXCLUDED.r3_id,
    r4_id = EXCLUDED.r4_id,
    r5_id = EXCLUDED.r5_id,
    r6_id = EXCLUDED.r6_id,
    r7_id = EXCLUDED.r7_id,
    r8_id = EXCLUDED.r8_id,
    updated_by = @username,
    updated_at = NOW()
RETURNING id;
--end

--name: createLibrary
--connection: MonomerDB:mdb1
--input: string name
--input: string username
--input: string friendlyName {nullable: true}
--input: string description {nullable: true}
--input: string libraryType {nullable: true}
--input: string source {nullable: true}
--input: string version {nullable: true}
--output: int id
INSERT INTO mdb1.monomer_libraries (
    name,
    friendly_name,
    description,
    library_type,
    source,
    version,
    created_by,
    updated_by
)
VALUES (
    @name,
    @friendlyName,
    @description,
    @libraryType,
    @source,
    @version,
    @username,
    @username
)
ON CONFLICT (name)
DO UPDATE SET
    friendly_name = EXCLUDED.friendly_name,
    description = EXCLUDED.description,
    library_type = EXCLUDED.library_type,
    source = EXCLUDED.source,
    version = EXCLUDED.version,
    updated_by = @username,
    updated_at = NOW()
RETURNING id;
--end

--name: upsertMonomerProperty
--connection: MonomerDB:mdb1
--input: int monomerId
--input: string name
--input: string username
--input: string propertyType
--input: string valueString {nullable: true}
--input: double valueNumber {nullable: true}
--input: bool valueBoolean {nullable: true}
--output: int id
INSERT INTO mdb1.monomer_properties (
    monomer_id,
    name,
    property_type,
    value_string,
    value_number,
    value_boolean,
    created_by,
    updated_by
)
VALUES (
    @monomerId,
    @name,
    @propertyType,
    @valueString,
    @valueNumber,
    @valueBoolean,
    @username,
    @username
)
ON CONFLICT (monomer_id, name)
DO UPDATE SET
    property_type = EXCLUDED.property_type,
    value_string = EXCLUDED.value_string,
    value_number = EXCLUDED.value_number,
    value_boolean = EXCLUDED.value_boolean,
    updated_by = @username,
    updated_at = NOW()
RETURNING id;
--end

--name: getMonomersByLibrary
--connection: MonomerDB:mdb1
--input: string libraryName
--meta.searchPattern: "get all monomers from library ${libraryName}"
SELECT 
    m.id,
    m.symbol,
    m.name,
    m.monomer_type AS "monomerType",
    m.polymer_type AS "polymerType",
    m.smiles,
    m.molblock AS "molfile",
    m.capped_smiles AS "cappedSmiles",
    m.natural_analog AS "naturalAnalog",
    m.created_by AS "author",
    m.created_at AS "createDate",
    -- R-groups array
    (
        SELECT COALESCE(json_agg(
            json_build_object(
                'label', rg.label,
                'capGroupSMILES', rg.cap_group_smiles || '[*:' || SUBSTRING(rg.label FROM 2) || ']',
                'capGroupName', rg.cap_group_name,
                'alternateId', rg.label || '-' || COALESCE(rg.cap_group_name, ''),
                'description', rg.description
            )
            ORDER BY rg.label
        )::text, '[]')
        FROM (
            SELECT r1.id, r1.cap_group_smiles, r1.cap_group_name, r1.description, 'R1' as label 
            FROM mdb1.r_groups r1 WHERE r1.id = m.r1_id
            UNION ALL
            SELECT r2.id, r2.cap_group_smiles, r2.cap_group_name, r2.description, 'R2' as label 
            FROM mdb1.r_groups r2 WHERE r2.id = m.r2_id
            UNION ALL
            SELECT r3.id, r3.cap_group_smiles, r3.cap_group_name, r3.description, 'R3' as label 
            FROM mdb1.r_groups r3 WHERE r3.id = m.r3_id
            UNION ALL
            SELECT r4.id, r4.cap_group_smiles, r4.cap_group_name, r4.description, 'R4' as label 
            FROM mdb1.r_groups r4 WHERE r4.id = m.r4_id
            UNION ALL
            SELECT r5.id, r5.cap_group_smiles, r5.cap_group_name, r5.description, 'R5' as label 
            FROM mdb1.r_groups r5 WHERE r5.id = m.r5_id
            UNION ALL
            SELECT r6.id, r6.cap_group_smiles, r6.cap_group_name, r6.description, 'R6' as label 
            FROM mdb1.r_groups r6 WHERE r6.id = m.r6_id
            UNION ALL
            SELECT r7.id, r7.cap_group_smiles, r7.cap_group_name, r7.description, 'R7' as label 
            FROM mdb1.r_groups r7 WHERE r7.id = m.r7_id
            UNION ALL
            SELECT r8.id, r8.cap_group_smiles, r8.cap_group_name, r8.description, 'R8' as label 
            FROM mdb1.r_groups r8 WHERE r8.id = m.r8_id
        ) rg
    ) AS "rgroups",
    -- Meta object containing all properties
    (
        SELECT COALESCE(json_build_object(
            'properties', (
                SELECT COALESCE(json_object_agg(
                    mp.name,
                    CASE 
                        WHEN mp.property_type = 'number' THEN to_jsonb(mp.value_number)
                        WHEN mp.property_type = 'boolean' THEN to_jsonb(mp.value_boolean)
                        ELSE to_jsonb(mp.value_string)
                    END
                ), '{}'::json)
                FROM mdb1.monomer_properties mp
                WHERE mp.monomer_id = m.id
            )
        )::text, '{}')
    ) AS "meta"
FROM mdb1.monomers m
JOIN mdb1.monomer_libraries ml ON m.library_id = ml.id
WHERE ml.name = @libraryName
ORDER BY m.symbol;
--end

--name: listLibraries
--connection: MonomerDB:mdb1
SELECT name FROM mdb1.monomer_libraries
--end

