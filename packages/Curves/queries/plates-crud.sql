-- name: getPlates
-- connection: Curves:plates
SELECT * FROM plates.plates;
-- end


-- name: getWellLevelProperties
-- connection: Curves:plates
SELECT DISTINCT p.*
FROM plates.properties p
WHERE p.scope = 'well'
ORDER BY name;
-- end


-- name: getPlateLevelProperties
-- connection: Curves:plates
SELECT DISTINCT p.*
FROM plates.properties p
WHERE p.scope = 'plate'
ORDER BY name;
-- end

-- name: getPropertyNames
-- connection: Curves:plates
SELECT name FROM plates.properties;
-- end


-- name: getPlateTypes
-- connection: Curves:plates
SELECT * FROM plates.plate_types;
-- end


-- name: getPlateTemplates
-- connection: Curves:plates
SELECT * FROM plates.templates;
-- end

-- name: getProperties
-- connection: Curves:plates
-- output: dataframe result
SELECT * FROM plates.properties;
-- end


-- name: getWellRoles
-- connection: Curves:plates
SELECT pav.id, pav.value_string AS name
FROM plates.property_allowed_values pav
JOIN plates.properties p ON pav.property_id = p.id
WHERE p.name = 'Well Role';
-- end


-- name: getWellValuesByBarcode
-- connection: Curves:plates
-- input: string barcode
select v.row, v.col, v.value_num, v.value_string, v.value_bool, v.property_id from plates.plate_well_values v
join plates.plates p on v.plate_id = p.id
where p.barcode = @barcode
order by property_id, row, col
-- end


-- name: getWellValuesById
-- connection: Curves:plates
-- input: int id
select v.row, v.col, v.value_num, v.value_string, v.value_bool, v.property_id from plates.plate_well_values v
join plates.plates p on v.plate_id = p.id
where p.id = @id
order by property_id, row, col
-- end


-- name: getAllowedValues
-- connection: Curves:plates
-- input: string propertyName { choices: getPropertyNames() }
SELECT pav.id, pav.value_string AS name
FROM plates.property_allowed_values pav
JOIN plates.properties p ON pav.property_id = p.id
WHERE p.name = @propertyName;
-- end


-- name: getUniquePlatePropertyValues
-- connection: Curves:plates
SELECT DISTINCT p.name, pd.value_string
FROM plates.plate_details pd
JOIN plates.properties p ON pd.property_id = p.id
WHERE p.type = 'string';
-- end


-- name: getUniqueWellPropertyValues
-- connection: Curves:plates
SELECT DISTINCT p.name, pwv.value_string
FROM plates.plate_well_values pwv
JOIN plates.properties p ON pwv.property_id = p.id
WHERE p.type = 'string';
-- end


-- name: createProperty
-- connection: Curves:plates
-- input: string propertyName
-- input: string valueType
-- input: int templateId { nullable: true }
-- input: string scope
-- input: string choices { nullable: true }
-- input: double min { nullable: true }
-- input: double max { nullable: true }
-- input: int originPlateId { nullable: true }
-- output: int propertyId
INSERT INTO plates.properties(name, type, template_id, scope, choices, min, max, origin_plate_id)
VALUES(@propertyName, @valueType, @templateId, @scope, @choices, @min, @max, @originPlateId)
RETURNING id;
-- end


-- name: createTemplate
-- connection: Curves:plates
-- input: string name
-- input: string description
-- output: int templateId
INSERT INTO plates.templates(name, description)
VALUES(@name, @description)
RETURNING id;
-- end

-- name: getTemplateProperties
-- connection: Curves:plates
-- input: int templateId
-- output: dataframe result
SELECT * FROM plates.properties
WHERE template_id = @templateId;
-- end


-- name: getPlateByBarcode
-- connection: Curves:Plates
-- input: string barcode
-- output: dataframe result
SELECT * FROM plates.plates WHERE barcode = @barcode;
-- end


-- name: createAnalysisRun
-- MODIFIED: Aligns with the new analysis_runs schema.
-- connection: Curves:plates
-- input: int plateId
-- input: string analysisType
-- input: list<string> groups
-- output: int runId
INSERT INTO plates.analysis_runs(plate_id, analysis_type, groups)
VALUES (@plateId, @analysisType, @groups)
RETURNING id;
-- end


-- name: saveAnalysisRunParameter
-- NEW: Saves a single parameter for an analysis run.
-- connection: Curves:plates
-- input: int analysisRunId
-- input: int propertyId
-- input: string valueString {nullable: true}
-- input: double valueNum {nullable: true}
-- input: bool valueBool {nullable: true}
-- input: string valueJsonb {nullable: true}
INSERT INTO plates.analysis_run_parameters(
    analysis_run_id, property_id, value_string, value_num, value_bool, value_jsonb
) VALUES (
    @analysisRunId, @propertyId, @valueString, @valueNum, @valueBool, CAST(@valueJsonb AS jsonb)
);
-- end


-- name: saveAnalysisResult
-- NEW: Saves a single result value for an analysis run.
-- connection: Curves:plates
-- input: int analysisRunId
-- input: list<string> groupCombination
-- input: int propertyId
-- input: string valueString {nullable: true}
-- input: double valueNum {nullable: true}
-- input: bool valueBool {nullable: true}
-- input: string valueJsonb {nullable: true}
INSERT INTO plates.analysis_results(
    analysis_run_id, group_combination, property_id,
    value_string, value_num, value_bool, value_jsonb
) VALUES (
    @analysisRunId, @groupCombination, @propertyId,
    @valueString, @valueNum, @valueBool, CAST(@valueJsonb AS jsonb)
);
-- end
