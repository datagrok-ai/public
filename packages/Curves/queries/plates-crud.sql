-- name: getPlates
-- connection: Admin:Plates
SELECT * FROM plates.plates;
-- end


-- name: getWellLevelProperties
-- connection: Admin:Plates
SELECT DISTINCT p.id, p.name, p.value_type, p.unit
FROM plates.properties p
JOIN plates.plate_well_values pd ON p.id = pd.property_id
ORDER BY p.name;
-- end


-- name: getPlateLevelProperties
-- connection: Admin:Plates
SELECT DISTINCT p.id, p.name, p.value_type, p.unit
FROM plates.properties p
JOIN plates.plate_details pd ON p.id = pd.property_id
ORDER BY p.name;
-- end


-- name: getPropertyNames
-- connection: Admin:Plates
SELECT name FROM plates.properties;
-- end


-- name: getPlateTypes
-- connection: Admin:Plates
SELECT * FROM plates.plate_types;
-- end


-- name: getPlateTemplates
-- connection: Admin:Plates
SELECT * FROM plates.templates;
-- end


-- name: getWellRoles
-- connection: Admin:Plates
SELECT pav.id, pav.value_string AS name
FROM plates.property_allowed_values pav
JOIN plates.properties p ON pav.property_id = p.id
WHERE p.name = 'Well Role';
-- end


-- name: getWellValuesByBarcode
-- connection: Admin:Plates
-- input: string barcode
select v.row, v.col, v.value_num, v.value_string, v.value_bool, v.property_id from plates.plate_well_values v
join plates.plates p on v.plate_id = p.id
where p.barcode = @barcode
order by property_id, row, col
-- end


-- name: getAllowedValues
-- connection: Admin:Plates
-- input: string propertyName { choices: getPropertyNames() }
SELECT pav.id, pav.value_string AS name
FROM plates.property_allowed_values pav
JOIN plates.properties p ON pav.property_id = p.id
WHERE p.name = @propertyName;
-- end


-- name: getUniquePlatePropertyValues
-- connection: Admin:Plates
-- input: int propertyId propertyName { choices: getPropertyNames() }
SELECT DISTINCT p.name, pd.value_string
FROM plates.plate_details pd
JOIN plates.properties p ON pd.property_id = p.id
WHERE p.value_type = 'string';
-- end


-- name: getUniqueWellPropertyValues
-- connection: Admin:Plates
-- input: int propertyId propertyName { choices: getPropertyNames() }
SELECT DISTINCT p.name, pwv.value_string
FROM plates.plate_well_values pwv
JOIN plates.properties p ON pwv.property_id = p.id
WHERE p.value_type = 'string';
-- end


-- name: createProperty
-- connection: Admin:Plates
-- input: string propertyName
-- input: string valueType
-- output: int propertyId
INSERT INTO plates.properties(name, value_type)
VALUES(@propertyName, @valueType)
RETURNING id;
-- end


-- name: createTemplate
-- connection: Admin:Plates
-- input: string name
-- input: string description
-- output: int templateId
INSERT INTO plates.templates(name, description)
VALUES(@name, @description)
RETURNING id;
-- end
