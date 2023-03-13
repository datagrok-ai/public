CREATE TABLE JSON_TYPE(json_type json);

INSERT INTO JSON_TYPE(json_type) VALUES ('{"key1": "value1", "key2": "value2"}');
INSERT INTO JSON_TYPE(json_type) VALUES ('{ "phones":[ {"type": "mobile", "phone": "001001"} , {"type": "fix", "phone": "002002"} ] }');
INSERT INTO JSON_TYPE(json_type) VALUES ('{"reading": 1.230e-5}');
