CREATE TABLE JSON_DATA (
    data JSON
);

INSERT INTO JSON_DATA(data)
VALUES ('{ "phones":[ {"type": "mobile", "phone": "001001"} , {"type": "fix", "phone": "002002"} ] }');

INSERT INTO JSON_DATA(data) VALUES ('{"bar": "baz", "balance": 7.77, "active":false}');

INSERT INTO JSON_DATA(data) VALUES ('{"reading": 1.230e-5}');