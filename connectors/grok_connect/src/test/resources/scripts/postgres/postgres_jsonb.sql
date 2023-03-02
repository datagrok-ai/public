CREATE TABLE JSONB_DATA (
    data JSONB
);

INSERT INTO JSONB_DATA(data)
VALUES ('{ "phones":[ {"type": "mobile", "phone": "001001"} , {"type": "fix", "phone": "002002"} ] }');

INSERT INTO JSONB_DATA(data) VALUES ('{"bar": "baz", "balance": 7.77, "active":false}');

INSERT INTO JSONB_DATA(data) VALUES ('{"reading": 1.230e-5}');
