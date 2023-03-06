CREATE TABLE null_safety(b BINARY, d DATE, t TIMESTAMP, tm TIME, g GEOGRAPHY, num NUMBER,
                         num10 NUMBER(10,1), fl1 FLOAT, fl2 FLOAT, i NUMBER(9), array1 ARRAY,
                         variant1 VARIANT,
                         object1 OBJECT, first_name VARCHAR(50),  bool BOOLEAN);

INSERT INTO null_safety VALUES (NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL);
