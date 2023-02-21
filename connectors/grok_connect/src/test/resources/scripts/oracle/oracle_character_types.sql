CREATE TABLE character_type (
                   ch CHAR(10),
                   varch VARCHAR2(10),
                   nch NCHAR(10),
                   nvarch NVARCHAR2(50)
);

INSERT INTO character_type(ch, varch, nch, nvarch) VALUES ('Hello', 'World', 'Datagrok', 'Groking');