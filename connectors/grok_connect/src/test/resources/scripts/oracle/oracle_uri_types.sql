CREATE TABLE uri_types (
    uri URIType
);

INSERT INTO uri_types(uri) VALUES (SYS.URIFACTORY.getURI('/home/oe/doc1.xml'));
INSERT INTO uri_types(uri) VALUES (SYS.URIFACTORY.getURI('/HR/EMPLOYEES/ROW[EMPLOYEE_ID=205]/SALARY'));
INSERT INTO uri_types(uri) VALUES (SYS.URIFACTORY.getURI('https://datagrok.ai'));
