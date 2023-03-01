CREATE TABLE xml_type (xml_data xml);
INSERT INTO xml_type(xml_data) VALUES (XML('<foo>Hello World!</foo>'));
INSERT INTO xml_type(xml_data) VALUES (XML('<?xml version="1.0"?><book><title>Manual</title><chapter>...</chapter></book>'));
