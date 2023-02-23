CREATE TABLE XML_DATA (
    data XMLType
);

INSERT INTO XML_DATA (data) VALUES ('<foo>Hello World!</foo>');
INSERT INTO XML_DATA (data) VALUES ('<book><title>Manual</title><chapter>...</chapter></book>');
