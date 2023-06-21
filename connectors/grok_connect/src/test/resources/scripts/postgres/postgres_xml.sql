CREATE TABLE XML_DATA (
    data XML
);

INSERT INTO XML_DATA (data) VALUES (XML('<foo>Hello World!</foo>'));
INSERT INTO XML_DATA (data) VALUES (XML('abc<foo>bar</foo><bar>foo</bar>'));
INSERT INTO XML_DATA (data) VALUES (XML('<?xml version="1.0"?><book><title>Manual</title><chapter>...</chapter></book>'));
