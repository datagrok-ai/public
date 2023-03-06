CREATE TABLE character_types (char_data char(20), varchar_data varchar(max), text_data text, nchar_data nchar(16),
nvarchar_data nvarchar(max), ntext_data ntext);

INSERT INTO character_types (char_data, varchar_data, text_data, nchar_data, nvarchar_data, ntext_data)
VALUES ('Datagrok', 'Hello, World!', 'Lorem ipsum', N'Hello, World!', N'Hello, Datagrok!', 'Hello, Datagrok!');

