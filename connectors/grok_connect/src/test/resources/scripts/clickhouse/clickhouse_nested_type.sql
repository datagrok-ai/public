CREATE TABLE NESTED_TYPE
(
    EventDate Date,
    UserID UInt64,
    Attrs Nested(
        Key String,
        Value String)
);

INSERT INTO NESTED_TYPE VALUES ('2016-01-01', 123, ['price', 'color'], ['high', 'red']);
