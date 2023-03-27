CREATE OR REPLACE TABLE test_number(num NUMBER,
                                    num10 NUMBER(10,1), fl1 FLOAT, fl2 FLOAT, i NUMBER(9));

INSERT INTO test_number (num, num10, fl1, fl2, i) VALUES (99999999999999999999999999999999999999, 255.5, 'NaN', 0.2, 100);
INSERT INTO test_number (num, num10, fl1, fl2, i) VALUES (-99999999999999999999999999999999999999, 15.5, 'inf', +1.34, 51200);
INSERT INTO test_number (num, num10, fl1, fl2, i) VALUES (9, 255.2, '-inf', 317, 7780);
INSERT INTO test_number (num, num10, fl1, fl2, i) VALUES (10000, 44444.4, 1.234E+2, 124.412, 12);
