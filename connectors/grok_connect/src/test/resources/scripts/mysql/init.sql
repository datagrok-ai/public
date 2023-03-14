CREATE TABLE BINARY_TYPES(binary_type BINARY(5), varbinary_type VARBINARY(10));

INSERT INTO BINARY_TYPES(binary_type, varbinary_type) VALUES (BINARY('Hello'), BINARY('datagrok'));

CREATE TABLE BIT_TYPE(bit_type1 bit(64), bit_type2 bit(1));

INSERT INTO BIT_TYPE(bit_type1, bit_type2) VALUES (b'10000000', b'1');

INSERT INTO BIT_TYPE(bit_type1, bit_type2) VALUES (b'1', b'0');

create table mock_data (
                           id bigint,
                           first_name VARCHAR(50),
                           last_name VARCHAR(50),
                           email VARCHAR(50),
                           gender VARCHAR(50),
                           ip_address VARCHAR(50),
                           bool BOOLEAN,
                           country VARCHAR(50),
                           date DATE,
                           some_number decimal(5, 2)
);

insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (1, 'Burk', 'Kemery', 'bkemery0@businesswire.com', 'Male', '249.64.22.121/32', true, 'China', '2017-09-20', 510.32);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (2, 'Nicholle', 'Karoly', 'nkaroly1@alexa.com', 'Female', '255.233.247.118/32', false, 'Poland', '2014-02-27', 864.09);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (3, 'Orlando', 'Westgate', 'owestgate2@dedecms.com', 'Polygender', '75.0.252.254/32', false, 'Netherlands', '2020-09-03', 822.7);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (4, 'Gothart', 'Cokayne', 'gcokayne3@plala.or.jp', 'Male', '196.83.12.163/32', true, 'Philippines', '2001-01-31', 251.05);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (5, 'Mitchell', 'Haglington', 'mhaglington4@indiegogo.com', 'Male', '209.93.181.190/32', true, 'Poland', '2020-10-09', 15.22);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (6, 'Jeromy', 'Twinn', 'jtwinn5@globo.com', 'Male', '25.13.2.132/32', true, 'Serbia', '2014-10-04', 378.4);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (7, 'Joela', 'Cornau', 'jcornau6@imgur.com', 'Female', '195.47.88.236/32', false, 'Indonesia', '2020-03-19', 349.11);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (8, 'Darren', 'Juares', 'djuares7@hexun.com', 'Male', '94.170.16.96/32', false, 'China', '2011-04-09', 631.89);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (9, 'Marlie', 'Mayze', 'mmayze8@google.com.au', 'Female', '68.41.25.65/32', false, 'France', '2011-11-10', 561.72);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (10, 'Scottie', 'Formilli', 'sformilli9@aol.com', 'Male', '101.241.191.228/32', false, 'Vietnam', '2003-01-04', 978.01);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (11, 'Lenci', 'Simecek', 'lsimeceka@cmu.edu', 'Agender', '252.190.171.190/32', false, 'Jamaica', '2001-09-23', 607.93);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (12, 'Faye', 'Elix', 'felixb@hatena.ne.jp', 'Female', '185.52.22.155/32', true, 'Peru', '2022-11-11', 972.79);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (13, 'Pail', 'Boxell', 'pboxellc@moonfruit.com', 'Genderqueer', '2.37.160.155/32', false, 'Indonesia', '2012-01-14', 73.47);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (14, 'Beverie', 'Jirousek', 'bjirousekd@arizona.edu', 'Female', '13.132.82.24/32', false, 'Indonesia', '2020-10-07', 950.04);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (15, 'Jecho', 'O''Garmen', 'jogarmene@woothemes.com', 'Male', '245.125.192.16/32', false, 'China', '2007-08-25', 257.19);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (16, 'Anni', 'Emmanueli', 'aemmanuelif@wikia.com', 'Female', '75.112.191.173/32', false, 'Saudi Arabia', '2015-08-11', 362.45);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (17, 'Evered', 'Marrow', 'emarrowg@tripadvisor.com', 'Male', '223.159.183.17/32', false, 'Niue', '2013-08-23', 418.18);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (18, 'Dell', 'Vondrak', 'dvondrakh@furl.net', 'Bigender', '83.89.160.155/32', true, 'Indonesia', '1999-03-06', 578.6);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (19, 'Renato', 'Swane', 'rswanei@scientificamerican.com', 'Genderfluid', '234.76.8.11/32', true, 'France', '2019-04-30', 80.81);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (20, 'Lucius', 'Edelmann', 'ledelmannj@bravesites.com', 'Male', '66.174.30.225/32', false, 'Brazil', '1999-06-22', 378.73);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (21, 'Gale', 'Norman', 'gnormank@skype.com', 'Female', '96.224.46.11/32', false, 'Russia', '2000-05-30', 152.93);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (22, 'Essie', 'McFfaden', 'emcffadenl@elpais.com', 'Female', '241.58.196.50/32', true, 'Greece', '2006-07-31', 75.77);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (23, 'Waly', 'Rogliero', 'wroglierom@berkeley.edu', 'Female', '122.90.196.231/32', true, 'Sweden', '2011-12-18', 147.69);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (24, 'Dillie', 'Iannazzi', 'diannazzin@biblegateway.com', 'Male', '112.79.17.198/32', true, 'Bangladesh', '2013-12-22', 699.62);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (25, 'Zolly', 'Wimmers', 'zwimmerso@hatena.ne.jp', 'Male', '123.12.225.114/32', false, 'Bosnia and Herzegovina', '2003-02-12', 217.18);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (26, 'Daryle', 'O''Shaughnessy', 'doshaughnessyp@com.com', 'Male', '204.107.16.207/32', false, 'Honduras', '2010-05-04', 983.03);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (27, 'Nick', 'Sings', 'nsingsq@boston.com', 'Male', '110.64.63.165/32', true, 'United States', '2011-03-17', 514.48);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (28, 'Ilsa', 'Huguenet', 'ihuguenetr@harvard.edu', 'Female', '147.1.198.181/32', false, 'China', '2014-05-11', 318.96);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (29, 'Grantham', 'Fayter', 'gfayters@desdev.cn', 'Male', '26.120.76.78/32', false, 'Sweden', '2009-10-02', 595.22);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (30, 'Bran', 'Longlands', 'blonglandst@tripod.com', 'Genderqueer', '14.92.3.30/32', false, 'France', '2016-07-10', 879.94);

CREATE TABLE CHARACTER_TYPES (char_type1 CHAR(8), char_type2 CHAR, varchar_type VARCHAR(255), text_type TEXT,
                              mediumtext_type MEDIUMTEXT, longtext_type LONGTEXT);

INSERT INTO CHARACTER_TYPES(char_type1, char_type2, varchar_type, text_type, mediumtext_type, longtext_type)
VALUES ('datagrok', 'A', 'Hello World', 'Hello Datagrok', 'HelloolleH', 'Datagrok');

CREATE TABLE DATE_TYPES(date_type date, time_type time, timestamp_type timestamp,
                        datetime_type datetime, year_type year);

INSERT INTO DATE_TYPES(date_type, time_type, timestamp_type, datetime_type, year_type)
VALUES ('9999-12-31', '18:59:59.000000', '1970-01-01 00:00:01.000000', '1000-01-01 00:00:00.000000', '2023');

create table dates_patterns (
    date DATE
);

insert into dates_patterns (date) values (CURDATE());
insert into dates_patterns (date) values (DATE_SUB(CURDATE(), INTERVAL 1 DAY)); -- yesterday

INSERT INTO dates_patterns(date) SELECT DATE(NOW() + INTERVAL (6 - WEEKDAY(NOW())) DAY) AS date
FROM dates_patterns WHERE (date = DATE(NOW() + INTERVAL (6 - WEEKDAY(NOW())) DAY)) HAVING COUNT(*) = 0; -- last day of week

insert into dates_patterns (date) values (DATE_SUB(CURDATE(), INTERVAL 150 DAY));
insert into dates_patterns (date) values ('2021-04-09');

CREATE TABLE FLOAT_TYPES(float_type float, double_type double, decimal_type decimal(65, 30));

INSERT INTO FLOAT_TYPES(float_type, double_type, decimal_type) VALUES (-3.402823466E+38, 1.7976931348623157E+308,
                                                                       999.9999);

INSERT INTO FLOAT_TYPES(float_type, double_type, decimal_type) VALUES (3.402823466E+38, -1.7976931348623157E+308,
                                                                       0.9999);

INSERT INTO FLOAT_TYPES(float_type, double_type, decimal_type) VALUES (-4.546544545, 4.457745745745457,
                                                                       23542363246234234234.46456456);

INSERT INTO FLOAT_TYPES(float_type, double_type, decimal_type) VALUES (-0.000124, 0.002,
                                                                       0.00001);

CREATE TABLE JSON_TYPE(json_type json);

INSERT INTO JSON_TYPE(json_type) VALUES ('{"key1": "value1", "key2": "value2"}');
INSERT INTO JSON_TYPE(json_type) VALUES ('{ "phones":[ {"type": "mobile", "phone": "001001"} , {"type": "fix", "phone": "002002"} ] }');
INSERT INTO JSON_TYPE(json_type) VALUES ('{"reading": 1.230e-5}');

CREATE TABLE INTEGER_TYPES (tinyint_type tinyint, smallint_type smallint, mediumint_type mediumint,
                            int_type int, bigint_type bigint);

INSERT INTO INTEGER_TYPES(tinyint_type, smallint_type, mediumint_type, int_type, bigint_type)
VALUES (12, 32000, 167772, 2147483647, 9223372036854775807);

INSERT INTO INTEGER_TYPES(tinyint_type, smallint_type, mediumint_type, int_type, bigint_type)
VALUES (0, 1212, -1000, -2147483648, -9223372036854775808);

CREATE TABLE GEOMETRY_TYPE (geometry_type GEOMETRY, point_type POINT);
INSERT INTO GEOMETRY_TYPE(geometry_type, point_type)VALUES (ST_GeomFromText('POLYGON((0 0,10 0,10 10,0 10,0 0),(5 5,7 5,7 7,5 7, 5 5))'),
                                                            ST_GeomFromText('POINT(1 1)'));

INSERT INTO GEOMETRY_TYPE(geometry_type, point_type)VALUES (ST_GeomFromText('GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 1,2 2,3 3,4 4))'),
                                                            ST_GeomFromText('POINT(1 0)'));


CREATE TABLE null_safety (binary_type BINARY(5), varbinary_type VARBINARY(10), bit_type1 bit(64), bit_type2 bit(1),
                          char_type1 CHAR(8), char_type2 CHAR, varchar_type VARCHAR(255), text_type TEXT,
                          mediumtext_type MEDIUMTEXT, longtext_type LONGTEXT, date_type date, time_type time, timestamp_type timestamp,
                          datetime_type datetime, year_type year, float_type float, double_type double, decimal_type decimal(65, 30),
                          json_type json, tinyint_type tinyint, smallint_type smallint, mediumint_type mediumint,
                          int_type int, bigint_type bigint, geometry_type GEOMETRY, point_type POINT);

INSERT INTO null_safety(binary_type, varbinary_type, bit_type1, bit_type2, char_type1, char_type2, varchar_type, text_type,
                        mediumtext_type, longtext_type, date_type, time_type, timestamp_type, datetime_type, year_type, float_type,
                        double_type, decimal_type, json_type, tinyint_type, smallint_type, mediumint_type, int_type, bigint_type,
                        geometry_type, point_type) VALUES (NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                                           NULL, NULL, NULL, NULL, NULL, NULL);
