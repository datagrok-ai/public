create table mock_data (
                           id UInt64,
                           first_name String,
                           last_name String,
                           email String,
                           gender String,
                           ip_address String,
                           bool Bool,
                           country String,
                           date Date,
                           some_number Float64
)
    ENGINE = MergeTree()
    ORDER BY tuple();

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

create table dates_patterns (
    date Date
) ENGINE = MergeTree()
    ORDER BY tuple();

insert into dates_patterns (date) values (today());
insert into dates_patterns (date) values (yesterday());

INSERT INTO dates_patterns(date) SELECT (today() + (7 - toDayOfWeek(today())))
FROM dates_patterns WHERE (date = (today() + (7 - toDayOfWeek(today())))) HAVING COUNT(*) = 0;

insert into dates_patterns (date) values (today() - 150);
insert into dates_patterns (date) values ('2021-04-09');

CREATE TABLE uuid_type (uuid UUID)
    ENGINE = MergeTree()
    ORDER BY tuple();
INSERT INTO uuid_type(uuid) VALUES ('417ddc5d-e556-4d27-95dd-a34d84e46a50');

CREATE TABLE SIGNED_INTEGER_TYPES (int8_type Int8, int16_type Int16, int32_type Int32, int64_type Int64,
                                   int128_type Int128, int256_type Int256) ENGINE = MergeTree()
    ORDER BY tuple();

INSERT INTO SIGNED_INTEGER_TYPES(int8_type, int16_type, int32_type, int64_type, int128_type, int256_type)
VALUES (-128, -32768, -2147483648, -9223372036854775808, -170141183460469231731687303715884105728,
        -57896044618658097711785492504343953926634992332820282019728792003956564819968);

INSERT INTO SIGNED_INTEGER_TYPES(int8_type, int16_type, int32_type, int64_type, int128_type, int256_type)
VALUES (127, 32767, 2147483647, 9223372036854775807, 170141183460469231731687303715884105727,
        57896044618658097711785492504343953926634992332820282019728792003956564819967);

CREATE TABLE UNSIGNED_INTEGER_TYPES (uint8_type UInt8, uint16_type UInt16, uint32_type UInt32, uint64_type UInt64,
                                     uint128_type UInt128, uint256_type UInt256) ENGINE = MergeTree()
    ORDER BY tuple();

INSERT INTO UNSIGNED_INTEGER_TYPES(uint8_type, uint16_type, uint32_type, uint64_type, uint128_type, uint256_type)
VALUES (255, 65535, 4294967295, 18446744073709551615, 340282366920938463463374607431768211455,
        115792089237316195423570985008687907853269984665640564039457584007913129639935);

CREATE TABLE FLOAT_TYPES(float32_type Float32, float64_type Float64, decimal32_type Decimal32(4),
                         decimal64_type Decimal64(9), decimal128_type Decimal128(18), decimal256_type Decimal256(26)) ENGINE = MergeTree()
    ORDER BY tuple();

INSERT INTO FLOAT_TYPES(float32_type, float64_type, decimal32_type, decimal64_type, decimal128_type, decimal256_type)
VALUES (0.5 / 0, 1.79E+308, 999.9999, 999999999.999999, 211111111111.22989999,
        999999999999999999.99999999999999 );

INSERT INTO FLOAT_TYPES(float32_type, float64_type, decimal32_type, decimal64_type, decimal128_type, decimal256_type)
VALUES (-0.5 / 0,  0 / 0, -999.9999, -999999999.999999, -21111111111.122988999,
        -999999999999999999999.999999999999 );

CREATE TABLE ARRAY_TYPE(array_type Array(UInt32), array_array_type Array(Array(UInt32)))
    ENGINE = MergeTree()
    ORDER BY tuple();
INSERT INTO ARRAY_TYPE (array_type, array_array_type) VALUES ([1, 2, 3, 4], [[24], [421, 12, 4], []]);
