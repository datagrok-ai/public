CREATE TABLE sal_emp (
                         name            text,
                         pay_by_quarter  integer[],
                         schedule        text[][]
);

INSERT INTO sal_emp
VALUES ('Bill',
        '{10000, 10000, 10000, 10000}',
        '{{"meeting", "lunch"}, {"training", "presentation"}}');

INSERT INTO sal_emp
VALUES ('Carol',
        '{20000, 25000, 25000, 25000}',
        '{{"breakfast", "consulting"}, {"meeting", "lunch"}}');

create table MOCK_DATA (
                           id bigint,
                           first_name VARCHAR(50),
                           last_name VARCHAR(50),
                           email VARCHAR(50),
                           gender VARCHAR(50),
                           ip_address cidr,
                           bool boolean,
                           country VARCHAR(50),
                           date DATE,
                           some_number numeric(5,2)
);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (1, 'Burk', 'Kemery', 'bkemery0@businesswire.com', 'Male', '249.64.22.121', true, 'China', '2017-09-20', 510.32);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (2, 'Nicholle', 'Karoly', 'nkaroly1@alexa.com', 'Female', '255.233.247.118', false, 'Poland', '2014-02-27', 864.09);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (3, 'Orlando', 'Westgate', 'owestgate2@dedecms.com', 'Polygender', '75.0.252.254', false, 'Netherlands', '2020-09-03', 822.7);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (4, 'Gothart', 'Cokayne', 'gcokayne3@plala.or.jp', 'Male', '196.83.12.163', true, 'Philippines', '2001-01-31', 251.05);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (5, 'Mitchell', 'Haglington', 'mhaglington4@indiegogo.com', 'Male', '209.93.181.190', true, 'Poland', '2020-10-09', 15.22);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (6, 'Jeromy', 'Twinn', 'jtwinn5@globo.com', 'Male', '25.13.2.132', true, 'Serbia', '2014-10-04', 378.4);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (7, 'Joela', 'Cornau', 'jcornau6@imgur.com', 'Female', '195.47.88.236', false, 'Indonesia', '2020-03-19', 349.11);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (8, 'Darren', 'Juares', 'djuares7@hexun.com', 'Male', '94.170.16.96', false, 'China', '2011-04-09', 631.89);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (9, 'Marlie', 'Mayze', 'mmayze8@google.com.au', 'Female', '68.41.25.65', false, 'France', '2011-11-10', 561.72);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (10, 'Scottie', 'Formilli', 'sformilli9@aol.com', 'Male', '101.241.191.228', false, 'Vietnam', '2003-01-04', 978.01);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (11, 'Lenci', 'Simecek', 'lsimeceka@cmu.edu', 'Agender', '252.190.171.190', false, 'Jamaica', '2001-09-23', 607.93);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (12, 'Faye', 'Elix', 'felixb@hatena.ne.jp', 'Female', '185.52.22.155', true, 'Peru', '2022-11-11', 972.79);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (13, 'Pail', 'Boxell', 'pboxellc@moonfruit.com', 'Genderqueer', '2.37.160.155', false, 'Indonesia', '2012-01-14', 73.47);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (14, 'Beverie', 'Jirousek', 'bjirousekd@arizona.edu', 'Female', '13.132.82.24', false, 'Indonesia', '2020-10-07', 950.04);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (15, 'Jecho', 'O''Garmen', 'jogarmene@woothemes.com', 'Male', '245.125.192.16', false, 'China', '2007-08-25', 257.19);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (16, 'Anni', 'Emmanueli', 'aemmanuelif@wikia.com', 'Female', '75.112.191.173', false, 'Saudi Arabia', '2015-08-11', 362.45);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (17, 'Evered', 'Marrow', 'emarrowg@tripadvisor.com', 'Male', '223.159.183.17', false, 'Niue', '2013-08-23', 418.18);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (18, 'Dell', 'Vondrak', 'dvondrakh@furl.net', 'Bigender', '83.89.160.155', true, 'Indonesia', '1999-03-06', 578.6);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (19, 'Renato', 'Swane', 'rswanei@scientificamerican.com', 'Genderfluid', '234.76.8.11', true, 'France', '2019-04-30', 80.81);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (20, 'Lucius', 'Edelmann', 'ledelmannj@bravesites.com', 'Male', '66.174.30.225', false, 'Brazil', '1999-06-22', 378.73);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (21, 'Gale', 'Norman', 'gnormank@skype.com', 'Female', '96.224.46.11', false, 'Russia', '2000-05-30', 152.93);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (22, 'Essie', 'McFfaden', 'emcffadenl@elpais.com', 'Female', '241.58.196.50', true, 'Greece', '2006-07-31', 75.77);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (23, 'Waly', 'Rogliero', 'wroglierom@berkeley.edu', 'Female', '122.90.196.231', true, 'Sweden', '2011-12-18', 147.69);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (24, 'Dillie', 'Iannazzi', 'diannazzin@biblegateway.com', 'Male', '112.79.17.198', true, 'Bangladesh', '2013-12-22', 699.62);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (25, 'Zolly', 'Wimmers', 'zwimmerso@hatena.ne.jp', 'Male', '123.12.225.114', false, 'Bosnia and Herzegovina', '2003-02-12', 217.18);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (26, 'Daryle', 'O''Shaughnessy', 'doshaughnessyp@com.com', 'Male', '204.107.16.207', false, 'Honduras', '2010-05-04', 983.03);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (27, 'Nick', 'Sings', 'nsingsq@boston.com', 'Male', '110.64.63.165', true, 'United States', '2011-03-17', 514.48);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (28, 'Ilsa', 'Huguenet', 'ihuguenetr@harvard.edu', 'Female', '147.1.198.181', false, 'China', '2014-05-11', 318.96);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (29, 'Grantham', 'Fayter', 'gfayters@desdev.cn', 'Male', '26.120.76.78', false, 'Sweden', '2009-10-02', 595.22);
insert into MOCK_DATA (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (30, 'Bran', 'Longlands', 'blonglandst@tripod.com', 'Genderqueer', '14.92.3.30', false, 'France', '2016-07-10', 879.94);


CREATE TABLE test1 (a BIT(3), b BIT VARYING(5));
INSERT INTO test1 VALUES (B'101', B'0011');
INSERT INTO test1 VALUES (B'001', B'101');

CREATE TABLE test2 (c BIT(1), d BIT VARYING(5));
INSERT INTO test2 VALUES (B'1', B'101');
INSERT INTO test2 VALUES (B'0', B'0');

CREATE TABLE BYTEA_DATA (data BYTEA);

INSERT INTO BYTEA_DATA VALUES (pg_read_binary_file('/etc/file.txt'));

CREATE TYPE inventory_item AS (
    name            text,
    supplier_id     integer,
    price           numeric
    );

CREATE TABLE on_hand (
                         item      inventory_item,
                         count     integer
);

INSERT INTO on_hand VALUES (ROW('fuzzy dice', 42, 1.99), 1000);

CREATE TABLE DATES(date DATE, time TIME, stamp TIMESTAMP, interval INTERVAL);

INSERT INTO DATES("date", "time", stamp, "interval") VALUES ('1999-01-08', '04:05:06.789', '1999-01-08 04:05:06',
                                                             '1 year 5 months 5 days');

create table DATES_PATTERNS (
    date DATE
);

insert into DATES_PATTERNS (date) values (current_date);
insert into DATES_PATTERNS (date) values (current_date - 1); -- yesterday
insert into DATES_PATTERNS (date)
select (date_trunc('week', current_date::timestamp) + '6 days'::interval)::date
    where not exists (select date from dates_patterns where date = (date_trunc('week', current_date::timestamp) + '6 days'::interval)::date); -- last day of this week
insert into DATES_PATTERNS (date) values (current_date - 150);
insert into DATES_PATTERNS (date) values ('2021-04-09');

CREATE TABLE JSONB_DATA (
    data JSONB
);

INSERT INTO JSONB_DATA(data)
VALUES ('{ "phones":[ {"type": "mobile", "phone": "001001"} , {"type": "fix", "phone": "002002"} ] }');

INSERT INTO JSONB_DATA(data) VALUES ('{"bar": "baz", "balance": 7.77, "active":false}');

INSERT INTO JSONB_DATA(data) VALUES ('{"reading": 1.230e-5}');


CREATE TABLE null_safety (array_type integer[], varchar_type varchar(50), cidr_type cidr, bool_type bool,
                          bit_type bit(3), bit_var_type bit varying(5), bytea_type bytea, date_type date, time_type time, stamp_type timestamp,
                          interval_type interval, json_type jsonb, numeric_type numeric, double_type double precision, real_type real,
                          bigint_type bigint, uuid_type uuid, xml_type xml);


CREATE TABLE NUMERIC_DATA (
    data NUMERIC
);

CREATE TABLE DOUBLES (
    data DOUBLE PRECISION
);

CREATE TABLE REALS (
    data REAL
);

CREATE TABLE BIGINT_DATA (
    data BIGINT
);

CREATE TABLE NUMERIC_DATA_PRECISION (
    data NUMERIC(5)
);

CREATE TABLE NUMERIC_DATA_PRECISION_SCALE (
    data NUMERIC(4, 2)
);

INSERT INTO NUMERIC_DATA_PRECISION (data) VALUES (2.2221);
INSERT INTO NUMERIC_DATA_PRECISION (data) VALUES (2222.2);

INSERT INTO NUMERIC_DATA_PRECISION_SCALE (data) VALUES (22.22);
INSERT INTO NUMERIC_DATA_PRECISION_SCALE (data) VALUES (55.55);

INSERT INTO BIGINT_DATA (data) VALUES (-9223372036854775808); -- min value for big int
INSERT INTO BIGINT_DATA (data) VALUES (9223372036854775807); -- max value for big int

INSERT INTO REALS (data) VALUES (1E-37); -- min value for real
INSERT INTO REALS (data) VALUES (1E+37); -- max value for float

INSERT INTO DOUBLES (data) VALUES ('-Infinity'::float); --edge case
INSERT INTO DOUBLES (data) VALUES ('+Infinity'::float); ----edge case
INSERT INTO DOUBLES (data) VALUES (1E-307); -- close zero value
INSERT INTO DOUBLES (data) VALUES (1E+308); -- max value for double

INSERT INTO NUMERIC_DATA (data) VALUES ('NaN'::numeric); --edge case
INSERT INTO NUMERIC_DATA (data) VALUES (2.22222222222);
INSERT INTO NUMERIC_DATA (data) VALUES (222222222222222222.2);

CREATE TABLE operators (id bigint, json_data jsonb, path_data path, circle_data circle);

INSERT INTO operators (id, json_data, path_data, circle_data) VALUES (1, '{"reading": 1.230e-5}', path('((0,0),(1,0))'), circle('((0,0),10)'));

INSERT INTO operators (id, json_data, path_data, circle_data) VALUES (2, '{"bar": "baz", "balance": 7.77, "active":false}', path(
        '(15.878137629895164,47.08306448089695),
         (15.56169808311181,47.219041634920686),
         (15.267442604782124,47.4201665137259),
         (15.092631384557304,47.71366328136526),
         (15.234428926980286,47.95865145177352)'
    ), circle('((0,0),15)'));

CREATE TABLE CARS_SMALL (
                            id SMALLSERIAL PRIMARY KEY,
                            brand VARCHAR NOT NULL
);

CREATE TABLE CARS (
                      id SERIAL PRIMARY KEY,
                      brand VARCHAR NOT NULL
);

CREATE TABLE CARS_BIG (
                          id BIGSERIAL PRIMARY KEY,
                          brand VARCHAR NOT NULL
);

INSERT INTO CARS_SMALL(brand)
VALUES('Fiat');

INSERT INTO CARS_SMALL(id, brand)
VALUES(DEFAULT, 'Honda');

INSERT INTO CARS(brand)
VALUES('Fiat');

INSERT INTO CARS(id, brand)
VALUES(DEFAULT, 'Honda');

INSERT INTO CARS_BIG(brand)
VALUES('Fiat');

INSERT INTO CARS_BIG(id, brand)
VALUES(DEFAULT, 'Honda');

CREATE TABLE UUID_DATA (
    id uuid
);

INSERT INTO UUID_DATA (id) VALUES ('6ecd8c99-4036-403d-bf84-cf8400f67836');
INSERT INTO UUID_DATA (id) VALUES ('c81d4e2e-bcf2-11e6-869b-7df92533d2db');
INSERT INTO UUID_DATA (id) VALUES ('237e9877-e79b-12d4-a765-321741963000');

CREATE TABLE XML_DATA (
    data XML
);

INSERT INTO XML_DATA (data) VALUES (XML('<foo>Hello World!</foo>'));
INSERT INTO XML_DATA (data) VALUES (XML('abc<foo>bar</foo><bar>foo</bar>'));
INSERT INTO XML_DATA (data) VALUES (XML('<?xml version="1.0"?><book><title>Manual</title><chapter>...</chapter></book>'));
