CREATE KEYSPACE datagrok WITH REPLICATION = {
'class' : 'SimpleStrategy',
'replication_factor' : 1
};

CREATE TABLE datagrok.mock_data (
id int,
first_name text,
last_name text,
email text,
gender text,
ip_address text,
bool boolean,
country text,
date date,
some_number double,
PRIMARY KEY (id));

CREATE CUSTOM INDEX mock_data_email_idx  ON datagrok.mock_data ("email") USING 'org.apache.cassandra.index.sasi.SASIIndex' WITH OPTIONS = {'mode': 'CONTAINS', 'analyzer_class': 'org.apache.cassandra.index.sasi.analyzer.NonTokenizingAnalyzer', 'case_sensitive': false};
CREATE CUSTOM INDEX mock_data_first_name_idx  ON datagrok.mock_data ("first_name") USING 'org.apache.cassandra.index.sasi.SASIIndex' WITH OPTIONS = {'mode': 'CONTAINS', 'analyzer_class': 'org.apache.cassandra.index.sasi.analyzer.NonTokenizingAnalyzer', 'case_sensitive': false};

insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (1, 'Burk', 'Kemery', 'bkemery0@businesswire.com', 'Male', '249.64.22.121/32', true, 'China', '2017-09-20', 510.32);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (2, 'Nicholle', 'Karoly', 'nkaroly1@alexa.com', 'Female', '255.233.247.118/32', false, 'Poland', '2014-02-27', 864.09);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (3, 'Orlando', 'Westgate', 'owestgate2@dedecms.com', 'Polygender', '75.0.252.254/32', false, 'Netherlands', '2020-09-03', 822.7);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (4, 'Gothart', 'Cokayne', 'gcokayne3@plala.or.jp', 'Male', '196.83.12.163/32', true, 'Philippines', '2001-01-31', 251.05);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (5, 'Mitchell', 'Haglington', 'mhaglington4@indiegogo.com', 'Male', '209.93.181.190/32', true, 'Poland', '2020-10-09', 15.22);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (6, 'Jeromy', 'Twinn', 'jtwinn5@globo.com', 'Male', '25.13.2.132/32', true, 'Serbia', '2014-10-04', 378.4);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (7, 'Joela', 'Cornau', 'jcornau6@imgur.com', 'Female', '195.47.88.236/32', false, 'Indonesia', '2020-03-19', 349.11);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (8, 'Darren', 'Juares', 'djuares7@hexun.com', 'Male', '94.170.16.96/32', false, 'China', '2011-04-09', 631.89);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (9, 'Marlie', 'Mayze', 'mmayze8@google.com.au', 'Female', '68.41.25.65/32', false, 'France', '2011-11-10', 561.72);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (10, 'Scottie', 'Formilli', 'sformilli9@aol.com', 'Male', '101.241.191.228/32', false, 'Vietnam', '2003-01-04', 978.01);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (11, 'Lenci', 'Simecek', 'lsimeceka@cmu.edu', 'Agender', '252.190.171.190/32', false, 'Jamaica', '2001-09-23', 607.93);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (12, 'Faye', 'Elix', 'felixb@hatena.ne.jp', 'Female', '185.52.22.155/32', true, 'Peru', '2022-11-11', 972.79);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (13, 'Pail', 'Boxell', 'pboxellc@moonfruit.com', 'Genderqueer', '2.37.160.155/32', false, 'Indonesia', '2012-01-14', 73.47);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (14, 'Beverie', 'Jirousek', 'bjirousekd@arizona.edu', 'Female', '13.132.82.24/32', false, 'Indonesia', '2020-10-07', 950.04);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (15, 'Jecho', 'O''Garmen', 'jogarmene@woothemes.com', 'Male', '245.125.192.16/32', false, 'China', '2007-08-25', 257.19);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (16, 'Anni', 'Emmanueli', 'aemmanuelif@wikia.com', 'Female', '75.112.191.173/32', false, 'Saudi Arabia', '2015-08-11', 362.45);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (17, 'Evered', 'Marrow', 'emarrowg@tripadvisor.com', 'Male', '223.159.183.17/32', false, 'Niue', '2013-08-23', 418.18);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (18, 'Dell', 'Vondrak', 'dvondrakh@furl.net', 'Bigender', '83.89.160.155/32', true, 'Indonesia', '1999-03-06', 578.6);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (19, 'Renato', 'Swane', 'rswanei@scientificamerican.com', 'Genderfluid', '234.76.8.11/32', true, 'France', '2019-04-30', 80.81);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (20, 'Lucius', 'Edelmann', 'ledelmannj@bravesites.com', 'Male', '66.174.30.225/32', false, 'Brazil', '1999-06-22', 378.73);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (21, 'Gale', 'Norman', 'gnormank@skype.com', 'Female', '96.224.46.11/32', false, 'Russia', '2000-05-30', 152.93);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (22, 'Essie', 'McFfaden', 'emcffadenl@elpais.com', 'Female', '241.58.196.50/32', true, 'Greece', '2006-07-31', 75.77);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (23, 'Waly', 'Rogliero', 'wroglierom@berkeley.edu', 'Female', '122.90.196.231/32', true, 'Sweden', '2011-12-18', 147.69);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (24, 'Dillie', 'Iannazzi', 'diannazzin@biblegateway.com', 'Male', '112.79.17.198/32', true, 'Bangladesh', '2013-12-22', 699.62);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (25, 'Zolly', 'Wimmers', 'zwimmerso@hatena.ne.jp', 'Male', '123.12.225.114/32', false, 'Bosnia and Herzegovina', '2003-02-12', 217.18);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (26, 'Daryle', 'O''Shaughnessy', 'doshaughnessyp@com.com', 'Male', '204.107.16.207/32', false, 'Honduras', '2010-05-04', 983.03);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (27, 'Nick', 'Sings', 'nsingsq@boston.com', 'Male', '110.64.63.165/32', true, 'United States', '2011-03-17', 514.48);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (28, 'Ilsa', 'Huguenet', 'ihuguenetr@harvard.edu', 'Female', '147.1.198.181/32', false, 'China', '2014-05-11', 318.96);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (29, 'Grantham', 'Fayter', 'gfayters@desdev.cn', 'Male', '26.120.76.78/32', false, 'Sweden', '2009-10-02', 595.22);
insert into datagrok.mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (30, 'Bran', 'Longlands', 'blonglandst@tripod.com', 'Genderqueer', '14.92.3.30/32', false, 'France', '2016-07-10', 879.94);

CREATE TABLE datagrok.date_types (
id bigint,
date_type date,
time_type time,
stamp_type timestamp,
duration_type duration,
PRIMARY KEY (id));

insert into datagrok.date_types (id, date_type, time_type, stamp_type, duration_type) values(1, '2011-02-03', '08:12:54.111',
'2011-02-03 04:05:01+0100', 89h4m48s);
insert into datagrok.date_types (id, date_type, time_type, stamp_type, duration_type) values(2, '2023-04-27', '12:50:54.333',
'2011-02-03 04:05:08+0300', 9m51s);

CREATE TABLE datagrok.dates_patterns(
id bigint,
date date,
PRIMARY KEY (id));

insert into datagrok.dates_patterns (id, date) values (1, currentDate());
insert into datagrok.dates_patterns (id, date) values(2, currentDate() - 1d);

CREATE OR REPLACE FUNCTION datagrok.lastWeekDay()
CALLED ON NULL INPUT
RETURNS DATE
LANGUAGE java
AS '
  java.time.LocalDate now =  java.time.LocalDate.now();
  int dayOfWeek = now.getDayOfWeek().getValue();
  return com.datastax.driver.core.LocalDate.fromDaysSinceEpoch(Math.toIntExact(now.plusDays(7 - dayOfWeek).toEpochDay()));
';

insert into datagrok.dates_patterns (id, date) values(3, lastWeekDay()) if not exists;
insert into datagrok.dates_patterns (id, date) values(4, currentDate() - 150d);
insert into datagrok.dates_patterns (id, date) values(5, '2021-04-09');

CREATE TABLE datagrok.integers(
tinyint_type tinyint,
smallint_type smallint,
int_type int,
bigint_type bigint,
varint_type varint,
PRIMARY KEY (tinyint_type));

insert into datagrok.integers(tinyint_type, smallint_type, int_type, bigint_type, varint_type) values(127, 32767,
2147483647, 999999999999999999, 999999999999999999);
insert into datagrok.integers(tinyint_type, smallint_type, int_type, bigint_type, varint_type) values(-128, -32768,
-2147483648, -999999999999999999, -9999999999999999999);

CREATE TABLE datagrok.float_types(
decimal_type decimal,
float_type float,
double_type double,
PRIMARY KEY (decimal_type));

insert into datagrok.float_types(decimal_type, float_type, double_type) values (123003.33211, 0.000123, 3.13457599923384753929348);
insert into datagrok.float_types(decimal_type, float_type, double_type) values (-20.01001, -9990.99221, -3.13457599923384753929348);

CREATE TYPE datagrok.phone (
country_code int,
number text,
);

CREATE TABLE datagrok.complex_types(
id bigint,
map_type map<text, text>,
set_type set<int>,
list_type frozen<list<varchar>>,
udt_type phone,
tuple_type tuple<int, text>,
PRIMARY KEY (id));

insert into datagrok.complex_types(id, map_type, set_type, list_type, udt_type, tuple_type) values(1, { 'fruit' : 'Apple', 'band' : 'Beatles' },
{ 3, 15, 16 }, ['hello', 'world'], {country_code: 380, number: '+3800900909'}, (3, 'hours'));
