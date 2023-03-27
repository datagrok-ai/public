CREATE EXTERNAL TABLE `mock_data`(
  `id` bigint,
  `first_name` varchar(50),
  `last_name` varchar(50),
  `email` varchar(50),
  `gender` varchar(50),
  `ip_address` varchar(50),
  `bool` boolean,
  `country` varchar(50),
  `date` date,
  `some_number` float)
ROW FORMAT SERDE
  'org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe'
STORED AS INPUTFORMAT
  'org.apache.hadoop.mapred.TextInputFormat'
OUTPUTFORMAT
  'org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat'
LOCATION
  's3://datagrok-provider-test/tables/mock-data'
TBLPROPERTIES (
  'transient_lastDdlTime'='1678111870')

insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (1, 'Burk', 'Kemery', 'bkemery0@businesswire.com', 'Male', '249.64.22.121/32', true, 'China', DATE'2017-09-20', 510.32);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (2, 'Nicholle', 'Karoly', 'nkaroly1@alexa.com', 'Female', '255.233.247.118/32', false, 'Poland', DATE'2014-02-27', 864.09);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (3, 'Orlando', 'Westgate', 'owestgate2@dedecms.com', 'Polygender', '75.0.252.254/32', false, 'Netherlands', DATE'2020-09-03', 822.7);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (4, 'Gothart', 'Cokayne', 'gcokayne3@plala.or.jp', 'Male', '196.83.12.163/32', true, 'Philippines', DATE'2001-01-31', 251.05);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (5, 'Mitchell', 'Haglington', 'mhaglington4@indiegogo.com', 'Male', '209.93.181.190/32', true, 'Poland', DATE'2020-10-09', 15.22);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (6, 'Jeromy', 'Twinn', 'jtwinn5@globo.com', 'Male', '25.13.2.132/32', true, 'Serbia', DATE'2014-10-04', 378.4);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (7, 'Joela', 'Cornau', 'jcornau6@imgur.com', 'Female', '195.47.88.236/32', false, 'Indonesia', DATE'2020-03-19', 349.11);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (8, 'Darren', 'Juares', 'djuares7@hexun.com', 'Male', '94.170.16.96/32', false, 'China', DATE'2011-04-09', 631.89);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (9, 'Marlie', 'Mayze', 'mmayze8@google.com.au', 'Female', '68.41.25.65/32', false, 'France', DATE'2011-11-10', 561.72);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (10, 'Scottie', 'Formilli', 'sformilli9@aol.com', 'Male', '101.241.191.228/32', false, 'Vietnam', DATE'2003-01-04', 978.01);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (11, 'Lenci', 'Simecek', 'lsimeceka@cmu.edu', 'Agender', '252.190.171.190/32', false, 'Jamaica', DATE'2001-09-23', 607.93);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (12, 'Faye', 'Elix', 'felixb@hatena.ne.jp', 'Female', '185.52.22.155/32', true, 'Peru', DATE'2022-11-11', 972.79);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (13, 'Pail', 'Boxell', 'pboxellc@moonfruit.com', 'Genderqueer', '2.37.160.155/32', false, 'Indonesia', DATE'2012-01-14', 73.47);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (14, 'Beverie', 'Jirousek', 'bjirousekd@arizona.edu', 'Female', '13.132.82.24/32', false, 'Indonesia', DATE'2020-10-07', 950.04);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (15, 'Jecho', 'O''Garmen', 'jogarmene@woothemes.com', 'Male', '245.125.192.16/32', false, 'China', DATE'2007-08-25', 257.19);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (16, 'Anni', 'Emmanueli', 'aemmanuelif@wikia.com', 'Female', '75.112.191.173/32', false, 'Saudi Arabia', DATE'2015-08-11', 362.45);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (17, 'Evered', 'Marrow', 'emarrowg@tripadvisor.com', 'Male', '223.159.183.17/32', false, 'Niue', DATE'2013-08-23', 418.18);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (18, 'Dell', 'Vondrak', 'dvondrakh@furl.net', 'Bigender', '83.89.160.155/32', true, 'Indonesia', DATE'1999-03-06', 578.6);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (19, 'Renato', 'Swane', 'rswanei@scientificamerican.com', 'Genderfluid', '234.76.8.11/32', true, 'France', DATE'2019-04-30', 80.81);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (20, 'Lucius', 'Edelmann', 'ledelmannj@bravesites.com', 'Male', '66.174.30.225/32', false, 'Brazil', DATE'1999-06-22', 378.73);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (21, 'Gale', 'Norman', 'gnormank@skype.com', 'Female', '96.224.46.11/32', false, 'Russia', DATE'2000-05-30', 152.93);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (22, 'Essie', 'McFfaden', 'emcffadenl@elpais.com', 'Female', '241.58.196.50/32', true, 'Greece', DATE'2006-07-31', 75.77);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (23, 'Waly', 'Rogliero', 'wroglierom@berkeley.edu', 'Female', '122.90.196.231/32', true, 'Sweden', DATE'2011-12-18', 147.69);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (24, 'Dillie', 'Iannazzi', 'diannazzin@biblegateway.com', 'Male', '112.79.17.198/32', true, 'Bangladesh', DATE'2013-12-22', 699.62);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (25, 'Zolly', 'Wimmers', 'zwimmerso@hatena.ne.jp', 'Male', '123.12.225.114/32', false, 'Bosnia and Herzegovina', DATE'2003-02-12', 217.18);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (26, 'Daryle', 'O''Shaughnessy', 'doshaughnessyp@com.com', 'Male', '204.107.16.207/32', false, 'Honduras', DATE'2010-05-04', 983.03);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (27, 'Nick', 'Sings', 'nsingsq@boston.com', 'Male', '110.64.63.165/32', true, 'United States', DATE'2011-03-17', 514.48);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (28, 'Ilsa', 'Huguenet', 'ihuguenetr@harvard.edu', 'Female', '147.1.198.181/32', false, 'China', DATE'2014-05-11', 318.96);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (29, 'Grantham', 'Fayter', 'gfayters@desdev.cn', 'Male', '26.120.76.78/32', false, 'Sweden', DATE'2009-10-02', 595.22);
insert into mock_data (id, first_name, last_name, email, gender, ip_address, bool, country, date, some_number) values (30, 'Bran', 'Longlands', 'blonglandst@tripod.com', 'Genderqueer', '14.92.3.30/32', false, 'France', DATE'2016-07-10', 879.94);
