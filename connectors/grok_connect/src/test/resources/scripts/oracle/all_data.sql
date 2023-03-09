-- ORACLE CONTAINER REQUIRES ADDITIONAL WORK TO
-- EXECUTE NEW SCRIPTS DURING RUN, DECIDED TO INJECT ALL DATA AT THE BEGINNING

CREATE TABLE MOCK_DATA (
                           id NUMBER(18),
                           first_name VARCHAR2(50),
                           last_name VARCHAR2(50),
                           email VARCHAR2(50),
                           gender VARCHAR2(50),
                           ip_address VARCHAR2(50),
                           country VARCHAR2(50),
                           dat DATE,
                           some_number NUMBER(5,2)
);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (1, 'Burk', 'Kemery', 'bkemery0@businesswire.com', 'Male', '249.64.22.121/32',  'China', TO_DATE('2017-09-20', 'YYYY-MM-DD'), 510.32);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (2, 'Nicholle', 'Karoly', 'nkaroly1@alexa.com', 'Female', '255.233.247.118/32',  'Poland', TO_DATE('2014-02-27', 'YYYY-MM-DD'), 864.09);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (3, 'Orlando', 'Westgate', 'owestgate2@dedecms.com', 'Polygender', '75.0.252.254/32',  'Netherlands', TO_DATE('2020-09-03', 'YYYY-MM-DD'), 822.7);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (4, 'Gothart', 'Cokayne', 'gcokayne3@plala.or.jp', 'Male', '196.83.12.163/32',  'Philippines', TO_DATE('2001-01-31', 'YYYY-MM-DD'), 251.05);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (5, 'Mitchell', 'Haglington', 'mhaglington4@indiegogo.com', 'Male', '209.93.181.190/32',  'Poland', TO_DATE('2020-10-09', 'YYYY-MM-DD'), 15.22);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (6, 'Jeromy', 'Twinn', 'jtwinn5@globo.com', 'Male', '25.13.2.132/32',  'Serbia', TO_DATE('2014-10-04', 'YYYY-MM-DD'), 378.4);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (7, 'Joela', 'Cornau', 'jcornau6@imgur.com', 'Female', '195.47.88.236/32',  'Indonesia', TO_DATE('2020-03-19', 'YYYY-MM-DD'), 349.11);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (8, 'Darren', 'Juares', 'djuares7@hexun.com', 'Male', '94.170.16.96/32',  'China', TO_DATE('2011-04-09', 'YYYY-MM-DD'), 631.89);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (9, 'Marlie', 'Mayze', 'mmayze8@google.com.au', 'Female', '68.41.25.65/32',  'France', TO_DATE('2011-11-10', 'YYYY-MM-DD'), 561.72);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (10, 'Scottie', 'Formilli', 'sformilli9@aol.com', 'Male', '101.241.191.228/32',  'Vietnam', TO_DATE('2003-01-04', 'YYYY-MM-DD'), 978.01);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (11, 'Lenci', 'Simecek', 'lsimeceka@cmu.edu', 'Agender', '252.190.171.190/32',  'Jamaica', TO_DATE('2001-09-23', 'YYYY-MM-DD'), 607.93);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (12, 'Faye', 'Elix', 'felixb@hatena.ne.jp', 'Female', '185.52.22.155/32',  'Peru', TO_DATE('2022-11-11', 'YYYY-MM-DD'), 972.79);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (13, 'Pail', 'Boxell', 'pboxellc@moonfruit.com', 'Genderqueer', '2.37.160.155/32',  'Indonesia', TO_DATE('2012-01-14', 'YYYY-MM-DD'), 73.47);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (14, 'Beverie', 'Jirousek', 'bjirousekd@arizona.edu', 'Female', '13.132.82.24/32',  'Indonesia', TO_DATE('2020-10-07', 'YYYY-MM-DD'), 950.04);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (15, 'Jecho', 'O''Garmen', 'jogarmene@woothemes.com', 'Male', '245.125.192.16/32',  'China', TO_DATE('2007-08-25', 'YYYY-MM-DD'), 257.19);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (16, 'Anni', 'Emmanueli', 'aemmanuelif@wikia.com', 'Female', '75.112.191.173/32',  'Saudi Arabia', TO_DATE('2015-08-11', 'YYYY-MM-DD'), 362.45);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (17, 'Evered', 'Marrow', 'emarrowg@tripadvisor.com', 'Male', '223.159.183.17/32',  'Niue', TO_DATE('2013-08-23', 'YYYY-MM-DD'), 418.18);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (18, 'Dell', 'Vondrak', 'dvondrakh@furl.net', 'Bigender', '83.89.160.155/32',  'Indonesia', TO_DATE('1999-03-06', 'YYYY-MM-DD'), 578.6);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (19, 'Renato', 'Swane', 'rswanei@scientificamerican.com', 'Genderfluid', '234.76.8.11/32',  'France', TO_DATE('2019-04-30', 'YYYY-MM-DD'), 80.81);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (20, 'Lucius', 'Edelmann', 'ledelmannj@bravesites.com', 'Male', '66.174.30.225/32',  'Brazil', TO_DATE('1999-06-22', 'YYYY-MM-DD'), 378.73);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (21, 'Gale', 'Norman', 'gnormank@skype.com', 'Female', '96.224.46.11/32',  'Russia', TO_DATE('2000-05-30', 'YYYY-MM-DD'), 152.93);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (22, 'Essie', 'McFfaden', 'emcffadenl@elpais.com', 'Female', '241.58.196.50/32',  'Greece', TO_DATE('2006-07-31', 'YYYY-MM-DD'), 75.77);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (23, 'Waly', 'Rogliero', 'wroglierom@berkeley.edu', 'Female', '122.90.196.231/32',  'Sweden', TO_DATE('2011-12-18', 'YYYY-MM-DD'), 147.69);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (24, 'Dillie', 'Iannazzi', 'diannazzin@biblegateway.com', 'Male', '112.79.17.198/32',  'Bangladesh', TO_DATE('2013-12-22', 'YYYY-MM-DD'), 699.62);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (25, 'Zolly', 'Wimmers', 'zwimmerso@hatena.ne.jp', 'Male', '123.12.225.114/32',  'Bosnia and Herzegovina', TO_DATE('2003-02-12', 'YYYY-MM-DD'), 217.18);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (26, 'Daryle', 'O''Shaughnessy', 'doshaughnessyp@com.com', 'Male', '204.107.16.207/32',  'Honduras', TO_DATE('2010-05-04', 'YYYY-MM-DD'), 983.03);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (27, 'Nick', 'Sings', 'nsingsq@boston.com', 'Male', '110.64.63.165/32',  'United States', TO_DATE('2011-03-17', 'YYYY-MM-DD'), 514.48);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (28, 'Ilsa', 'Huguenet', 'ihuguenetr@harvard.edu', 'Female', '147.1.198.181/32',  'China', TO_DATE('2014-05-11', 'YYYY-MM-DD'), 318.96);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (29, 'Grantham', 'Fayter', 'gfayters@desdev.cn', 'Male', '26.120.76.78/32',  'Sweden', TO_DATE('2009-10-02', 'YYYY-MM-DD'), 595.22);
INSERT INTO MOCK_DATA (id, first_name, last_name, email, gender, ip_address, country, dat, some_number) VALUES (30, 'Bran', 'Longlands', 'blonglandst@tripod.com', 'Genderqueer', '14.92.3.30/32',  'France', TO_DATE('2016-07-10', 'YYYY-MM-DD'), 879.94);



CREATE TABLE character_type (
                                ch CHAR(10),
                                varch VARCHAR2(10),
                                nch NCHAR(10),
                                nvarch NVARCHAR2(50)
);

INSERT INTO character_type(ch, varch, nch, nvarch) VALUES ('Hello', 'World', 'Datagrok', 'Groking');


CREATE TABLE dates_patterns (
    dat DATE
);

INSERT INTO dates_patterns(dat) VALUES (SYSDATE); --TODAY
INSERT INTO dates_patterns(dat) VALUES (SYSDATE - 1); --YESTERDAY
INSERT INTO dates_patterns(dat) SELECT TRUNC(SYSDATE + 6, 'DAY') FROM dual
WHERE NOT EXISTS (SELECT * FROM dates_patterns WHERE to_date(dat) = to_date(TRUNC(SYSDATE + 6, 'DAY')));
INSERT INTO dates_patterns(dat) VALUES (SYSDATE - 150);
INSERT INTO dates_patterns(dat) VALUES (TO_DATE('2021-04-09', 'YYYY-MM-DD'));

CREATE TABLE dates_type (
                            dat DATE,
                            stamp TIMESTAMP,
                            zoned_stamp TIMESTAMP WITH TIME ZONE,
                            interval1 INTERVAL YEAR TO MONTH,
                            interval2 INTERVAL DAY TO SECOND
);

INSERT INTO dates_type(dat, stamp, zoned_stamp, interval1, interval2)
VALUES (TO_DATE('01-01-2023', 'DD-MM-YYYY'), '03-AUG-17 11:20:30.45 AM',
        TO_TIMESTAMP_TZ ('21-FEB-2009 18:00:00 -5:00', 'DD-MON-YYYY HH24:MI:SS TZH:TZM'),
        INTERVAL '10-2' YEAR TO MONTH, INTERVAL '4 5:12:10.222' DAY TO SECOND);

CREATE TABLE JSON_DATA (
    data JSON
);

INSERT INTO JSON_DATA(data)
VALUES ('{"phones":[{"type": "mobile", "phone": "001001"}, {"type": "fix", "phone": "002002"}]}');

INSERT INTO JSON_DATA(data) VALUES ('{"bar": "baz", "balance": 7.77, "active":false}');

INSERT INTO JSON_DATA(data) VALUES ('{"reading": 1.230e-5}');

CREATE TABLE numeric_type (
                              number_value NUMBER(6, 2),
                              small_value NUMBER(9),
                              float_value FLOAT(126),
                              binary_float_value BINARY_FLOAT,
                              binary_double_value BINARY_DOUBLE

);

INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (1.987, 1123, 1.2e-4, 1.17549E-38F, 2.22507485850720E-308);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (1.9, 1, 0.55, binary_float_infinity, 0.2222);
INSERT INTO numeric_type(number_value, small_value, float_value, binary_float_value, binary_double_value)
VALUES (13.0, 1244124, 12.123, -binary_float_infinity, 2.2222);

CREATE TABLE uri_types (
    uri URIType
);

INSERT INTO uri_types(uri) VALUES (SYS.URIFACTORY.getURI('/home/oe/doc1.xml'));
INSERT INTO uri_types(uri) VALUES (SYS.URIFACTORY.getURI('/HR/EMPLOYEES/ROW[EMPLOYEE_ID=205]/SALARY'));
INSERT INTO uri_types(uri) VALUES (SYS.URIFACTORY.getURI('https://datagrok.ai'));

CREATE OR REPLACE TYPE mem_type IS VARRAY(10) of VARCHAR2(15);

CREATE TABLE varrays (members mem_type);

INSERT INTO varrays(members) VALUES (mem_type('Brenda','Richard'));

CREATE TABLE XML_DATA (
    data XMLType
);

INSERT INTO XML_DATA (data) VALUES ('<foo>Hello World!</foo>');
INSERT INTO XML_DATA (data) VALUES ('<book><title>Manual</title><chapter>...</chapter></book>');

CREATE TABLE lobs (blob_type BLOB, clob_type CLOB, nclob_type NCLOB);

INSERT INTO lobs (blob_type, clob_type, nclob_type) VALUES (UTL_RAW.CAST_TO_RAW('Grok'), TO_CLOB('Grok'), TO_NCLOB('Grok'));

CREATE TABLE NULL_SAFETY(ch CHAR(10),
                        varch VARCHAR2(10),
                        nch NCHAR(10),
                        nvarch NVARCHAR2(50), dat DATE, stamp TIMESTAMP,
                        zoned_stamp TIMESTAMP WITH TIME ZONE,
                        interval1 INTERVAL YEAR TO MONTH,
                        interval2 INTERVAL DAY TO SECOND, data JSON, number_value NUMBER(6, 2),
                        small_value NUMBER(9),
                        float_value FLOAT(126),
                        binary_float_value BINARY_FLOAT,
                        binary_double_value BINARY_DOUBLE,  uri URIType, members mem_type, xml_data XMLType,
                        blob_type BLOB, clob_type CLOB, nclob_type NCLOB);

INSERT INTO NULL_SAFETY(ch, varch, nch, nvarch, dat, stamp, zoned_stamp, interval1, interval2, data, number_value,
                        small_value, float_value, binary_float_value, binary_double_value, uri,
                        members, xml_data, blob_type, clob_type, nclob_type)
VALUES (NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL);
