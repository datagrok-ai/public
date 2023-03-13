CREATE TABLE geom (geometry_type GEOMETRY, point_type POINT);
INSERT INTO geom(g) VALUES (ST_GeomFromText('POLYGON((0 0,10 0,10 10,0 10,0 0),(5 5,7 5,7 7,5 7, 5 5))'),
                            ST_GeomFromText('POINT(1 1)'));

INSERT INTO geom(g) VALUES (ST_GeomFromText('GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 1,2 2,3 3,4 4))'),
                            ST_GeomFromText('POINT(1 0)'));

-- SELECT ST_AsText(geometry_type), ST_AsText(point_type) FROM geom;
