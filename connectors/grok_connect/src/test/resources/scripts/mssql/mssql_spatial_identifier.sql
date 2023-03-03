CREATE TABLE SpatialTable1( id int IDENTITY (1,1),
  GeomCol1 geometry,
  GeomCol2 AS GeomCol1.STAsText() );

INSERT INTO SpatialTable1 (GeomCol1) VALUES (geometry::STGeomFromText('LINESTRING (100 100, 20 180, 180 180)', 0));

INSERT INTO SpatialTable1 (GeomCol1) VALUES (geometry::STGeomFromText('POLYGON ((0 0, 150 0, 150 150, 0 150, 0 0))', 0));

CREATE TABLE SpatialTable2
( id int IDENTITY (1,1),
  GeogCol1 geography,
  GeogCol2 AS GeogCol1.STAsText() );

INSERT INTO SpatialTable2 (GeogCol1)
VALUES (geography::STGeomFromText('LINESTRING(-122.360 47.656, -122.343 47.656 )', 4326));

INSERT INTO SpatialTable2 (GeogCol1)
VALUES (geography::STGeomFromText('POLYGON((-122.358 47.653 , -122.348 47.649, -122.348 47.658, -122.358 47.658, -122.358 47.653))', 4326));
