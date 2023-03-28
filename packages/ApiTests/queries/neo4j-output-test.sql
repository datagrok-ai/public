//name: Neo4jAllTypes
//connection: Neo4jApiTests
//test: ApiTests:expectTable(Neo4jAllTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/neo4j/neo4j_output_all.d42'))
MATCH(t:PropertyTest) RETURN t.datetime as datetime, t.localdatetime as localdatetime,
    t.localtime as localtime, t.time as time, t.float as float,
    t.integer as integer, t.point as point;
//end
