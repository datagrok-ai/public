//name: Categories
//friendlyName: Product categories provided by each supplier
//connection: Neo4jNorthwind

MATCH (s:suppliers)-->(:products)-->(c:categories)
RETURN s.companyname AS company, collect(DISTINCT c.categoryname) AS category;

//end

//name: Suppliers
//friendlyName: Suppliers by category
//connection: Neo4jNorthwind
//input: string category = "Produce" {pattern: string}

MATCH (c:categories)<--(:products)<--(s:suppliers)
  WHERE c.categoryname = @category

RETURN DISTINCT s.companyname AS suppliername;

//end

//name: Purchases
//friendlyName: Total purchased products by customer
//connection: Neo4jNorthwind
//input: string category = "Produce" {pattern: string}

MATCH (cust:customers)-[:purchased]->(:orders)-[i:includes]->(p:products),
      (p)-[:partof]->(c:categories)
  WHERE c.categoryname = @category

RETURN DISTINCT cust.contactname AS customername, sum(i.quantity) AS totalproductspurchased;

//end
