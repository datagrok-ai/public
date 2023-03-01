---
title: "Substructure search"
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

Substructure search is a task of finding the query substructure pattern among other structures.  
You can set the query substructure either with Sketcher or just selecting a molecule in your dataset.

## Data sources

Substructure search supports several data sources: uploaded datasets, public databases, and relational databases.

### Uploaded datasets

In uploaded datasets, you can run substructure search in three ways:

```mdx-code-block
<Tabs>
<TabItem value="First way" default>
```

1. Right click on a molecule
2. In context menu, click **Current Value** > **Use as filter**

![Substructure search](substructure-search-uploaded-dataset-3.gif "Substructure search")

```mdx-code-block
</TabItem>
<TabItem value="Second way" default>
```

1. On the left sidebar, hover over **Tables** section > <i class='fa fa-filter'></i> **Filter icon**.
2. Find molecules' column > **Click to edit**.
3. Sketch a query substructure. Appropriate rows get filtered with highlighted substructures.

![Substructure search](substructure-search-uploaded-dataset-1.gif "Substructure search")

```mdx-code-block
</TabItem>
<TabItem value="Third way" default>
```

1. On the top menu, click **Chem** > **Sketcher**.
2. Sketch a query substructure
3. Click **Filter** or **Select**

![Substructure search](substructure-search-uploaded-dataset-2.gif "Substructure search")  

```mdx-code-block
</TabItem>
</Tabs>
```

### Public databases

  We support substructure search for 3 public databases: PubChem, Drugbank, and ChEMBL. ChEMBL, in turn, can be accessed by either API or our internal backup.

  Right in info panels you can run substructure search and see the results.

<!-- GIF, but not all the public DBs work -->

### Relational databases

  Chemical cartridges on our platform let you run substructure search in relational databases.

  ![DB Substructure and Similarity Search](../../uploads/gifs/db-substructure-similarity-search.gif "DB Substructure and Similarity Search")

## Advanced options

Datagrok extended substructure search functionality with multiple advanced options:

* Sketcher-oriented structural adjustment of the found molecules  

  When you sketch a molecule in Sketcher, molecules in grid adjust to how molecule is oriented in Sketcher.

* Explicit hydrogen search support

* Extended toolset for aromaticity search
* Multicolumn filtering
