---
title: "Sticky meta"
keywords:
  - metadata
  - tags
  - tags
---

:::note

This feature is in Beta

:::

You can tag anything in Datagrok with custom metadata ("sticky meta"). This sticky meta can serve different purposes. For example, you can annotate a molecule or add arbitrary information about your users. The information about tagged objects and associated metadata are stored in a [Postgres database](../../develop/under-the-hood/infrastructure.md#1-core-components) and seamlessly integrated across the platform, including searching and filtering capabilities based on metadata parameters.

## Configuring metadata

### Metadata schemas

The primary way to use sticky meta is through schemas. Schemas serve as containers for metadata parameters, grouping them into logical categories. For example, a "Rating" schema may include parameters like "Average score" and "Number of votes". 

Metadata schemas are important for several reasons. They: 

* Link metadata parameters with [metadata objects](#metadata-objects)
* Enforce uniformity across your company's metadata
* Reduce errors and accelerate data entry
* Enable advanced search and filtering.

<details>
<summary>Setting up a schema</summary>

To set up a schema:

1. Go to **Browse** > **Platform** > **Sticky Meta** > **Schemas** and click the
   **CREATE A NEW SCHEMA...** button. A dialog opens.
1. In the dialog:
   * Enter the name for your schema.
   * In the **Associated with:** field, link the schema to the relevant metadata object.
   * Add **Properties** to your schema. Define each property with a name and data type. The available data types are `string`, `int`, `bool`, `double`, and `datetime`.
   * Click **OK**.

To preview, edit, and manage schemas, use the [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel).

</details>

### Metadata objects

Metadata objects can be either Datagrok [entities](../../datagrok/concepts/objects.md) or a custom class of objects you create (an "entity type"). For example, to annotate a molecule, create an entity type "Molecule"; to add information about an experiment, create an entity type "Experiment".

<details>
<summary>Creating an entity type</summary>

To create an entity type:

1. Go to **Browse** > **Platform** > **Sticky Meta** > **Types** and click the
   **CREATE A NEW ENTITY TYPE...** button. A dialog opens.
1. In the dialog:
   * Enter the name of your entity type (e.g., `Molecule`).
   * Create a **Matching expression** that defines metadata objects using [tags](tags.md). For tabular data, you would typically use [semantic types](../../govern/catalog/semantic-types.md). For example, to set up a matching expression for an entity type "Molecule", enter `semtype=molecule`. If needed, include several conditions, separating them with commas (e.g., `type=id,belongs=molecule,private=true`).
   * Click **OK**.

</details>

## Applying metadata

Once metadata is configured, Datagrok shows parameters for the selected object in the [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel) under **Sticky Meta**. Here, you can enter or edit metadata as needed. All metadata entries are saved to the [Postgres database](../../develop/under-the-hood/infrastructure.md#1-core-components) in real time.

### Working with dataframes

Datagrok automatically detects columns that match specified [entity types](#metadata-objects) and marks data points within these columns with a blue circle in the top right corner. A dark blue marker means a data point has metadata; a light blue marker suggests you can add metadata.

To add metadata to a data point, click it and enter the details in the **Context Panel** under **Sticky Meta**. From here, you can also create a "sticky column" with parameter values for your entire dataset. To do so, click the **plus (+)** icon next to the parameter or schema name.

Sticky columns are marked by a blue circle in the header. You can add or edit metadata directly in these columns. To add metadata to multiple objects at once, use the [batch edit](../../transform/batch-edit.md) option.

Sticky columns in a dataframe act like standard columns. For example, you can use them to filter data, visualize information using viewers, and so on. However, sticky columns are unique because they link directly to the metadata database. This means removing a sticky column from a dataset doesn't delete its metadata. What's more, metadata updates anywhere on the platform are shared instantly across contexts. For example, commenting on a molecule in a sticky column instantly updates the molecule's tooltip, and other users can see your comment in their datasets too.

## Security and governance

To view a complete list of metadata objects along with their associated metadata, go to **Browse** > **Sticky Meta** > **Types** and click the relevant entity type. Depending on your access rights, you can edit and manage metadata from here.   



