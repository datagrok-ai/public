# HitTriage

The HitTriage package is a powerful tool designed for molecule analysis and campaign management within the Datagrok environment. It consists of two applications: HitTriage and HitDesign. This README provides an overview of the package's functionalities, and subsequent readmes will dive deeper into each application's usage.

## Features

Common Workflow

1. **Template Creation:**

    Define a template specifying the data source for molecules, name, key, additional needed information and compute functions. This source can be a file upload or a query in any other package tagged with `HitTriageDataSource` tag.
    The Compute functions are collected from any package with a tag `HitTriageFunction`.

    ![template](files/images/template.png)

2. **Campaign Building**:

    Create campaigns based on the template.
    Provide a campaign name, select the data source, provide additional information and initiate the campaign.
    During the campaign run, the specified compute functions are executed, and their results are appended to the dataframe. For example, you can compute molecular descriptors, toxicity risks, structural alerts and more.

    ![template](files/images/campaign.png)

After running a campaign, you can submit the dataframe to any chosen function or query. or
save the campaign for later use. Saved campaigns can be reloaded and run again by any user on the platform usign the link or the campaigns table on the first page.

