---
title: "Release notes policy"
position: 5 # float position is supported
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

It is important to keep track of any significant changes made in each version of software. This not only helps users understand the updates, but also encourages them to use the latest version available. 

As a software developer, it is your responsibility to maintain a clean and clear release notes for your project. To do so, please refer to our [guideline](#changelog-guideline) on how to create a changelog.

## Changelog guideline

1. Every Datagrok product, including [Datagrok package](https://datagrok.ai/help/develop/#packages) and [Datagrok library](https://github.com/datagrok-ai/public/tree/master/libraries), should have a `CHANGELOG.md` file in its root directory.
2. The changelog should be kept up-to-date with every new version release. The latest version should be listed first in the changelog along with the date of release, a brief description of what has been updated, and datagrok api dependency version.
   <details>
   <summary>Example</summary>
   
    ```
    ## 97.89.83 (2023-06-31)
   
    This release focuses on improving data access speed and convenience, 
    new visualization and usability features, and ensuring platform stability.
   
    *Dependency: datagrok-api >= Y.Y.Y*
    ```

   <details>
   <summary>Result</summary>

   <!-- markdownlint-disable heading-start-left -->

   ## 97.89.83 (2023-06-31)

   This release focuses on improving data access speed and convenience,
   new visualization and usability features, and ensuring platform stability.

   *Dependency: datagrok-api >= Y.Y.Y*

   <!-- markdownlint-enable heading-start-left -->

   </details>

   </details>
3. Group the changes by type, such as Features (all improvements also go here), Bug Fixes, Breaking change. 
   <details>
   <summary>Example</summary>

   ```
   ### Features

   * Added [logger for packages](https://datagrok.ai/help/develop/advanced/debugging#logger) to report debug records to the server
   * [#1988](https://github.com/datagrok-ai/public/issues/1988): Improved the ability to resize a legend on [Trellis Plot](https://datagrok.ai/help/visualize/viewers/trellis-plot)
   
   ### Bug Fixes

   * [#1984](https://github.com/datagrok-ai/public/issues/1984): Filter's missing values settings are not properly synced between different tabs/views

   ### Breaking change
   
   * Removed Grok Connect from Datagrok image. The host in connectors host should be changed to grok_connect instead of localhost.
   ```

   <details>
   <summary>Result</summary>

   <!-- markdownlint-disable heading-start-left -->

   ### Features

   * Added [logger for packages](https://datagrok.ai/help/develop/advanced/debugging#logger) to report debug records to the server
   * [#1988](https://github.com/datagrok-ai/public/issues/1988): Improved the ability to resize a legend on [Trellis Plot](https://datagrok.ai/help/visualize/viewers/trellis-plot)

   ### Bug Fixes

   * [#1984](https://github.com/datagrok-ai/public/issues/1984): Filter's missing values settings are not properly synced between different tabs/views

   ### Breaking change

   * Removed Grok Connect from Datagrok image. The host in connectors host should be changed to grok_connect instead of localhost.

   <!-- markdownlint-enable heading-start-left -->

   </details>

   </details>

4. When writing the changelog, keep in mind that its primary purpose is to be read by humans in order to understand the changes. Avoid adding insignificant details and commit listing, focus on noteworthy changes.
   <details>
   <summary>Example</summary>

   **Bad**:
   * #1282: fixed molecule size when drawing on canvas
   * Implement Logger
   * Fix error message, bump version
   
   **Good**:
   * [#1282](https://github.com/datagrok-ai/public/issues/1282) Fixed molecule size when drawing on canvas
   * Added [logger for packages](https://datagrok.ai/help/develop/advanced/debugging#logger) to report debug records to the server

   </details>
5. Remember to provide links to the relevant pages in [wiki](https://datagrok.ai/help), instead of including a lengthy explanation of the feature in the changelog. If a GitHub issue was fixed, make sure to mention it.
   <details>
   <summary>Example</summary>

   **Bad**:
   * Added Elemental Analysis which can analyze the elemental composition of a molecular structure and visualizes the results in a radar viewer. 
     To use Elemental Analysis: 
     In the Menu Ribbon, open the Chem menu and select Analyze structure > Elemental Analysis... A parameter input dialog opens.
     Select the source table and the molecular column that you want to analyze.
     Select the desired visualization option.
     Click OK to execute the analysis.

   **Good**:
   * Added [Elemental Analysis](https://datagrok.ai/help/datagrok/solutions/domains/chem/#elemental-analysis) to analyze the elemental composition of a molecular structure.

   </details>

## Changelog example

```markdown
# `<PACKAGE NAME>` changelog

## X.X.X (YYYY-MM-DD)

This release focuses on improving feature stability and usability.

*Dependency: datagrok-api >= Y.Y.Y*

### Features

* New [super cool feature](https://datagrok.ai/help/super-cool-feature)
* Moved Feature from the Context Pane to the Top Menu ( Feature > Calculate).

### Bug Fixes

* [#123](https://github.com/datagrok-ai/public/issues/123) Fixed very nasty bug in feature

### Breaking change

* Retired `copyCoolFeature` method in favor of `duplicateCoolFeature`

```

<!-- markdownlint-disable heading-start-left -->
<!-- markdownlint-disable no-duplicate-heading -->

<details>
<summary>Result</summary>

# `<PACKAGE NAME>` changelog

## X.X.X (YYYY-MM-DD)

This release focuses on improving feature stability and usability.

*Dependency: datagrok-api >= Y.Y.Y*

### Features

* New [super cool feature](https://datagrok.ai/help/super-cool-feature)
* Moved Feature from the Context Pane to the Top Menu ( Feature > Calculate).

### Bug Fixes

* [#123](https://github.com/datagrok-ai/public/issues/123) Fixed very nasty bug in feature

### Breaking change

* Retired `copyCoolFeature` method in favor of `duplicateCoolFeature`

</details>

<!-- markdownlint-enable heading-start-left -->
<!-- markdownlint-enable no-duplicate-heading -->
