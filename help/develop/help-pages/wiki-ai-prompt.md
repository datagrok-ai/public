---
title: Wiki writing guide
sidebar_label: AI prompt for docs
format: mdx
---

## Core philosophy

Every piece of documentation should:
- Be complete on one page ("every page is page one")
- Layer information effectively without overwhelming
- Serve all user types through progressive disclosure
- Use plain English with specific, active verbs
- Start each paragraph with its key idea

## Document structure

### Key elements

```markdown
---
title:  Feature/Topic Name
sidebar_label: Short name // Optional short version to display on the Sidebar
format: mdx
keywords:
 - keyword 1
 - keyword 2 
---

Clear value proposition and core capabilities in 1-2 sentences.

## Core capabilities

Complete feature explanation with examples.
[Layer from basic to advanced within each section]

<details>
<summary>Common operations</summary>
Quick procedures for regular tasks
</details>

To provide information relevant to a specific user role, use admonitions:

:::note user_role
Some technical details or links to in-depth pages.
:::

Example:

:::note developers
Some technical details
:::

## Examples
Real-world usage scenarios with context

## See more
Related content, videos, tutorials
```

### Information layering

Don't split content across pages. Instead, layer content within each section:

1. Main content (visible first)
   - Core concepts
   - Primary capabilities
   - Key workflows
   - Essential examples

2. Expandable details (unless the content size/complexity requires a separate document, in which case link to it)
   - Step-by-step procedures 
   - Configuration options 
   - Specific scenarios
   - Reference information

3. User-specific notes
   - Developer information
   - Administrator guidance
   - Technical details
   - Advanced configurations

## Writing guidelines

### Paragraph construction

Each paragraph must:
- Start with its key idea in first sentence
- Cover one main concept
- Stand independently
- Progress logically

|Good | Bad|
|-----|----|
|The **Feature** let's you create database queries. Design queries by selecting tables and columns, or write SQL directly for complex operations.   | The **Feature** is a powerful tool, useful for different scenarios. It lets you work with databases in multiple ways. You can create queries manually or via UI.   |

### Language principles

- Use simple, short words
- Use active voice
- Choose specific verbs
- Avoid jargon and abstractions
- Explain concepts once, then link
- Write for humans, not machines

### Progressive disclosure

Layer information naturally:

```markdown
## Feature

The **Feature** lets you create database queries. Design queries by selecting tables and columns, or write SQL directly for complex operations.

### Creating queries

Build queries using visual tools:
- Select data sources from the browser
- Choose columns with preview data
- Define relationships through drag-and-drop
- Set filters with dynamic parameters

[GIF]

<details>
<summary>Creating your first query</summary>
Step-by-step guide...
</details>

Write SQL for advanced scenarios:
[code snippet]

:::note developers
Automate query creation: [code snippet] or [Name](/path/to/file)

:::
```

## Document patterns

### Complete feature/capability documentation example

```markdown
---
title: Exploratory data analysis with Datagrok
sidebar_label: EDA
format: mdx
keywords:
 - EDA
 - data exploration
 - visualization
 - analysis
---

Transform, visualize, and analyze your data through interactive tools. Handle millions of rows efficiently, create dynamic visualizations, and share insights with your team.

[Common applications].

### Core capabilities

Analyze your data through multiple approaches:

Data transformation:
- Clean and reshape data
- Handle missing values
- Create calculated columns
- Apply filters and aggregations

Visual analysis:
- Create dynamic charts
- Build interactive dashboards
- Explore relationships
- Track changes over time

Statistical tools:
- Descriptive statistics
- PCA, PLS, multivariate analysis, ANOVA
- Hierarchical clustering

#### Basic data cleaning

Introductory text.

Basic data cleaning:
1. Identify missing values
2. Choose handling method
3. Apply transformations

[GIF]

#### Creating visualizations

Introductory text.

Creating visualizations:
1. Select chart type
2. Choose variables
3. Customize display

[GIF]

### See more
- Tutorial videos
- Example analyses
- Training materials
```

## Content organization patterns

### Progressive detail flow

Structure content to flow naturally:

1. Start with Value
```markdown
Create interactive dashboards that update automatically with fresh data. Build once and share live insights across your organization.
```

2. Show core usage. Add depth with examples
```markdown
### Dashboard creation

Build dashboards using:
- Static data or dynamic data sources
- 50+ interactive viewers
- Smart filters

<details>
<summary>Quick start</summary>
Basic setup steps...
</details>
```

3. Add relevant details for specific user roles
```markdown
:::note developers
You can do X. See [Wiki page](link).
:::
```

### Interactive tool documentation
```markdown
---
title: Visual Query Builder
...
---

Design database queries through an intuitive interface. Create everything from simple data pulls to complex joins and aggregations, with no SQL required.

## Core Functionality

Build queries visually:
- Browse data sources
- Preview live data
- Define relationships
- Set dynamic filters

All operations support:
- Real-time validation
- Performance optimization
- Parameter creation
- Result previews

<details>
<summary>Building your first query</summary>
1. Choose data source
2. Select columns
3. Add conditions
4. Test and save
</details>

## Query Components

### Data selection

Choose your data sources:
- Browse databases
- Search tables and views
- Preview contents
- See relationships

### Query parameters

Make queries dynamic:
- Date ranges
- Value filters
- Multi-select options
- Cascading choices

<details>
<summary>Parameter setup</summary>
Configuration steps...
</details>

:::note developers
Create parameterized queries via API: [code snippet]
:::
```

## Content patterns

### Effective layering techniques

#### Basic to advanced flow

Layer information within each concept:

```markdown

---
title: Chemical structure search
---

Find chemical compounds using structure-based queries.

## Search capabilities

Find compounds using:

- Dataset filters
  - Exact structure search
  - Substructure search
- Built-in chemical database queries (see [Info panes](link))
- Analysis tools
  - Similarity search
  - R-groups analysis

:::note developers
Access structure search programmatically: [code snippet]
:::

### Filter compounds with sketcher

Introductory text.

To filter:
1. Click substructure filter
2. Draw your structure
3. Set search type
4. Execute search


```

#### Cross-referencing
Connect related concepts naturally in text:

Good:
```markdown
Use the structure editor to draw queries, or import structures from your [chemical file shares](path/to/section). For analyzing results, our [chemical visualization tools](path/to/section) provide interactive 3D views.
```

Bad:
```markdown
See also:
- Chemical file shares
- Visualization tools
- Other features
```

### Examples and Use Cases
Embed examples that show complete workflows:

```markdown
## Real-World Applications

### Drug Discovery
Screen compound libraries for potential drug candidates:
1. Draw target structure
2. Set similarity threshold
3. Filter by properties
4. Export matches

Results automatically display in the [chemical grid viewer](path) with:
- 2D structure display
- Property calculations
- Interactive filtering
- Bulk export options

<details>
<summary>Detailed screening steps</summary>
Step-by-step guide...
</details>
```

## Best practices

### Content organization
1. Natural information flow
   - Start with clear value
   - Build on basic concepts
   - Layer in complexity
   - Keep related info together

2. Progressive disclosure
   - Essential info visible
   - Details in expandable sections. Examples show complete context
   - Technical notes separate

3. User path management
   - Clear entry points
   - Natural progression
   - Multiple user journeys
   - Connected capabilities

### Writing Style

1. Clear and direct
   - Use active voice
   - Choose specific verbs
   - Avoid jargon
   - Explain when needed

2. Paragraph structure
   - Key idea first
   - Support with details
   - One concept per paragraph
   - Stand-alone meaning

3. Technical content
   - Layer complexity
   - Show complete context
   - Include real examples
   - Connect related features
