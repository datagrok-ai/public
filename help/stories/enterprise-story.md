<!-- TITLE: Enterprise -->
<!-- SUBTITLE: -->

# Problem

Software is eating the world, and companies that best manage to understand their data better 
and derive actionable insights from it will be tomorrow’s winners. Big organizations are 
naturally well-positioned to take advantage of the data they collect for the following reasons:

1. They have a lot of data accumulated already, both in terms of volume and variety
2. This data spans different domains (pharma: discovery / clinical / manufacturing / RWE). 
   Being able to connect the dots between these different areas has a lot of value. Big 
   pharma has an advantage here, since smaller companies typically focus on one area
3. There is a lot of internal talent and in-depth understanding of corresponding domains

However, there are plenty of naturally occurring obstacles in big organizations preventing 
them from unlocking the full potential of their data. The data lives in multiple, often isolated, 
silos. The very structure of the organization often defines barriers for data, initiatives, 
and competencies, which hampers collaboration efforts. Legacy critical systems are not going 
anywhere. Disparate sets of tools that do not work together are being used in different 
departments for different purposes. Ad-hoc, built-for-purpose solutions do not play along 
with current IT policies. Commonly used BI tools are inadequate for dealing with scientific data. 
Once a dataset is prepared and passed around, its history is lost. Quite often, the most common 
technology denominator is an Excel file passed in an email.

# Solution

In order to address all of these issues, we developed Grok, a unique, next-generation data 
platform. Out-of-the-box, it offers the most popular common data-related functionality (storage, ingestion, 
transformation, analysis, visualization, publishing, reporting, predicting, etc). It does it in 
a centralized, coherent manner, with all parts based on the same solid infrastructure, 
similar design and usage patterns, and complete audit trail. 

Grok is a web application with auto updates and zero installation costs. The user interface 
can be tailored for particular user roles. 
Cross-cutting concepts like performance, security, and extensibility were thought-through
and carefully integrated into the platform since the very beginning. 
All functionality can be extended or customized, and both server and client 
APIs are exposed for integration purposes. In a way, Grok is an operating system for data 
and algorithms, with customizations, extensions and applications built on top of it. It is 
designed in such a way that no matter who creates new data, scripts, visualization, or other 
artifacts, the whole organization benefits from it.

While the platform offers a solution for most of the typical data-related needs, it doesn’t force you 
to move all your data to a new system. Out-of-the-box, it immediately creates value by providing unmatched data 
exploratory analysis capabilities. From there, you can grow the ecosystem by adding new data 
connections, pipelines, transformations, publishing dashboards, training predictive models, 
connecting with scripts developed in R, Python, and Julia, building applications, and a lot more. 

## Data Retrieval and Governance

If your organization already has ten different databases, the last thing you want is another, 
all-encompassing database #11. Instead, Grok serves as an intelligence layer, and seamlessly 
and efficiently works with all the popular databases (30+ connectors). Connections and queries are 
centrally managed and are subject to our role-based privilege management system. Queries can be 
used as functions in scripts. All objects can be annotated with user-defined metadata. Built-in 
data provenance is used for tools like data lineage, or impact analysis. Database schemas can 
be visually navigated. We are big proponents of the FAIR data principles. Pipelines, data jobs, 
and alerts can be set up. User-defined functions can be used anywhere. We support multiple 
popular file formats, including scientific formats like sdf, edf, and mat; files can simply be 
dragged-and-dropped into the system.

## Data Visualization and Exploration

Our unique technology lets you explore datasets faster and more efficiently than ever, allowing
you to find patterns that were previously impossible to spot, resulting in the acceleration of 
data-driven decisions. For performance reasons, we heavily utilize our proprietary browser-based 
in-memory database technology that allows us to handle hundreds of millions of rows on the
client side, right in the browser.
 
25+ high-performance interactive viewers are available out-of-the box, and there are a few different 
mechanisms for defining your own custom viewers, ranging from editing custom SVGs to building 
JavaScript-based viewers with complete control over its behavior. Multiple viewers have built-in 
functionality not typically found in plotting software, such as regression lines and confidence 
intervals for scatter plots, p-values for box plots, or dendrograms for tree maps. All viewers 
have access to the R or Python, and it’s easy to integrate a R-, Python-, or Julia-based viewers 
into Grok (it will be less interactive, though). Powerful data transformation and editing 
capabilities are built-in. It is also possible to record macros by simply using the UI, and 
later apply them to new datasets. We have tools for quick inspection of data and drill-down 
capabilities. An AI-based system mines the platform’s historical usage to make intelligent 
suggestions regarding different ways to plot and explore the data; the more you use the system, 
the better it becomes.

Data augmentation is another powerful feature of the platform. It provides dynamically-generated
information for the data points you are looking at. This information can also be targeted. 
For instance, when looking at a molecule medicinal chemists will see predictions on the 
compound solubility (also generated on the fly), while people from the legal department will 
see patents related to that molecule. 
 
## Machine Learning and Data Science

It is estimated that data scientists spend about 80% of their time finding, cleaning, and 
organizing data - and this is exactly what our platform excels at. It is also easy to use scripts 
written in statistical languages, and turn them into functions that are available in any other
context (for instance, an ETL or a custom application).

A number of popular techniques are available out-of-the box as interactive tools: multivariate 
analysis, missing value imputation, clustering, and PCA/PLS.

Moreover, predictive modeling is a first-class citizen in the Grok ecosystem, which makes it easy 
even for non-technical users to train, assess, apply, or share predictive models. We support 
popular model building backends, including H2O and OpenCPU, as well as custom models developed in any 
supported languages.

Model deployment is another hot topic nowadays; in order to extract value from the model, it has 
to be used by the right people in the right contexts. Again, since the platform is used for data 
browsing and exploration, we have multiple options available, such as explicitly running models 
against the datasets (including automatic model suggestions), data augmentation, or using models 
in custom applications.

## Cheminformatics

Datagrok provides first-class support for small molecules, as well as most popular building 
blocks for cheminformatics. It understands several popular notations for representing chemical 
(sub)structures, such as SMILES and SMARTS. Molecules can be rendered in either 2D or 3D with 
different visualization options. They can be sketched as well. Chemical properties, descriptors, 
and fingerprints can be extracted. Substructure and similarity searches are supported for both 
in-memory and db datasets (via chemical cartridges). Predictive models that accept molecules 
as an input can be easily trained, assessed, executed, deployed, reused by other scientists, 
and used in pipelines or in info panels. A lot of chemical functionality is exposed in the 
form of functions, which can be used from anywhere - this includes MCS, chemical clustering, 
conversion of molecule identifiers, advanced rendering capabilities, reactions, scaffolds, etc. 
In addition, a few commonly used techniques, such as R-Group analysis, similarity analysis, or diversity 
analysis, are implemented as interactive apps. A number of core viewers are capable of 
rendering molecules where appropriate.

Cheminformatics functionality is built on top of the Grok core, and rather than being a 
stand-alone, separate application, it enriches the ecosystem. That means that by adding
cheminformatics-related functions, existing viewers get the ability to render molecules, 
new chemical functions are now available for use in data transformations, 
new molecule-related information panels will show up when user clicks on a 
molecule in the grid, etc.

## Application development

The platform provides all necessary building blocks for developing modern data-driven
solutions, and a way for a non-developer users to put it all together. However, in terms
of power and flexibility, nothing compares to developing applications using the tools 
and languages native the for the environment, which in our case is JavaScript.

Datagrok exposes a powerful JavaScript API that lets you control anything within the platform, and
a standard way to develop, debug, and publish your custom code. It is possible to
extend and enrich the platform in many ways, including introducing new functions, viewers, 
widgets, data sources, applications, and other objects. What is unique is that should you choose so,
your custom additions will be immediately accessible by other people and teams, therefore
multiplying the platform's value.

In addition to video lessons, samples gallery, and tutorials available to everyone 
we offer few other ways to help companies with the custom development, including 
on-site training for developers, or developing solutions in-house.  

## Enterprise Customizations

We offer different hosting options, including on-prem, cloud, or hybrid. A flexible role-based 
privilege mechanism lets IT set up corresponding roles for users, and tweak the UI accordingly. 
Backup and storage settings can be handled by IT as well. We support a few commonly used 
authentication systems, including Active Directory and OAuth. There is an option to setup 
flexible rules for audit, which might depend on users, roles, and data. There are a number of 
interactive tutorials to help familiarize users with the system, and an interactive forum 
and knowledge base. To keep your finger on the pulse of the ecosystem, we also include usage 
analysis system.
