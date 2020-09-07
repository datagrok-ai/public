<!-- TITLE: Public repository -->
<!-- SUBTITLE: -->

# Public Repository

Datagrok is an extensible data-agnostic data analysis platform. It can be extended 
with different scientific methods (functions) that are automatically discovered and suggested to user
when the relevant dataset is open. For instance, for small molecules we have implemented 
popular cheminformatics methods (such as R-group analysis or substructure search), 
and when a user clicks on a molecule they see all functions applicable to it, as well
as different [visualizations](../discover/data-augmentation.md).

These functions are registered in the platform in the form of [scripts](../compute/scripting.md)
developed primarily in Python or R. Scripts can be either registered in the platform ad-hoc,
or stored in the public repository: [https://github.com/datagrok-ai/public](https://github.com/datagrok-ai/public).

Each function can be either executed explicitly, or suggested by the platform whenever a user
opens a dataset relevant to that function. 

All of the packages are globally available to everyone for use from within the Datagrok platform, 
although for the [public environment](https://public.datagrok.ai)
there are some restrictions related to the server compute capacities.

Enterprises with the platform deployed [on-premises](../develop/admin/architecture.md#deployment) 
can also access the public repository.  

## How to contribute

We are interested in other fields, such as NLP, bionformatics, GIS, and others. We are looking for experts (PhD is a plus) 
who would be able to both identify these popular methods, and integrate them into the platform.
They should be proficient in R or Python, and enjoy working with and contributing to open-source projects.

If you are interested in contributing to the public repository, please reach out to us: [info@datagrok.ai](mailto:info@datagrok.ai)

## Academic license

Combined with the provision of a special license to academic institutions, it stands to incentivize academicians to use
the tool, and in the process of doing so, publish new methods in an auto-discoverable and machine-consumable format.
Moreover, the platform is capable of automatically suggesting applicable methods based on a currently open dataset
even for users not connected to the original research. If done right, this cross-pollination of knowledge could be
transformative within and across a broad range of scientific disciplines.

See also:
* [Info panels](../discover/info-panels.md)
* [Self-learning platform]()
* [Predictive modeling](../learn/predictive-modeling.md)
