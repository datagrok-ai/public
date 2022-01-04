# Data access

Understanding the data starts with retrieving that data in some way. This is where Datagrok helps
you. In short, Datagrok provides access to anything that is machine-readable. That includes all sort
of major databases, relational databases, non-relational and NoSQL databases, flat files from a
local computer, network file shares and cloud storages such as Dropbox.

Additionally, Datagrok can easily connect to web services. That could be done in multiple flavors.
If a web services with an [OpenAPI][] support is provided, Datagrok can connect to it automatically.
For anything else, writing some custom JavaScript is required to connect using REST API and other
transports.

Datagrok works natively with lots of existing popular data formats, such as CSV, TXT, JSON, and
others. It can also work with tons of scientific formats. The platform has a system that can
natively extend it to understand the very custom formats. This is very handy especially in the areas
of life sciences, such as bioinformatics.

For relational databases, Datagrok allows to connect to them visually, interrogate their schema,
browse the content of the database, and visually build queries (SQL or SparQL). Datagrok can also
connect to lots of publicly available databases, such as Chembl or Unichem, and this list is
growing.

Finally, everything could be automated. For instance, there is a macro that is being generated
automatically while you operate in the platform, and can be replayed at any time, so that all the
data ingestion and transformation steps that you are performing could be scripted and replayed
back at any time.

There is also a concept of a data job, so that the data transformation steps on the way to a
visual dashboard could be neatly arranged as a diagram and later managed and extended.

Generating data from arbitrary functions or scripts in Python and R, which are in turn functions
as well, is also possible.