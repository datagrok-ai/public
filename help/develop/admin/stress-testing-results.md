<!-- TITLE: Stress testing results -->
<!-- SUBTITLE: -->

# Stress testing

Stress testing is an important part of the [Quality Assurance](quality-assurance.md#stress-testing)
process. To make sure the platform is stable under heavy load (either many users working simultaneously, or executing
CPU-intensive computations), we perform automated stress testing regularly, as part of our continuous integration
process.

The real workflow is emulated by spawning many virtual machines
(in this document, we test for 1, 4, and 10 concurrent test sessions)
that have Selenium installed on them, and executing selected UI tests on them.

We execute all [tests](#tests) for all [environments](#environments), for different numbers of concurrent test sessions.
Note that "N concurrent test session"
means that tests are being executed on N client machines at the exact same time. In real world, 10 concurrent test
sessions could be equivalent to 100 users currently working with the system
(especially taking into account that Datagrok performs as much computations on the client as possible), or 1000
registered users. Of course it highly depends on usage patterns.

TODO: The results will be automatically published in the Datagrok platform as a table with the following columns: (date,
env, test, min duration, max duration, avg duration, error)

TODO: Add the "interpetation" section for each test results section

## Environments

For stress testing, we use a number of differently setup environments to test different aspects of the Datagrok
platform.

* **"Medium"**: single mid-range server, no auto scaling. Designed to demonstrate the platform's performance when
  running on a single server, and identify platform limitations. Suitable for medium-size companies.
* **"Big"**: powerful servers, auto-scaling capability. Suitable for enterprises. Designed to demonstrate scalability
  under load.
* **"Global"**: tens of thousands of concurrent users representing global community.

Currently, this document only contains results for the "medium" environment. There is an ongoing effort for
stress-testing, documenting, and interpreting results from the other environments.

The **test client** is the same for all environments. It is a Docker image containing the Chrome browser with
the [Selenium IDE CL Runner](https://selenium.dev/selenium-ide/docs/en/introduction/command-line-runner)
that we use for running stress tests.

### Environment: medium

This environment is a single mid-range server, no auto scaling. It is designed to demonstrate the platform's performance
when running on a single server, and identify platform limitations. Suitable for medium-size companies.

The **app server** is the following AWS ECS Cluster:

* c5.xlarge
* 4-core CPU (our app server takes advantage of all cores)
* 8GB RAM

The **database server** is a c5.xlarge AWS RDS Postgres

## Tests

Tests are designed how the platform works under different kind of loads. Some of the tests are parameterized (quite
often, the parameter is a size of the dataset the test works with)

### Test: db-query-n

This test queries the chembl database (ligand_eff table) for the first N rows. It uses three pre-created queries
called "10K_query", "100K_query", and "1M_query".

**Steps**:

1. Login
2. Open "Connect to data" view
3. Search "1M_Test" query in the Postgres provider, which will result in a dataframe containing N rows from the
   chembl/ligand_eff table
4. Run query by "N_Test" double-clicking on it
5. In case of successful query run, Selenium will check for the received dataframe and the number of table rows.
6. Otherwise, it will report an error and stop test-script execution.

## Results

Below are results for the db-query-N test. Each test is executed for 1, 4, and 10 concurrent test sessions.

For one user:

| Rows       | Duration, ms |
|------------|--------------|
| 10k        |          187 |
| 100k       |          461 |
| 1M         |        10596 |

For 4 users:

| Rows       | Duration (min), ms | Duration (max), ms | Duration (avg), ms |
|------------|--------------------|--------------------|--------------------|
| 10k        |                134 |                181 |                161 |
| 100k       |                442 |                636 |                545 |
| 1M         |              14004 |              18894 |              16553 |

For 10 users:

| Rows       | Duration (min), ms | Duration (max), ms | Duration (avg), ms |
|------------|--------------------|--------------------|--------------------|
| 10k        |                137 |                404 |                200 |
| 100k       |                769 |               2743 |               1516 |
| 1M         |              21207 |              53856 |              40302 |
