<!-- TITLE: Deployment on AWS EC2 -->
<!-- SUBTITLE: -->

# Deployment on AWS EC2

This document contains instructions to deploy Datagrok on an AWS EC2 instance with RDS as a database and S3 as a
storage.

## Prerequisites

1. Create t2.medium EC2 instance for Datagrok VM and c5.xlarge for Compute VM.
   *Public IP address association is required for the instances in case you do not use Application Load Balancer in
   front of the instances*
2. Create DNS records for Datagrok public IP and CVM public IP
3. Configure S3 bucket and RDS database, which should be available from Datagrok VM

## Setup Datagrok virtual machine

1. Create docker network for Datagrok: `docker network create datagrok`
2. Pull Datagrok images to docker

- `docker pull datagrok/datagrok:latest`
- `docker pull datagrok/grok_connect:latest`

3. Prepare string `GROK_PARAMETERS`:

```json
{
  "amazonStorageRegion": "us-east-2",                             # S3 region
  "amazonStorageBucket": "datagrok-test",                         # S3 bucket name
  "amazonStorageId": "ACCOUNTID",                                 # S3 credential ID, Datagrok will resolve EC2 role if empty
  "amazonStorageKey": "SECRETKEY",                                # S3 credential secret key, Datagrok will resolve EC2 role if empty
  "dbServer": "datagrok-db-1.abc.us-east-2.rds.amazonaws.com",    # RDS endpoint
  "db": "datagrok_docker",                                        # RDS new database name
  "dbLogin": "datagrok_docker",                                   # RDS new user name, Datagrok will use it to connect to Postgres database
  "dbPassword": "SoMeVeryCoMpLeXpAsSwOrD",                        # RDS new user password, Datagrok will use it to connect to Postgres database
  "dbAdminLogin": "postgres",                                     # RDS admin login
  "dbAdminPassword": "postgres"                                   # RDS admin password
}
```

4. Run Grok Connect image

```bash
docker run -it -d \
  --network datagrok \
  --network-alias grok_connect \
  --restart unless-stopped \
  datagrok/grok_connect:latest
```

5. Run Datagrok image

```bash
docker run -it -d \
  -e GROK_PARAMETERS="<GROK_PARAMETERS>" \
  --network datagrok \
  --network-alias datagrok \
  -p 80:8080 \
  --restart unless-stopped \
  datagrok/datagrok:latest
```

6. Check if Datagrok started successfully: `http://<DATAGROK_DNS>`, login to Datagrok using username "
   admin"
   and password "admin"

7. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.
    * Connectors
        * External Host: `grok_connect`

## Setup Compute Virtual Machine

1. Create docker network for CVM: `docker network create cvm`
2. Pull CVM images

- `docker pull datagrok/jupyter_kernel_gateway:latest`
- `docker pull datagrok/jupyter_notebook:latest`
- `docker pull datagrok/grok_compute:latest`
- `docker pull datagrok/h2o:latest`
- `docker pull datagrok/cvm_nginx:latest`

3. Run CVM images

```bash
# Grok Compute
docker run -it -d \
  -e GROK_COMPUTE_NUM_CORES=4 \
  --network cvm \
  --network-alias grok_compute \
  --restart unless-stopped \
  datagrok/grok_compute:latest
# H2O
docker run -it -d \
  --network cvm \
  --network-alias h2o \
  -p 54321:54321 
  -p 5005:5005 
  --restart unless-stopped \
  datagrok/h2o:latest
# JKG
docker run -it -d \
  --network cvm \
  --network-alias jupyter_kernel_gateway \
  --restart unless-stopped \
  datagrok/jupyter_kernel_gateway:latest
# JN
docker run -it -d \
  --network cvm \
  --network-alias jupyter_notebook \
  --restart unless-stopped \
  datagrok/jupyter_notebook:latest
# CVM Nginx
docker run -it -d \
  --network cvm \
  --network-alias cvm \
  -p 80:8090 
  --restart unless-stopped \
  datagrok/cvm_nginx:latest
```

4. Check if CVM started successfully: `http://<CVM_DNS>/jupyter/helper/info`

5. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.
   * Scripting:
      * CVM Url: `http://<CVM_DNS>`
      * CVM Url Client: `http://<CVM_DNS>`
      * H2o Url: `http://<CVM_DNS>:54321`
      * Api Url: `http://<DATAGROK_DNS>/api`
      * Cvm Split: `true`
   * Dev:
      * CVM Url: `http://<CVM_DNS>`
      * Cvm Split: `true`
      * Api Url: `http://<DATAGROK_DNS>/api`

See also:

* [Architecture](architecture.md#application)
* [Architecture Details](architecture-details.md)
* [Compute VM](compute-vm.md)
