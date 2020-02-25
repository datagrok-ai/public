<!-- TITLE: Datagrok Installation -->
<!-- SUBTITLE: -->

# Datagrok Installation

This document contains instructions for installation of the [DatagrokVM](architecture.md#datlas), 
the server that enables core functionality of the platform. See [GrokCompute Installation](compute-vm-install.md) for instructions
on deploying advanced computational services.  

## AWS bare EC2 installation

1. Get access to Datagrok docker images
2. [Install Docker](https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04)
3. [Install AWS console](https://linuxhint.com/install_aws_cli_ubuntu/) 
4. Login to Amazon `${aws ecr get-login --no-include-email --region us-east-2 --registry-ids 766822877060})`
5. Download Datagrok image `docker pull 766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.32-99f2d0a` 
6. Prepare environment variable
 ```GROK_PARAMETERS
{
"amazonStorageRegion": "us-east-2",                             # S3 region
"amazonStorageBucket": "datagrok-test",                         # S3 bucket name
"amazonStorageId": "ACCOUNTID",                                 # S3 credential ID, Datagrok will resolve EC2 role if empty
"amazonStorageKey": "SECRETKEY",                                # S3 credential secret key, Datagrok will resolve EC2 role if empty
"dbServer": "datagrok-db-1.abc.us-east-2.rds.amazonaws.com",    # RDS endpoint
"db": "datagrok_docker",                                        # RDS new database name
"dbLogin": "datagrok_docker",                                   # RDS new user name
"dbPassword": "SoMeCoMpLeXpAsSwOrD",                            # RDS new user password
"dbAdminLogin": "postgres",                                     # RDS admin login
"dbAdminPassword": "postgres"                                   # RDS admin password
}                                  
```     
7. Run Datagrok image in deploy mode      
`docker run -it -e GROK_PARAMETERS="..." -e GROK_MODE=deploy -p 80:80 766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.27-6cd4234`
8. Wait deploy process to complete
9. Run Datagrok image in regular mode
`docker run -it -e GROK_PARAMETERS="..." -e GROK_MODE=start -p 80:80 766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.27-6cd4234`


## AWS cluster installation

Setup: 
    - EC2: t2.medium x2 (2 vCPU and 4GiB)
    - 1 instance of Datagrok

1. Get access to Datagrok docker images
2. Setup Postgres 11.4 RDS database, allow inbound connections for VPS network
3. Setup Amazon S3 bucket
4. Create "Task Definition" with parameters:
    - Task Definition Name: datagrok
    - Task Role: ecsTaskExecutionRole
    - Network Mode: Host
    - Compatibilities: EC2
    - Task execution role: ecsTaskExecutionRole
    - Task memory (MiB): 2048
    - Task CPU (unit): 2048
    - Container Definitions:
        * Container name: datagrok
        * Image: datagrok:latest
        * Memory Limits (MiB): Soft Limit: 2048
        * Port mappings: 80
        * CPU units: 2048
        * Environment variables (remove comments and make single line): 
            - GROK_MODE = deploy
            - GROK_PARAMETERS = see below
 ```GROK_PARAMETERS
{
"amazonStorageRegion": "us-east-2",                             # S3 region
"amazonStorageBucket": "datagrok-test",                         # S3 bucket name
"amazonStorageId": "ACCOUNTID",                                 # S3 credential ID, Datagrok will resolve ECS task role if empty
"amazonStorageKey": "SECRETKEY",                                # S3 credential secret key, Datagrok will resolve ECS tas if empty
"dbServer": "datagrok-db-1.abc.us-east-2.rds.amazonaws.com",    # RDS endpoint
"db": "datagrok_docker",                                        # RDS new database name
"dbLogin": "datagrok_docker",                                   # RDS new user name
"dbPassword": "SoMeCoMpLeXpAsSwOrD",                            # RDS new user password
"dbAdminLogin": "postgres",                                     # RDS admin login
"dbAdminPassword": "postgres"# RDS admin password
}                                  
```           
           
5. Create cluster (EC2 Linux + Networking):
    - Name: datagrok
    - EC2 instance type: t2.medium
    - Number of instances: 1
    - Key pair: <your ssh key pair>
    - Container instance IAM role: ecsInstanceRole
6. Create service:
    - Launch type: EC2
    - Task Definition: datagrok
    - Cluster: datagrok
    - Service name: datagrok
    - Service type: REPLICA
    - Number of tasks: 1
    - Deployment type: Rolling update
    
    Datagrok starts to deploy database immediately, you can check status by accessing http://public-dns-node-name/api/info/server and
    http://public-dns-node-name/api/admin/deploy/log
    Usually, deploy takes up to 30-40 minutes
    
7. Create new revision for task definition
    - Container Definitions:
        * Modify environment variable
            - GROK_MODE = start
8. Update service, select new task definition
9. Go to the web browser, login to Datagrok using username "admin" and password "SM9ekKEkZuBDp5eD", open Tools | Settings. 
   In sections "Machine Learning" and "Scripting", change all hostnames to the CVM machine hostname. 
  (but leave the part after the hostname intact)
  
## Run docker image
  
1. Download latest docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images)
2. [Install Docker](https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04)
3. Prepare JSON string GROK_START_PARAMETERS:
 ```
{
\"dbServer\": \"host.docker.internal\",     # Postgres database server (use host.docker.internal to connect to localhost)
\"db\": \"datagrok\",                       # New database name
\"dbLogin\": \"datagrok\",                  # New user name
\"dbPassword\": \"SoMeCoMpLeXpAsSwOrD\",    # New user password
\"dbAdminLogin\": \"postgres\",             # Poatgres admin login
\"dbAdminPassword\": \"postgres\"           # Postgres admin password
}                                  
```  
4. Prepare local directory to store data: DATA_PATH   
5. Prepare local directory to store config files: CFG_PATH   
6. Run Datagrok image in deploy mode. Datagrok will create database automatically. 
`docker run -it -v GROK_DATA_PATH:/home/grok/data/prod -v GROK_CFG_PATH:/home/grok/cfg -e GROK_PARAMETERS="GROK_START_PARAMETERS" -e GROK_MODE=deploy -p 80:80 IMAGE_NAME`
7. Wait deploy process to complete
8. Run Datagrok image in regular mode. You can remove dbAdminLogin and dbAdminPassword parameters.
`docker run -it -v GROK_DATA_PATH:/home/grok/data/prod -v GROK_CFG_PATH:/home/grok/cfg -e GROK_PARAMETERS="GROK_START_PARAMETERS" -e GROK_MODE=start -p 80:80 IMAGE_NAME`


