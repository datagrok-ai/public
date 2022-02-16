<!-- TITLE: Deployment on a regular machine -->
<!-- SUBTITLE: -->

# Deployment on a regular machine

Datagrok consist of Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

Like a regular machine, any bare-metal server or virtual machine, including virtual machines in cloud providers, for
example, [AWS EC2](https://aws.amazon.com/ec2/), can be used.

As [database](infrastructure.md#database) Datagrok supports any PostgreSQL database out-of-the-box, including cloud
solutions for PostgreSQL database, for example [AWS RDS](https://aws.amazon.com/rds/).

For [persistent file storage](infrastructure.md#storage), Datagrok supports a lot of options, including cloud solutions,
for example [AWS S3](https://aws.amazon.com/s3/) and Local File System storage.

This document contains instructions to deploy Datagrok using [Docker Compose](https://docs.docker.com/compose/)
on [AWS EC2](https://aws.amazon.com/ecs/) virtual machines with [AWS RDS](https://aws.amazon.com/rds/) as database and
Local File System for persistent storage. This instruction does not cover load balancers creation, which is recommended
for production usage: one load balancer for Datagrok components and one for CVM components.

More information about Datagrok design and components:

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)

In case you want to jump-start using Datagrok with minimum manual effort on a local machine,
check [Local Deployment with Docker Compose](../../develop/admin/docker-compose.md).

## Prerequisites

1. We use native [Docker compose](https://docs.docker.com/compose/) commands to run applications on machines. It
   simplifies multi-container application development and deployment.
    1. Download and install the latest version of [Docker Compose](https://docs.docker.com/compose/install/) to your
       local machine
2. Additional components: instance, database, storage, etc., can be created using AWS CLI. To perform AWS CLI commands
   provided in the document
    1. [Install AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
    2. [Configure authorization for AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-quickstart.html)

## Preparations

The below example contains steps to create EC2 instances as a virtual machine with public IP association. In your case,
it can be any virtual machine. Also, Load Balancers for each VM can be used instead of public IP addresses.

1. Generate SSH key to access virtual machines

    ```shell
    ssh-keygen -t rsa -N '' -m PEM -C 'Datagrok SSH Key' -f ~/.ssh/datagrok-deploy.pem
    ```

2. [Import keypair to AWS](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#how-to-generate-your-own-key-and-import-it-to-aws)
   . Skip this stage if you do not use AWS EC2.

    ```shell
    aws ec2 import-key-pair --key-name datagrok-deploy --public-key-material fileb://~/.ssh/datagrok-deploy.pem.pub
    ```

3. [Create VPC](https://docs.aws.amazon.com/vpc/latest/userguide/vpc-getting-started.html#getting-started-create-vpc)
   for Datagrok EC2 Instances. Skip this stage if you do not use AWS EC2.
    1. [Create VPC](https://docs.aws.amazon.com/vpc/latest/userguide/working-with-vpcs.html#Create-VPC)

       ```shell
       aws ec2 create-vpc --cidr-block '10.0.0.0/17' --output text --query Vpc.VpcId
       ```

    2. [Create Subnet in VPC](https://docs.aws.amazon.com/vpc/latest/userguide/working-with-vpcs.html#AddaSubnet)

        ```shell
        aws ec2 create-subnet --vpc-id "<VPC_ID_FROM_1_STEP>" --cidr-block '10.0.0.0/24' --output text --query Subnet.SubnetId
        ```

    3. [Create and attach internet Gateway to Subnet in VPC](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_Internet_Gateway.html#Add_IGW_Attach_Gateway)

       ```shell
       aws ec2 create-internet-gateway --output text --query InternetGateway.InternetGatewayId
       aws ec2 attach-internet-gateway --vpc-id "<VPC_ID_FROM_1_STEP>" --internet-gateway-id "<IGW_ID_FROM_1_COMMAND>"
       ```

    4. [Create route table with a public route for Subnet in VPC](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_Internet_Gateway.html#Add_IGW_Routing)

        ```shell
        aws ec2 create-route-table --vpc-id "<VPC_ID_FROM_1_STEP>" --output text --query RouteTable.RouteTableId
        aws ec2 associate-route-table --subnet-id "<SUBNET_ID_FROM_2_STEP>" --route-table-id "<ROUTE_TABLE_ID_FROM_1_COMMAND>"
        aws ec2 create-route --route-table-id "<ROUTE_TABLE_ID_FROM_1_COMMAND>" --destination-cidr-block 0.0.0.0/0 --gateway-id "<IGW_ID_FROM_3_STEP>"
        ```

4. [Create Security Group](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html#creating-security-groups)
   for EC2 instances. Skip this stage if you do not use AWS EC2.

    1. [Create Security Group](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html#creating-security-groups)

       ```shell
       aws ec2 create-security-group --group-name datagrok-sg --description "Datagrok SG" --vpc-id <VPC_ID_FROM_4_STAGE>
       ```

    2. [Add a rule for inbound SSH traffic](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/authorizing-access-to-an-instance.html)
       . Note that this rule allows access worldwide

        ```shell
        aws ec2 authorize-security-group-ingress --group-id <SG_ID_FROM_1_STEP> --protocol tcp --port 22 --cidr 0.0.0.0/0
        ```

    3. Add a rule for inbound traffic for Datagrok (8080) and CVM (8090, 5005, 54321). Note that these rules allow
       access worldwide.

       ```shell
       aws ec2 authorize-security-group-ingress --group-id <SG_ID_FROM_1_STEP> --protocol tcp --port 8080  --cidr 0.0.0.0/0
       aws ec2 authorize-security-group-ingress --group-id <SG_ID_FROM_1_STEP> --protocol tcp --port 8090  --cidr 0.0.0.0/0
       aws ec2 authorize-security-group-ingress --group-id <SG_ID_FROM_1_STEP> --protocol tcp --port 5005  --cidr 0.0.0.0/0
       aws ec2 authorize-security-group-ingress --group-id <SG_ID_FROM_1_STEP> --protocol tcp --port 54321 --cidr 0.0.0.0/0
       ```

5. Create a virtual machine for Datagrok components. Requirements: 2 vCPU and 4 GB RAM.
    1. [Create EC2 instance](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html#ec2-launch-instance)
       for Datagrok VM.
        1. Choose AMI with any system which is preferable for you
        2. Press Next
        3. Choose an Instance Type `t3.medium`
        4. Press Next
        5. Choose created in 4th stage VPC for Network
        6. Choose created in 4th stage VPC for Subnet
        7. Auto-assign Public IP: Enable
        8. Press Next
        9. Set Size for Storage to 20 GiB
        10. Press Next
        11. We are okay with default Tags. Press Next
        12. Select an existing security group. And check the security group created in the 5th stage: datagrok-sg
        13. Review and Launch
        14. Launch
        15. Choose an existing key pair which we imported on the 3rd stage: datagrok-deploy

       ```shell
       aws ec2 run-instances --image-id ami-092cce4a19b438926 --block-device-mappings 'Ebs={VolumeSize=20}' --network-interfaces 'AssociatePublicIpAddress=true' --count 1 --instance-type t3.medium --key-name datagrok-deploy --security-group-ids <SG_ID_FROM_5_STAGE> --subnet-id <SUBNET_ID_FROM_4_STAGE>
       ```

6. Create a virtual machine for CVM components. Requirements: 4 vCPU and 8 GB RAM.
    1. [Create EC2 instance](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html#ec2-launch-instance)
       for Compute VM.
        1. Choose AMI with any Linux OS which is preferable for you
        2. Press Next
        3. Choose an Instance Type `c5.xlarge`
        4. Press Next
        5. Choose created in 4th stage VPC for Network
        6. Choose created in 4th stage VPC for Subnet
        7. Auto-assign Public IP: Enable
        8. Press Next
        9. Set Size for Storage to 100 GiB
        10. Press Next
        11. We are okay with default Tags. Press Next
        12. Select an existing security group. And check the security group created in the 5th stage: datagrok-sg
        13. Review and Launch
        14. Launch
        15. Choose an existing key pair which we imported on the 3rd stage: datagrok-deploy

       ```shell
       aws ec2 run-instances --image-id ami-092cce4a19b438926 --block-device-mappings 'Ebs={VolumeSize=100}' --network-interfaces 'AssociatePublicIpAddress=true' --count 1 --instance-type c5.xlarge --key-name datagrok-deploy --security-group-ids <SG_ID_FROM_5_STAGE> --subnet-id <SUBNET_ID_FROM_4_STAGE>
       ```

7. Configure virtual machines
    1. Log in to machines
        1. Use private key created on the first stage for EC2 instances
    2. [Install Docker](https://docs.docker.com/get-docker/) on virtual machines
    3. Add login user to docker group on virtual machines

       ```shell
       sudo usermod -a -G docker <login_user>
       ```

8. Create PostgreSQL 12 database for Datagrok
    1. [Create Security Group](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html#creating-security-groups)
       for EC2 instances. Skip this step if you do not use AWS RDS.

       ```shell
       aws ec2 create-security-group --group-name datagrok-rds-sg --description "Datagrok RDS SG" --vpc-id <VPC_ID_FROM_4_STAGE>
       ```

    2. [Create RDS instance](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/USER_CreateDBInstance.html) for
       Datagrok. Skip this step if you do not use AWS RDS.
        1. DB instance identifier: datagrok-rds
        2. Dev/Test Template
        3. Master username: postgres
        4. Master password and Confirm password: postgres
        5. DB instance class: Burstable classes: `db.t3.medium`
        6. Allocated storage: 50
        7. Enable storage autoscaling
        8. Maximum storage threshold: 100
        9. Do not create a standby instance
        10. Choose created in 4th stage VPC for Virtual private cloud (VPC)
        11. Create a new DB Subnet Group for the Subnet group
        12. Public access: No
        13. VPC security group: Choose existing: select security group created in the 1st step: datagrok-rds-sg

        ```shell
        aws rds create-db-subnet-group \
        --db-subnet-group-name "datagrok-rds" \
        --db-subnet-group-description "DB subnet group for datagrok-rds" \
        --subnet-ids "['<SUBNET_ID_FROM_4_STAGE>']"
        aws rds create-db-instance \
        --db-instance-identifier "datagrok-rds" \
        --db-name "datagrok" \
        --engine 'postgres' \
        --engine-version '12.9' \
        --auto-minor-version-upgrade \
        --allocated-storage 50 \
        --max-allocated-storage 100 \
        --db-instance-class 'db.t3.medium' \
        --master-username "postgres" \
        --master-user-password "postgres" \
        --port "5432" \
        --no-publicly-accessible \
        --storage-encrypted \
        --deletion-protection \
        --backup-retention-period 3 \
        --output text --query 'DBInstance.[DBInstanceIdentifier, DBInstanceStatus]'
        ```

    3. Copy Database address
        1. Copy RDS endpoint

           ```shell
           aws rds describe-db-instances --db-instance-identifier "datagrok-rds" --output text --query 'DBInstances[].[DBInstanceStatus, Endpoint.Address]'
           ```

9. Locally create docker context for virtual machines:

    ```shell
    docker context create --docker 'host=ssh://<DATAGROK_VM_IP_ADDRESS>:22' datagrok
    docker context create --docker 'host=ssh://<CVM_VM_IP_ADDRESS>:22' cvm
    ```

10. Download Docker Compose YAML
    file: [link](https://github.com/datagrok-ai/public/blob/master/docker/localhost.docker-compose.yaml).

## Setup Datagrok components

1. Switch to the datagrok context `docker context use datagrok`
2. Replace in `GROK_PARAMETERS` value with

    ```json
    {
      "dbServer": "<DATABASE_SERVER>",
      "dbPort": "5432",
      "db": "datagrok",
      "dbLogin": "datagrok",
      "dbPassword": "SoMeVeRyCoMpLeXpAsSwOrD",
      "dbAdminLogin": "postgres",
      "dbAdminPassword": "postgres"
    }
    ```

3. Run Datagrok deploy. Wait for the deployment process to complete.

    ```shell
    COMPOSE_PROFILES=datagrok docker-compose --project-name datagrok up -d
    ```

   > **_NOTE:_**  Datagrok provides demo databases with demo data for the full experience.
   > If you want to try datagrok with demo data run the following command instead.
   >
   > ```shell
   > COMPOSE_PROFILES=datagrok,demo docker-compose --project-name datagrok up -d
   > ```

4. Check if Datagrok started successfully: `http://<DATAGROK_VM_IP_ADDRESS>:8080`, login to Datagrok using the
   username "`admin`" and password "`admin`".
5. Switch back to default docker context:

    ```shell
    docker context use default
    ```

## Setup CVM components

1. Switch to the datagrok context `docker context use cvm`

2. Run Datagrok deploy. Wait for the deployment process to complete.

    ```shell
    COMPOSE_PROFILES=cvm docker-compose --project-name cvm up -d
    ```

3. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.
    * Scripting:
        * CVM Url: `http://<CVM_VM_IP_ADDRESS>:8090`
        * CVM URL Client: `http://<CVM_VM_IP_ADDRESS>:8090`
        * H2o Url: `http://<CVM_VM_IP_ADDRESS>:54321`
        * API Url: `http://<DATAGROK_VM_IP_ADDRESS>:8080/api`
        * Cvm Split: `true`
    * Dev:
        * CVM Url: `http://<CVM_VM_IP_ADDRESS>:8090`
        * Cvm Split: `true`
        * API Url: `http://<DATAGROK_VM_IP_ADDRESS>:8080/api`

4. Switch back to default docker context:

    ```shell
    docker context use default
    ```

## Users access

Both [Compute](#setup-cvm-components) and [Datagrok](#setup-datagrok-components) engines should be accessible by users.
The easiest way is to create DNS endpoints pointing to public IPs or load balancers in front of the
services: `datagrok.example`
and `cvm.example`.
