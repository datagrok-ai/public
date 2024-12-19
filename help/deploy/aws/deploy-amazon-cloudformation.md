---
title: "AWS CloudFormation"
sidebar_position: 1
format: 'mdx'
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

The deployment consists of a few docker containers, [database](../../develop/under-the-hood/infrastructure.md#database) for storing metadata,
and [persistent file storage](../../develop/under-the-hood/infrastructure.md#storage) for storing files

This document contains instructions to deploy Datagrok using [CloudFormation](https://aws.amazon.com/cloudformation/)
on [AWS ECS cluster](https://aws.amazon.com/ecs/) with [AWS RDS](https://aws.amazon.com/rds/)
and [AWS S3](https://aws.amazon.com/s3/).

We considered a lot of typical security nuances during the CloudFormation template development. As a result, you will
create a Datagrok infrastructure in AWS that applies to all standard security policies.

More information about Datagrok design and components:

* [Architecture](../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../develop/under-the-hood/infrastructure.md)

## Prerequisites

1. Check that you
   have [required permissions](https://github.com/datagrok-ai/public/blob/master/help/deploy/iam.list)
   on AWS account to perform CloudFormation deployment to ECS.

## Deploy Datagrok components

We prepared specific template for every need of our customers, answer simple questions below to use the right one for you.

### Would you like to use an existing VPC in your AWS account?

<Tabs groupId="net" queryString>
<TabItem value="vpc" label="Yes">

    Datagrok stand will be put in an existing VPC you choose upon creation.

    #### Do you use Route53 as DNS provider?

    <Tabs groupId="dns" queryString>
    <TabItem value="r53" label="Yes" default>
    
    ##### Requirements
    
    1. Create [Route53 public hosted zone](https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/CreatingHostedZone.html)
    
    ##### How to deploy

    1. Use the [link](https://console.aws.amazon.com/cloudformation/home#/stacks/quickcreate?templateURL=https%3A%2F%2Fdatagrok-data.s3.us-east-2.amazonaws.com%2Fdeployment%2Fvpc-fargate-r53-basic.yml&stackName=datagrok) 
       to open CloudFormation template and fill all required parameters.

        1. [Specify stack name](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html). 
           To meet AWS naming requirements, name must be shorter than _10 symbols_ and correspond [S3 Bucket naming rules](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucketnamingrules.html).
           We use 'datagrok' by default, but you may prefer to also specify env in the stack name.

    2. Wait until AWS completes the deployment. The stack status will be 'CREATE_COMPLETE.' 
       The script created datagrok stand in existing VPC using existing Route53 hosted zone. 
       Your Datagrok instance is now ready to use.

       If you see one of the following statuses then something went wrong: CREATE_FAILED, ROLLBACK_IN_PROGRESS, ROLLBACK_COMPLETE, ROLLBACK_FAILED. [Check the stack events](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/troubleshooting.html#basic-ts-guide) for more information about error.

    3. Enter the platform `datagrok.<subdomain>` using the `admin` user. To get the password:
        1. Go to stack Outputs. Find _DatagrokAdminPassword_ and click on the link to open AWS Secret Manager.
        2. Click '_Retrieve secret value_' and copy _password_. It is a generated password for the first admin login.
        3. To increase security, [change the password for the admin user](../complete-setup/configure-auth.md) on the first login. Datagrok will ignore the admin password from secrets on subsequent restarts.

    4. Complete the initial setup in platform and you are ready to use Datagrok.


    </TabItem>
    <TabItem value="dns" label="No">
    
    Our CloudFormation scripts support external DNS providers, however, it will require a few manual steps to configure the endpoint. 

    ##### Requirements

    1. Come up with two endpoints: `DATAGROK_DNS`, `CVM_DNS`. Datagrok requires two endpoints: `DATAGROK_DNS` and `CVM_DNS`.
       Users will use `DATAGROK_DNS` to access Datagrok Web UI, and requests `CVM_DNS` will be sent automatically by
       Datagrok Client.

    2. Create RSA SSL certificate for `DATAGROK_DNS` and `CVM_DNS`.
    
        * If you use AWS ACM service for SSL certificates
            1. [Generate ACM certificate in AWS](https://docs.aws.amazon.com/acm/latest/userguide/gs-acm-request-public.html)
               which will be valid for both endpoints: `DATAGROK_DNS`, `CVM_DNS`.
            2. Copy AWS ARN for the created certificate. It should look like
               this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.
        * If you do not use AWS ACM service for SSL certificates, you can create a certificate for `DATAGROK_DNS`
          , `CVM_DNS` endpoints any way you are already using. Wildcard certificate also suffices.
            1. [Upload certificate to AWS ACM](https://docs.aws.amazon.com/acm/latest/userguide/import-certificate-api-cli.html)
            2. Copy AWS ARN for the created certificate(s). It should look like
               this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.
    
    ##### How to deploy

    1. Use the [link](https://console.aws.amazon.com/cloudformation/home#/stacks/quickcreate?templateURL=https%3A%2F%2Fdatagrok-data.s3.us-east-2.amazonaws.com%2Fdeployment%2Fvpc-fargate-dns-basic.yml&stackName=datagrok) 
       to open CloudFormation template and fill all required parameters.

        1. [Specify stack name](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html). 
           To meet AWS naming requirements, name must be shorter than _10 symbols_ and correspond [S3 Bucket naming rules](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucketnamingrules.html).
           We use 'datagrok' by default, but you may prefer to also specify env in the stack name.
        2. `DatagrokArnSSLCertificate`: Specify AWS ACM ARN for `DATAGROK_DNS` and `CVM_DNS` from the 2nd prerequisites step.

    2. Wait until AWS completes the deployment. The stack status will be 'CREATE_COMPLETE.' 
       The script created datagrok stand with all required infrastructure from scratch using external DNS service and existing AWS ACM certificate. 
       Your Datagrok instance is now ready to use.

       If you see one of the following statuses then something went wrong: CREATE_FAILED, ROLLBACK_IN_PROGRESS, ROLLBACK_COMPLETE, ROLLBACK_FAILED. [Check the stack events](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/troubleshooting.html#basic-ts-guide) for more information about error.

    3. As you chose the fulfillment option with external DNS, you need to create CNAME DNS records for CVM and Datagrok Load Balancers.
       To get the Load Balancer endpoints for DNS record:

        1. Go to stack Outputs. Copy values for _DatagrokLoadBalancerDNSName_ and _CvmLoadBalancerDNSName_.
        2. Use copied DNS names to create CNAME DNS records, for example
            * Host: `DATAGROK_DNS`, Target: _DatagrokLoadBalancerDNSName_
            * Host: `CVM_DNS`, Target: _CvmLoadBalancerDNSName_
    
    4. Enter the platform `DATAGROK_DNS` using `admin` user. To get the password:
        1. Go to stack Outputs. Find _DatagrokAdminPassword_ and click on the link to open AWS Secret Manager.
        2. Click '_Retrieve secret value_' and copy _password_. It is a generated password for the first admin login.
        3. To increase security, [change the password for the admin user](../complete-setup/configure-auth.md) on the first login. Datagrok will ignore the admin password from secrets on subsequent restarts.
    
    5. Complete the initial setup in platform and you are ready to use Datagrok. 

    </TabItem>
    </Tabs>

</TabItem>
<TabItem value="network" label="No" default>

    Datagrok stand will create VPC and all required network resources itself.

    #### Do you use Route53 as DNS provider?

    <Tabs groupId="dns" queryString>
    <TabItem value="r53" label="Yes" default>
    
    ##### Requirements
    
    1. Create [Route53 public hosted zone](https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/CreatingHostedZone.html)
    
    ##### How to deploy

    1. Use the [link](https://console.aws.amazon.com/cloudformation/home#/stacks/quickcreate?templateURL=https%3A%2F%2Fdatagrok-data.s3.us-east-2.amazonaws.com%2Fdeployment%2Fnetwork-fargate-r53-basic.yml&stackName=datagrok)
    to open CloudFormation template and fill all required parameters.
    
        1. [Specify stack name](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html). 
           To meet AWS naming requirements, name must be shorter than _10 symbols_ and correspond [S3 Bucket naming rules](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucketnamingrules.html).
           We use 'datagrok' by default, but you may prefer to also specify env in the stack name.

    2. Wait until AWS completes the deployment. The stack status will be 'CREATE_COMPLETE.' 
       The script created datagrok stand with all required infrastructure from scratch using existing Route53 hosted zone. 
       Your Datagrok instance is now ready to use.

       If you see one of the following statuses then something went wrong: CREATE_FAILED, ROLLBACK_IN_PROGRESS, ROLLBACK_COMPLETE, ROLLBACK_FAILED. [Check the stack events](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/troubleshooting.html#basic-ts-guide) for more information about error.


    3. Enter the platform `datagrok.<subdomain>` using the `admin` user. To get the password:
        1. Go to stack Outputs. Find _DatagrokAdminPassword_ and click on the link to open AWS Secret Manager.
        2. Click '_Retrieve secret value_' and copy _password_. It is a generated password for the first admin login.
        3. To increase security, [change the password for the admin user](../complete-setup/configure-auth.md) on the first login. Datagrok will ignore the admin password from secrets on subsequent restarts.

    5. Complete the initial setup in platform and you are ready to use Datagrok.

    </TabItem>
    <TabItem value="dns" label="No">
    
    Our CloudFormation scripts support external DNS providers, however, it will require a few manual steps to configure the endpoint.

    ##### Requirements

    1. Come up with two endpoints: `DATAGROK_DNS`, `CVM_DNS`. Datagrok requires two endpoints: `DATAGROK_DNS` and `CVM_DNS`.
       Users will use `DATAGROK_DNS` to access Datagrok Web UI, and requests `CVM_DNS` will be sent automatically by
       Datagrok Client.

    2. Create RSA SSL certificate for `DATAGROK_DNS` and `CVM_DNS`.
    
        * If you use AWS ACM service for SSL certificates
            1. [Generate ACM certificate in AWS](https://docs.aws.amazon.com/acm/latest/userguide/gs-acm-request-public.html)
               which will be valid for both endpoints: `DATAGROK_DNS`, `CVM_DNS`.
            2. Copy AWS ARN for the created certificate. It should look like
               this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.
        * If you do not use AWS ACM service for SSL certificates, you can create a certificate for `DATAGROK_DNS`
          , `CVM_DNS` endpoints any way you are already using. Wildcard certificate also suffices.
            1. [Upload certificate to AWS ACM](https://docs.aws.amazon.com/acm/latest/userguide/import-certificate-api-cli.html)
            2. Copy AWS ARN for the created certificate(s). It should look like
               this: `arn:aws:acm:<region>:<account_id>:certificate/<certificate_id>`.
    
    ##### How to deploy

    1. Use the [link](https://console.aws.amazon.com/cloudformation/home#/stacks/quickcreate?templateURL=https%3A%2F%2Fdatagrok-data.s3.us-east-2.amazonaws.com%2Fdeployment%2Fvpc-fargate-dns-basic.yml&stackName=datagrok) 
       to open CloudFormation template and fill all required parameters.

        1. [Specify stack name](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-using-console-create-stack-parameters.html). 
           To meet AWS naming requirements, name must be shorter than _10 symbols_ and correspond [S3 Bucket naming rules](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucketnamingrules.html).
           We use 'datagrok' by default, but you may prefer to also specify env in the stack name.
        2. `DatagrokArnSSLCertificate`: Specify AWS ACM ARN for `DATAGROK_DNS` and `CVM_DNS` from the 2nd prerequisites step.
        3. 

    2. Wait until AWS completes the deployment. The stack status will be 'CREATE_COMPLETE.' 
       The script created datagrok stand with all required infrastructure from scratch using external DNS service and existing AWS ACM certificate. 
       Your Datagrok instance is now ready to use.

       If you see one of the following statuses then something went wrong: CREATE_FAILED, ROLLBACK_IN_PROGRESS, ROLLBACK_COMPLETE, ROLLBACK_FAILED. [Check the stack events](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/troubleshooting.html#basic-ts-guide) for more information about error.

    3. As you chose the fulfillment option with external DNS, you need to create CNAME DNS records for CVM and Datagrok Load Balancers.
       To get the Load Balancer endpoints for DNS record:

        1. Go to stack Outputs. Copy values for _DatagrokLoadBalancerDNSName_ and _CvmLoadBalancerDNSName_.
        2. Use copied DNS names to create CNAME DNS records, for example
            * Host: `DATAGROK_DNS`, Target: _DatagrokLoadBalancerDNSName_
            * Host: `CVM_DNS`, Target: _CvmLoadBalancerDNSName_
    
    4. Enter the platform `DATAGROK_DNS` using the `admin` user. To get the password:
        1. Go to stack Outputs. Find _DatagrokAdminPassword_ and click on the link to open AWS Secret Manager.
        2. Click '_Retrieve secret value_' and copy _password_. It is a generated password for the first admin login.
        3. To increase security, [change the password for the admin user](../complete-setup/configure-auth.md) on the first login. Datagrok will ignore the admin password from secrets on subsequent restarts.

    5. Complete the initial setup in platform and you are ready to use Datagrok. 

    </TabItem>
    </Tabs>

</TabItem>
</Tabs>

<!--

## Optional: Cost reduction stand

AWS stack uses `FARGATE` instances for deployment by default. To reduce
infrastructure costs, you can use EC2 instances. To do so, follow the instructions above with additions below in
the [prerequisites](#ec2-prerequisites) and [parameters](#ec2-parameters).

### EC2 Prerequisites

1. Before deploying the Datagrok Stand in addition to [Prerequisites](#prerequisites), create RSA key pair. It is
   required to get access to the instances that will be created, you need to have SSH key pair: a private key and a
   public key.

     * If you already have an RSA key pair, you can use the existing one.
     * If you have a Linux-based OS or macOS, type in terminal `ssh-keygen` and hit **Enter**.
       You'll be asked to enter a passphrase. Hit **Enter** to skip this step.
       It will create `id_rsa` and `id_rsa.pub` files in the `~/.ssh` directory.
     * If you have Windows open the Settings panel, then click Apps.
       Under the *Apps and Features* heading, click **Optional Features**.
       Scroll down the list to see if **OpenSSH Client** is listed.
       If it's not, click the plus-sign next to **Add a feature**.
       Scroll through the list to find and select **OpenSSH Client**. Finally, click **Install**.
       Press the **Windows key**, type **cmd** under *Best Match*, and right-click **Command Prompt**.
       Click **Run as Administrator**.
       If prompted, click **Yes** in the *Do you want to allow this app to make changes to your device?* pop-up.
       In the command prompt, type `ssh-keygen` and hit **Enter**.
       You'll be asked to enter a passphrase. Hit **Enter** to skip this step.
       By default, the system will save the keys to `C:\Users\your_username/.ssh/id_rsa`.

2. Copy the content of the public part of the key pair: `id_rsa.pub`. It will be placed in the EC2 instances
   using the `Ec2PublicKey` parameter to access machines.

### EC2 Parameters

1. Change the `LaunchType` parameter to `EC2`.

2. Paste your public key content from [EC2 prerequisites](#ec2-prerequisites) to
   the parameter `Ec2PublicKey`.

-->
