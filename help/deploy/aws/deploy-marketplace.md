---
title: "AWS Marketplace"
sidebar_position: 0
---

To provide our users with the easiest way to deploy Datagrok to AWS, we integrated Datagrok with [AWS Marketplace](https://aws.amazon.com/marketplace). Using AWS Marketplace, you can deploy all [Datagrok components](../../develop/under-the-hood/infrastructure.md) and the required AWS infrastructure from scratch with just a few clicks. If you want to customize the deployment or use your existing infrastructure, we suggest using [CloudFormation](deploy-amazon-cloudformation.md) or [Terraform](deploy-amazon-terraform.md) options instead.

<!-- markdownlint-disable no-bare-urls -->
:::info

Datagrok product on AWS Marketplace uses BYOL (Bring Your Own License) pricing model. To purchase a license directly from Datagrok before using the product, please contact info@datagrok.ai.

:::
<!-- markdownlint-enable no-bare-urls -->

## Deployment

1. Open [Datagrok product](https://aws.amazon.com/marketplace/pp/prodview-uqum2jw2yvp52) on AWS Marketplace > Continue to Subscribe > Accept terms. Wait until the subscription request is processed, and then Continue to Configuration.

   <!-- markdownlint-disable no-bare-urls -->
   :::info

   Datagrok product has yet to be publicly available on AWS Marketplace. To request early access, please contact support@datagrok.ai.

   :::
   <!-- markdownlint-enable no-bare-urls -->

2. Choose one of the fulfillment options from the list, and click Continue to launch.

   <details>
   <summary>Fulfillment options description</summary>

    * **ECS [Fargate](https://aws.amazon.com/fargate/) with manual DNS**. Requires to add SSL certificate to ACM and create DNS records manually. It creates all other infrastructure resources. 
    * **ECS [Fargate](https://aws.amazon.com/fargate/) with [Route53](https://aws.amazon.com/route53/) DNS**. Creates all infrastructure from scratch, including SSL certificates and DNS records

   </details>

3. To start deployment, click on the **CloudFormation** deployment template. Quick create stack page opens.

4. Fulfill all parameters for the stack; they can vary depending on the chosen fulfillment option. Click Create stack. The AWS starts the deployment, and the stack has 'CREATE_IN_PROGRESS' status. 

   <details>
   <summary>Parameters details</summary>

   * **Stack name.** It must be shorter than ten symbols to meet AWS naming requirements.

   </details>

5. Wait until AWS completes the deployment. The stack status will be 'CREATE_COMPLETE.' Your Datagrok instance is now ready to use. If you chose the fulfillment option with manual DNS, remember to create CNAME DNS records for CVM and Datagrok Load Balancers.

   :::note
   If you see one of the following statuses then something went wrong: CREATE_FAILED, ROLLBACK_IN_PROGRESS, ROLLBACK_COMPLETE, ROLLBACK_FAILED. Check the stack events for more information about error.
   :::

<!-- GIF with installation --->
