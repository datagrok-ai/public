---
title: "AWS S3 backup configuration"
---

This document contains instructions to configure AWS S3 bucket backup

## Configure with terraform AWS deployment

To configure AWS S3 bucket backup use [instructions](./deploy-amazon-terraform.md)

## Manual S3 bucket backup configuration

To configure manual AWS S3 bucket backup do the next steps:

1. Login to your AWS Console and navigate to the S3 service
2. Choose the S3 bucket that you use for datagrok
3. Go to `Properties --> Bucket Versioning` and press `Edit`
4. Set Bucket Versioning `Enable` and press `Save changes`
5. Go to IAM, press `Roles` and press `Create new role`.
6. Choose `AWS Service`, in `Use cases for other AWS services:` choose `AWS Backup` and press `next`
7. In step two choose the correct permissions:
   - Search and add `AWSBackupServiceRolePolicyForS3Backup` policy
   - Press `Create policy`, choose `JSON`, insert the following code, and press `next`

     ```
              "Statement": [
             {
                 "Action": [
                     "s3:ListBucket",
                     "s3:GetBucketLocation",
                     "s3:GetObjectVersion",
                     "s3:GetObjectVersionAcl",
                     "s3:GetObject",
                     "s3:ListBucketMultipartUploads",
                     "s3:*",
                     "backup:CreateBackupPlan",
                     "backup:CreateBackupSelection",
                     "backup:StartBackupJob",
                     "backup:ListBackupPlans",
                     "backup:ListBackupSelections",
                     "backup:ListBackupVaults",
                     "cloudwatch:GetMetricData"
                 ],
                 "Effect": "Allow",
                 "Resource": [
                     "arn:aws:s3:::<Your S3 bucket name>",
                     "arn:aws:s3:::<Your S3 bucket name>/*",
                     "*"
                 ]
             }
         ],
         "Version": "2012-10-17"
     ```

   - Set tag and policy name and save it
   - Add saved policy to IAM Role
   - Set IAM Role name and create it
8. Go to `AWS Backup --> Backup vaults`, press `create backup vault`
9. Choose Vault name and kms key and press `create backup vault`
10. Go to `AWS Backup --> Backup vaults`, press `create backup plan`
11. Choose `Start with a template`, select a template, and write backup plan name
12. Press `Add backup rule`, set backup rule name, choose your backup vault,
    set backup frequency and retention period, press `Add backup rule`
13. Press `Create backup plan`
14. Set `Resource assignment name`, choose your IAM role that was created earlier
15. Next you can include all resource types in `Define resource selection`,
    but we recommend to use include specific resource types, choose `S3`,
    choose your S3 bucket and press `Assign resources`

Now your S3 bucket backup is enabled.
