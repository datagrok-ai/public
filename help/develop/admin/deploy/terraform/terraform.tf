terraform {
  required_version = ">= 1.2.0"

  backend "s3" {
    bucket  = "<SET_YOUR_S3_BUCKET_FOR_TERRAFORM_STATE>"
    key     = "datagrok/terraform.tfstate"
    region  = "<SET_YOUR_REGION>"
    encrypt = true
  }
}

provider "aws" {
  region = "<SET_YOUR_REGION>"
}

variable "name" {
  type = string
}

variable "environment" {
  type = string
}

variable "domain_name" {
  type = string
}

variable "docker_hub_user" {
  type = string
}

variable "docker_hub_token" {
  type      = string
  sensitive = true
}

module "datagrok_core" {
  # We recommend to specify an exact tag as ref argument
  source = "git@github.com:datagrok-ai/tf-module-datagrok-core.git//aws?ref=main"

  name                = var.name
  environment         = var.environment
  domain_name         = var.domain_name
  docker_hub_user     = var.docker_hub_user
  docker_hub_password = var.docker_hub_token
}

module "datagrok_cvm" {
  # We recommend to specify an exact tag as ref argument
  source = "git@github.com:datagrok-ai/tf-module-datagrok-cvm.git//aws?ref=main"

  name                             = "${var.name}-cvm"
  environment                      = var.environment
  domain_name                      = var.domain_name
  vpc_id                           = module.datagrok_core.vpc_id
  public_subnet_ids                = module.datagrok_core.public_subnets
  private_subnet_ids               = module.datagrok_core.private_subnets
  docker_hub_secret_arn            = module.datagrok_core.docker_hub_secret
  create_route53_internal_zone     = false
  create_route53_external_zone     = false
  route53_internal_zone            = module.datagrok_core.route_53_internal_zone
  create_cloudwatch_log_group      = false
  cloudwatch_log_group_arn         = module.datagrok_core.cloudwatch_log_group_arn
  cloudwatch_log_group_name        = module.datagrok_core.cloudwatch_log_group_name
  service_discovery_namespace      = module.datagrok_core.service_discovery_namespace
  monitoring_email_alerts          = false
  monitoring_email_alerts_datagrok = false
  monitoring_sns_topic_arn         = module.datagrok_core.sns_topic
  log_bucket                       = module.datagrok_core.log_bucket
}
