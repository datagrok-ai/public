#! /bin/bash

./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o network-fargate-r53-marketplace.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o network-fargate-dns-marketplace.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o vpc-fargate-dns-marketplace.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o vpc-fargate-r53-marketplace.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o network-fargate-r53-basic.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o vpc-fargate-r53-basic.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o network-fargate-dns-basic.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o vpc-fargate-dns-basic.yml

# EKS variants (non-marketplace only; the EKS + ECS templates share logical IDs
# for RDS and S3 so customers can update an ECS stack to EKS in-place without
# losing data — see deploy-amazon-eks.mdx).
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/eks.template.yml -i templates/r53.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o network-eks-r53-basic.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/eks.template.yml -i templates/r53.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o vpc-eks-r53-basic.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/eks.template.yml -i templates/dns.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o network-eks-dns-basic.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/eks.template.yml -i templates/dns.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o vpc-eks-dns-basic.yml
