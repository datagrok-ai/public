#! /bin/bash

./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o network-fargate-r53-marketplace.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o network-fargate-dns-marketplace.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o vpc-fargate-dns-marketplace.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/marketplace.template.yml -i templates/versions.template.yml -o vpc-fargate-r53-marketplace.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o network-fargate-r53-basic.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/r53.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o vpc-fargate-r53-basic.yml
./templates/cf_template_merge.py -i templates/create_network.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o network-fargate-dns-basic.yml
./templates/cf_template_merge.py -i templates/choose_vpc.template.yml -i templates/fargate.template.yml -i templates/dns.template.yml -i templates/basic.template.yml -i templates/versions.template.yml -o vpc-fargate-dns-basic.yml
