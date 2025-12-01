#! /bin/bash

rain fmt --write ./templates/*.yml

# Merge templates using yq (preserves CloudFormation tags), example:
yq -n 'load("templates/create_network.template.yml") * load("templates/fargate.template.yml") * load("templates/r53.template.yml") * load("templates/marketplace.template.yml") * load("templates/versions.template.yml")' > network-fargate-r53-marketplace.yml
yq -n 'load("templates/create_network.template.yml") * load("templates/fargate.template.yml") * load("templates/dns.template.yml") * load("templates/marketplace.template.yml") * load("templates/versions.template.yml")' > network-fargate-dns-marketplace.yml
yq -n 'load("templates/choose_vpc.template.yml") * load("templates/fargate.template.yml") * load("templates/dns.template.yml") * load("templates/marketplace.template.yml") * load("templates/versions.template.yml")' > vpc-fargate-dns-marketplace.yml
yq -n 'load("templates/choose_vpc.template.yml") * load("templates/fargate.template.yml") * load("templates/r53.template.yml") * load("templates/marketplace.template.yml") * load("templates/versions.template.yml")' > vpc-fargate-r53-marketplace.yml
yq -n 'load("templates/create_network.template.yml") * load("templates/fargate.template.yml") * load("templates/r53.template.yml") * load("templates/basic.template.yml") * load("templates/versions.template.yml")' > network-fargate-r53-basic.yml
yq -n 'load("templates/choose_vpc.template.yml") * load("templates/fargate.template.yml") * load("templates/r53.template.yml") * load("templates/basic.template.yml") * load("templates/versions.template.yml")' > vpc-fargate-r53-basic.yml
yq -n 'load("templates/create_network.template.yml") * load("templates/fargate.template.yml") * load("templates/dns.template.yml") * load("templates/basic.template.yml") * load("templates/versions.template.yml")' > network-fargate-dns-basic.yml
yq -n 'load("templates/choose_vpc.template.yml") * load("templates/fargate.template.yml") * load("templates/dns.template.yml") * load("templates/basic.template.yml") * load("templates/versions.template.yml")' > vpc-fargate-dns-basic.yml

rain fmt --write *.yml
