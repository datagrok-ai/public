  dep="compute-utils"
  cd $(awk -F'=' '{print $2}' <<<$dep)
  unpublished_deps_deep="$(jq -r '. | to_entries | map(select(.value | match("\\.\\./.*")))[] | "\(.key)=\(.value)"' <<<$dependencies  | tr '\n' ' ')"
  echo $unpublished_deps_deep
  for dep_deep in $(echo -e ${unpublished_deps_deep}); do
    crnt_deep=$(pwd)
    echo "Install dependencies for $(awk -F'=' '{print $2}' <<<$dep_deep)"
    cd $(awk -F'=' '{print $2}' <<<$dep_deep)
    npm install
    npm run build
    cd $crnt_deep
  done
  echo "Install dependencies for $(awk -F'=' '{print $2}' <<<$dep)"
  npm install
  npm run build
