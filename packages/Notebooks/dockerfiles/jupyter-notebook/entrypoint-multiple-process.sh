#!/bin/bash

while [[ $# -gt 1 ]]
do
    key="$1"
    case $key in
        --main)
        main=$2
        shift
        ;;
        --helper)
        helper=$2
        shift
        ;;
        *)
        echo "option \'$1\' is not understood!"
        exit -1
        break
        ;;
    esac
    shift
done

if [ -z "$main" ] || [ -z "$helper" ]
then
  echo 'You need to specify both main and helper processes'
  exit 1
fi

# Start the second process
$helper &
helper_pid=$!
sleep 2
helper_status=$(ps -p $helper_pid -o pid=)
if [ -z "$helper_status" ]
then
  echo "$helper has already exited."
  exit 1
fi

# Start the first process
$main &
main_pid=$!
sleep 2
main_status=$(ps -p $main_pid -o pid=)
if [ -z "$main_status" ]
then
  echo "$main has already exited."
  exit 1
fi

# Naive check runs checks once a minute to see if either of the processes exited.
# This illustrates part of the heavy lifting you need to do if you want to run
# more than one service in a container. The container exits with an error
# if it detects that either of the processes has exited.
# Otherwise it loops forever, waking up every 60 seconds

while sleep 10; do
  MAIN_STATUS=$(ps -p $main_pid -o pid=)
  HELPER_STATUS=$(ps -p $helper_pid -o pid=)
  # If the greps above find anything, they exit with 0 status
  # If they are not both 0, then something is wrong
  if [ -z "$MAIN_STATUS" -o -z "$HELPER_STATUS" ]
  then
    echo "One of the processes has already exited."
    exit 1
  fi
done
