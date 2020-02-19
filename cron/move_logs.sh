#!/bin/bash

if [ -n "$(ls -A /krummellab/data1/arrao/logs)" ]
then
  mkdir /krummellab/data1/arrao/dated_logs/$(date "+%Y_%m_%d") && \
  mv /krummellab/data1/arrao/logs/* /krummellab/data1/arrao/dated_logs/$(date "+%Y_%m_%d")/
  ln -sf -T /krummellab/data1/arrao/dated_logs/$(date "+%Y_%m_%d") /krummellab/data1/arrao/dated_logs/MOST_RECENT
fi

