#!/bin/bash
# Add this to crontab with 
# */15 * * * * bash /krummellab/data1/arrao/scripts/cron/check_job_completion.sh

jobids=(1796021 1796115 1796116 1796117 1796118 1796119 1796120 1796121 1796122 1796123 1796124 1796125 1796126 1796127 1796128 1796129 1796130 1796131)
for jid in ${jobids[@]}
do
  if ! /bin/grep -q ${jid} <(/usr/bin/qstat -u arrao)
  then
    if [ ! -f /tmp/arrao/job_completed_${jid}.txt ]
    then
      cat << EOF > /tmp/arrao/job_completed_${jid}.txt
Subject: Job ${jid} completed

Job completed: ${jid}
EOF
      /usr/sbin/sendmail arjunarkal.rao@ucsf.edu < /tmp/arrao/job_completed_${jid}.txt
    fi
  fi
done
