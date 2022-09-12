#!/bin/bash

##################################
#  This is an automated script that
# runs every night on ccls15 using 
# CRON. DO NOT MODIFY THIS FILE
#
# --Arjun
##################################

email_file=/tmp/ihg_library_status_$(date "+%Y_%m_%d_%H_%M_%S").txt
echo -e "Subject: [TEST] New libraries identified and downloaded from ihg/${ihg_username}\n\n" >> ${email_file}
echo -e "***** This is a test! ***** \n\n" >> ${email_file}
for ihg_username in krummelm combesa immunox
do
    cd /home/${USER}/ihg_library_status
    rm -rf ihg-client.ucsf.edu/${ihg_username}/
    export WGETRC=/data1/immunox/staging/wgetrcs/${ihg_username}
    wget --spider \
         -r \
         -l 1 \
         -R "index.html*" \
         --no-parent \
         --no-check-certificate \
         https://ihg-client.ucsf.edu/${ihg_username}/ &> /dev/null
        
    if [ -f ${ihg_username}_files.list ]
    then
        auto_lib_file=/data1/immunox/staging/ihg_auto/${ihg_username}_$(date "+%Y_%m_%d_%H_%M_%S").list
        # Get a list of NEW folders on ihg. If IHG has removed any libraries, they won't show up with this command
        # Diff would have shown those.
        new_libraries=( $(grep -vf ${ihg_username}_files.list <(ls ihg-client.ucsf.edu/${ihg_username})) )
        if [ ${#new_libraries[@]} -gt 0 ]
        then
            echo "Libraries for username: ${ihg_username}" >> ${email_file}
            for new_lib in ${new_libraries[@]}
            do
                echo "https://ihg-client.ucsf.edu/${ihg_username}/${new_lib}" >> ${email_file}
                echo "https://ihg-client.ucsf.edu/${ihg_username}/${new_lib}" >> ${auto_lib_file}
            done
            # Now actually download the libraries
            #cd /data1/immunox/staging
            #wget --no-check-certificate \
            #     -m \
            #     -R "index.html*" \
            #     --no-parent \
            #     -c \
            #     ${auto_lib_file}
        fi
    fi
    # Now update the list of files for next
    ls ihg-client.ucsf.edu/${ihg_username} > ${ihg_username}_files.list
done

if [[ $(wc -l ${email_file} | cut -d " " -f 1) -gt 1 ]]
then
    echo "\n\n Remember: Data is always downloaded to /krummellab/data1/immunox/staging/ihg-client.ucsf.edu" >> ${email_file}
    /usr/sbin/sendmail dscolab@ucsf.edu < ${email_file}
fi

