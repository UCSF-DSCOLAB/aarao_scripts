cd /home/${USER}/ihg_library_status
for ihg_username in krummelm
do
    rm -rf ihg-client.ucsf.edu/${ihg_username}/
    IHGUSERNAME=${ihg_username} source /krummellab/data1/immunox/staging/auth_source.sh
    wget --spider \
         --no-check-certificate \
         -r \
         -l 1 \
         -R "index.html*" \
         --no-parent \
         --user=${ihg_username} \
         --password=${IHGPASSWORD} \
         https://ihg-client.ucsf.edu/${ihg_username}/ &> /dev/null
    
    if [ -f ${ihg_username}_files.list ]
    then
        # Get a list of NEW folders on ihg. If IHG has removed any libraries, they won't show up with this command
        # Diff would have shown those.
        new_libraries=( $(grep -vf ${ihg_username}_files.list <(ls ihg-client.ucsf.edu/${ihg_username})) )
        if [ ${#new_libraries[@]} -gt 0 ]
        then
            echo -e "Subject: New libraries identified in ihg/${ihg_username}\n\n" > /tmp/${USER}/ihg_library_status_${ihg_username}.txt
            for new_lib in ${new_libraries[@]}
            do
                echo "https://ihg-client.ucsf.edu/${ihg_username}/${new_lib}" >> /tmp/${USER}/ihg_library_status_${ihg_username}.txt
            done
            /usr/sbin/sendmail arjunarkal.rao@ucsf.edu < /tmp/${USER}/ihg_library_status_${ihg_username}.txt
        fi
    fi
    # Now update the list of files
    ls ihg-client.ucsf.edu/${ihg_username} > ${ihg_username}_files.list
done
