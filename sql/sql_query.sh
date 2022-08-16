#!/bin/bash

date=${date:-${1}}

while read query
do
  mysql -A -h darwinsark2.sftp.wpengine.com -u ${{ secrets.DB_USERNAME }} -p ${{ secrets.DB_PASSWORD }} -P 13306 --ssl-ca=${{ secrets.DB_CERTIFICATE }} wp_darwinsark2 < ./sql/query_${query}.sql > ./raw/${date}/${query}.tsv
done < ./sql/query_list.txt
