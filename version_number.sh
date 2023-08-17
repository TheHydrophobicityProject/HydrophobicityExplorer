#!/bin/bash

setup_file="setup.py"
init_file="hydrophobicity_explorer/__init__.py"

setup_version=$(awk -F "\"" '/version/ {print $2}' "$setup_file")
init_version=$(awk -F "\"" '/version/ {print $2}' "$init_file")
# echo "$setup_version"
# echo "$init_version"

if [[ "$setup_version" == "$init_version" ]]
then
	versionsEqual=1
else
	versionsEqual=0
fi

i=1
while [ $versionsEqual -lt 1 ]
do
	setup_dec=$(echo "$setup_version" | cut -d "." -f $i)
	init_dec=$(echo "$init_version" | cut -d "." -f $i)

	# echo comparing "$setup_dec" to "$init_dec"

	if [[ $setup_dec > $init_dec ]]
	then
		sed -i "s/$init_version/$setup_version/" $init_file 
		versionsEqual=1
	elif [[ $setup_dec < $init_dec ]]
	then
		sed -i "s/$setup_version/$init_version/" $setup_file
		versionsEqual=1
	else #they are equal
		(( i++ ))
	fi
done

echo version number in both files is $(awk -F "\"" '/version/ {print $2}' "$setup_file")

