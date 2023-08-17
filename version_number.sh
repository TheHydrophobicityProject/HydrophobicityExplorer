#!/bin/bash

read_version_number(){
	awk -F "\"" '/version/ {print $2}' "$1"
}

is_equal_vn(){
	if [[ $1 == $2 ]]
	then
		echo 1
	else
		echo 0
	fi
}

conform_to_greater_vn(){
	i=1
	versionsEqual=$(is_equal_vn $1 $2)
	while [ $versionsEqual -lt 1 ]
	do
		first_dec=$(echo "$1" | cut -d "." -f $i)
		second_dec=$(echo "$2" | cut -d "." -f $i)

		# echo comparing "$first_dec" to "$second_dec"

		if [[ $first_dec > $second_dec ]]
		then
			sed -i "s/$2/$1/" $4
			versionsEqual=1
		elif [[ $first_dec < $second_dec ]]
		then
			sed -i "s/$1/$2/" $3
			versionsEqual=1
		else #they are equal
			(( i++ ))
		fi
	done
}

setup_file="setup.py"
init_file="hydrophobicity_explorer/__init__.py"

setup_version=$(read_version_number "$setup_file")
init_version=$(read_version_number "$init_file")

# echo "$setup_version"
# echo "$init_version"


conform_to_greater_vn "$setup_version" "$init_version" "$setup_file" "$init_file"

echo version number in both files is $(read_version_number "$setup_file")

