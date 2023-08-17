#!/bin/bash

while getopts :hn: opt
do
	case $opt in
	h)
		echo use no arguments for the versions in __init_.py and setup.py to be updated to the highest version between the two.
		echo use \-n \<version.number\> to update both to the perscribed version number.
		exit 0
		;;
	n)
		perscribed_version=$OPTARG		
		echo version numbers will be updated to "$perscribed_version"
		;;
	\?)
		echo "unrecognized flag used: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done

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

		if [[ $first_dec > $second_dec && $4 != "FAIL" ]]
		then
			sed -i "s/$2/$1/" $4
			versionsEqual=1
		elif [[ $first_dec < $second_dec && $3 != "FAIL" ]]
		then
			sed -i "s/$1/$2/" $3
			versionsEqual=1
		elif [[ $first_dec == $second_dec ]]
		then
			(( i++ ))
		else # FAIL
			echo FAILED: Perscribed version number is LESS THAN existing version number >&2
			exit 1
		fi
	done
}

setup_file="setup.py"
init_file="hydrophobicity_explorer/__init__.py"

setup_version=$(read_version_number "$setup_file")
init_version=$(read_version_number "$init_file")

# echo "$setup_version"
# echo "$init_version"
# echo "$perscribed_version"

if [[ $perscribed_version != "" ]]
then
	if [[ $perscribed_version == $setup_version || $perscribed_version == $init_version ]]
	then
		echo FAILED: Perscribed version number is EQUAL TO at least one existing version number >&2
		exit 1
	else
		conform_to_greater_vn "$perscribed_version" "$init_version" "FAIL" "$init_file"
		conform_to_greater_vn "$setup_version" "$perscribed_version" "$setup_file" "FAIL"
	fi
else
	conform_to_greater_vn "$setup_version" "$init_version" "$setup_file" "$init_file"
fi

echo version number in both files is $(read_version_number "$setup_file")
