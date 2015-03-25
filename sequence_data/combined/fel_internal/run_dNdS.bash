location=`pwd`
echo $location
sed "s|change_path|$location|g" automate_fel_temp.bf > automate_fel.bf

HYPHYMP automate_fel.bf
