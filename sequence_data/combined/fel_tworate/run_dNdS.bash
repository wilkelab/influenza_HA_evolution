location=`pwd`
echo $location
sed "s|change_path|$location|g" automate_fel_temp.bf > automate_fel.bf

rm errors.log
rm run.log

HYPHYMP automate_fel.bf > run.log 
rm messages.log
