#repairs broken links by changing a targeting directory

for i in $(find -lname '*');
do
echo "Found link..." "$i"
liol=$(readlink $i)
linw=$(echo $liol | sed 's/gems_sata/dhts_gems_sata/')
#bn=$(basename $i)
echo "Location" $liol
echo "New location" $linw

ln -nsf $linw $i
 
done
