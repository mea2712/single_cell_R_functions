#!/bin/bash -l

<<'//_COMMENT_//'
#### File name [no extension tar.gz] in wharf
args1="filename" 
D1="date"

### Arguments array
var_=$( echo "$(compgen -v | grep 'args' -)" )
read -d " " -a var_array <<< "$var_"

### Get compressed files and compare to home directories
idx=$(#var_array[*])
for i in $(eval echo "{0..$idx}")
do
	lastFile="proj/sens2018122/nobackup/wharf/mararc/mararc-sens2018122/${var_array[i]}_D${i}.tar.gz"
	tmp_dir="proj/sens2018122/nobackup/wharf/mararc/mararc-sens2018122/tmp"
	diff_file="home/mararc/rackham_transfers.d/diff.d/${var_array[i]}_D${i}"
	echo "Comparing:" > $diff_file
	echo "1) $lastFile" >> $diff_file
	echo "2) home/mararc/${var_array[i]}" >> $diff_file
	echo " " >> $diff_file
	mkdir -p $tmp_dir 
	tar --extract -zp --file=$lastFile --directory $tmp_dir
	diff -qr $tmp_dir/${var_array[i]}_D${i} home/mararc/${var_array[i]} >> $diff_file
	
	### Check for confirmation on what to do
	diff_file="home/mararc/rackham_transfers.d/diff.d/${var_array[i]}_D${i}"
	check_diff=$( wc -l < $diff_file)
	if [[check_diff gt 1 ]]
	then
		echo "$lastFile and directory home/mararc/${var_array[i]} ARE equal"
	else
		echo "$lastFile and directory home/mararc/${var_array[i]} ARE NOT equal"
	fi
	
	### Prompt to make user decide whether to continue or not
	while true
	do
		read -p "CONTINUE? (y/n)" choice
		case $choice in
		[nN]* ) exit;;		
		[yY]* ) ### Move file from wharf to home and backup
				mv $lastFile home/mararc/rackham_transfers.d/;
				mv $tmp_dir/${var_array[i]}_D${i} home/mararc/${var_array[i]};
		    
				### Update absolute paths
				grep -rlz 'crex[2]/proj/sllstore2017016/maria.d' /home/mararc/${var_array[i]} | xargs -0 sed -i.bck 's/crex[2]\/proj\/sllstore2017016\/maria.d/home\/mararc/g';
		    
				### Update symbolic links
				link_name=$(find home/mararc/${var_array[i]}/data.d/ -type l); 
				read -d " " -a link_array <<< $link_name;

				target_name=$( find home/mararc/${var_array[i]}/data.d/ -type l -ls | \
							grep -o 'crex[2]/proj/sllstore2017016/maria.d.*$' | \
							awk '{gsub("crex/proj/sllstore2017016/maria.d","home/mararc")}1');
				IFS' ' -a -s target_array <<< $target_name;

				idx1=${#link_array[@]};
				for i in $(eval echo "{0..$idx1}")
				do
					ln -sfn ${target_array[i]} ${link_array[i]}
				done
				exit;;
		* ) echo "Just please enter Y or N. Easy";;
		esac
	done
done
unset i
//_COMMENT_//

# Ask for target directory
read -p "List target directory to update absolute paths: " dir
args=$dir

### Update absolute paths
grep -rlz -I --exclude-dir={output.d,log.d,routs.d} 'crex[2]/proj/sllstore2017016/maria.d' /castor/project/proj/maria.d/${args}/code.d | \
xargs sed -i 's/crex\/proj\/sllstore2017016\/maria.d/castor\/project\/proj\/maria.d/g;s/crex2\/proj\/sllstore2017016\/maria.d/castor\/project\/proj\/maria.d/g';

grep -rlz -I --exclude-dir={output.d,log.d,routs.d} 'crex[2]/proj/sllstore2017016/maria.d' /castor/project/proj/maria.d/${args}/data.d/*/README | \
xargs sed -i 's/crex\/proj\/sllstore2017016\/maria.d/castor\/project\/proj\/maria.d/g;s/crex2\/proj\/sllstore2017016\/maria.d/castor\/project\/proj\/maria.d/g';		    

grep -rlz -I --exclude-dir={output.d,log.d,routs.d} 'snic2019-3-462' /castor/project/proj/maria.d/${args}/code.d/submitters.d/ | \
xargs sed -i 's/snic2019-3-462/sens2018122/g';		    

#mkdir /castor/project/proj/maria.d/${args}/code.d/bck.d
#echo "### DIRECTORY HOLDING BACK UP SCRIPT FROM RUNNING rackham_transfer.sh ###" > /castor/project/proj/maria.d/${args}/code.d/bck.d/README
#mv /castor/project/proj/maria.d/${args}/code.d/*.bck* /castor/project/proj/maria.d/${args}/code.d/bck.d/

###mkdir -p /castor/project/proj/maria.d/${args}/code.d/submitters.d/*/bck.d
###mv castor/project/proj/maria.d/${args}/code.d/submitters.d/*/*.bck* castor/project/proj/code.d/submitters.d/*/bck.d/
