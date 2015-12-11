###########################################################################
# CompleteTCR - repertoire determination from RNAseq of paired TCR chains #
# Copyright (C) 2015  Kami E. Chiotti chiotti@ohsu.edu                    #
#                                                                         #
# This program is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
###########################################################################

#!/usr/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V

#Usage: completeTCR.sh <manifest> <outdir> <Path/To/CompleteTCR> <opt_TF_tosavetmps>
#Manifest format: Full/Path/To/file_of_alpha_fastq:Full/Path/To/file_of_beta_fastq:Prefix_for_paired_output
###

mani=$1
outdir=$2
tooldir=$3
savetmp=$4

if [[ -z "$mani" || -z "$outdir"  || -z "$tooldir" ]]
then
	echo "ERROR: Missing arguments. Run aborted"
	exit 1
fi

if [ ! -f $mani ]
then
	echo "ERROR: Cannot locate manifest file $manii. Run aborted."
	exit 1
elif [ ! -d $outdir ]
then
	echo "ERROR: $outdir does not exist. Run aborted"
	exit 1
elif [ ! -d $tooldir ] 
then
	echo "ERROR: Cannot locate CompleteClone. Run aborted"
	exit 1
else
	echo "~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~"
fi

for line in $(cat $mani)
do
		inputA=`echo $line | cut -d: -f1`
		inputB=`echo $line | cut -d: -f2`
		name=`echo $line | cut -d: -f3`

		if [[ -z "$inputA" || -z "$inputB" || -z "$name" ]]
		then
			echo "ERROR: Improper manifest format. Run aborted."
			exit 1
		fi

		if [ ! -f $inputA ]
		then
			echo "$name error: $inputA does not exist."
		elif [ ! -f $inputB ]
		then
			echo "$name error: $inputB does not exist."
		elif [ ! -d "$outdir/$name" ]
		then	
			mkdir $outdir/$name
			cd $outdir/$name
			
			echo "$name initiated $(date)"
			java -Xmx10g -jar $tooldir/completeTCR.jar -pset flex -gene TRA $inputA $name.ccResult_alpha.txt
			mv descToID_mapping.txt $name.descToID_alpha.txt
			java -Xmx10g -jar $tooldir/completeTCR.jar -pset flex $inputB $name.ccResult_beta.txt
			mv descToID_mapping.txt $name.descToID_beta.txt	

			rsltA=`echo $name.ccResult_alpha.txt`
			mapA=`echo $name.descToID_alpha.txt`
			rsltB=`echo $name.ccResult_beta.txt`
			mapB=`echo $name.descToID_beta.txt`

			rPath=`which R`

			$rPath CMD BATCH --no-save "--args $outdir/$name $rsltA $mapA $rsltB $mapB $name" $tooldir/completeTCR.R $name.Rout
		
			if [[ $savetmp = "F" || -z $savetmp ]]
			then
				rm $name.descToID_alpha.txt
				rm $name.descToID_beta.txt
				rm $name.ctResult_alpha.txt
				rm $name.ctResult_beta.txt
				rm $name.Rout
			fi

			echo "$name complete: $(date)"
			cd ..
		else
			echo "ERROR: Sample name $name already exists. Please provide a unique sample name."
		fi	
done

