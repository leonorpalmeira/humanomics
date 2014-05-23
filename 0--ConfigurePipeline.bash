#!/bin/bash

## Default parameters ##
output="pipeline.conf"
fastqsuffix="_R1.fastq.gz"
ref="hg19"
emailtype="END"
########################

## Functions ######
get_response() {
    echo -n $1
    read response
    if [ -n "$response" ]; then
	PARAM=$response
    fi
}

check_file() {
    if [ ! -f $1 ] ; then
	echo "File does not exist."
	exit 1
    # else
    # 	echo "Using file: [$1]"
    fi
}

check_output() {
    if [ -f $1 ]; then
	echo -n "Output file exists. Overwrite? [y/n] > "
	read response
	if [ "$response" != "y" ]; then
            echo "Exiting program."
            exit 1
	# else
	#     echo "Overwriting file: [$output]"
	fi
    fi
}

check_directory() {
    if [ ! -d $1 ] ; then
	echo "Directory does not exist."
	echo "Creating directory."
	mkdir -p $1
    elif [ "$1" == "" ] ; then
	echo "Please provide a directory."
	echo "Exiting program."
	exit 1
    # else
    # 	echo "Using existing directory: [$1]"
    fi
}
###################

unset PARAM
get_response "Enter name of output configuration file [$output] > "
if [ -n "$PARAM" ]; then
    output=$PARAM
fi
check_output $output

unset PARAM
get_response "Enter the full path of the raw data directory [no default] > "
if [ -n "$PARAM" ]; then
    datadir=$PARAM
fi
check_directory $datadir

unset PARAM
get_response "Enter the full path of the directory for pipeline results [no default] > "
if [ -n "$PARAM" ]; then
    resdir=$PARAM
fi
check_directory $resdir

unset PARAM
get_response "Enter name of raw Fastq.gz suffix (paired-end expected) [$fastqsuffix] > " # WTF?? Very strange behaviour when a file named 'a' is present in the same directory as this script...???
echo $fastqsuffix
if [ -n "$PARAM" ]; then
    fastqsuffix=$PARAM
fi
# echo "Using suffix: [$fastqsuffix]"

unset PARAM
get_response "Enter name of reference assembly [$ref] > "
if [ -n "$PARAM" ]; then
    ref=$PARAM
fi
if [ "$ref" != "hg19" ]; then
    echo "Only hg19 is implemented yet. For other assemblies, you should change the pipeline configuration file by hand and at your own risk."
fi
# echo "Using suffix: [$ref]"

unset PARAM
get_response "Enter full path of target file [no default] > "
if [ -n "$PARAM" ]; then
    targets=$PARAM
fi
check_file $targets

unset PARAM
get_response "Enter baitNames [no default] > "
if [ -n "$PARAM" ]; then
    baitNames=$PARAM
fi

unset PARAM
get_response "Enter full path of baits (Picard) file [no default] > "
if [ -n "$PARAM" ]; then
    baitsPicard=$PARAM
fi
check_file $baitsPicard

unset PARAM
get_response "Enter full path of target (Picard) file [no default] > "
if [ -n "$PARAM" ]; then
    targetsPicard=$PARAM
fi
check_file $targetsPicard

unset PARAM
get_response "Email address for SLURM [no default] > "
if [ -n "$PARAM" ]; then
    email=$PARAM
fi

unset PARAM
get_response "Email type for SLURM [$emailtype] > "
if [ -n "$PARAM" ]; then
    emailtype=$PARAM
fi

unset PARAM
get_response "Log directory for SLURM [no default] > "
if [ -n "$PARAM" ]; then
    slurmlogdir=$PARAM
fi
check_directory $slurmlogdir

echo "Estimating array value for SLURM ..."
slurmarray="1-"`\ls $datadir/*R1* | wc -l`
echo "Done."


# # Read Groups should follow:

# # They are assumed to have been run on a hiseq





## Preparing output ##
str="Raw data directory\t$datadir\n"
str+="Results directory\t$resdir\n"
str+="Fastq.gz suffix (paired-end expected)\t$fastqsuffix\n"
str+="Reference assembly\t$ref\n"
str+="Target file\t$targets\n"
str+="baitNames file\t$baitNames\n"
str+="Baits file (Picard)\t$baitsPicard\n"
str+="Target file (Picard)\t$targetsPicard\n"
str+="SLURM Email address\t$email\n"
str+="SLURM Email type\t$emailtype\n"
str+="SLURM Log directory\t$slurmlogdir\n"
str+="SLURM array\t$slurmarray\n"

echo ""
echo "######### Pipeline configuration ########"
echo -en $str
echo "#########################################"
echo ""

echo -en $str > $output
######################

echo "This pipeline configuration has been written to: '$output'"
echo "WARNING: currently, Read Groups have to be modified manually"
exit 1
