# first parameter: subject ID for the desired trial (e.g. control/CF1)
# second parameter: ibifolder selection based on the method 0-> IBIHEALPEAK 1-> ibi.txt
# third parameter: 0 -> EMBC method 1 -> new HEAL-T method
# fourth  parameter: sampling frequency (FS) value (synchronize BVP with Accel) [Hz].
# fifth parameter: Fp1 (Hz) low cut-off first filter BHW [Hz].
# sixth parameter: Fz1 (Hz) high cut-off first filter BHW [Hz].
# seventh parameter: Fp2 (Hz) low cut-off second filter BHW [Hz].
# eight parameter:  Fz2 (Hz) high cut-off second filter BHW [Hz].
# 9th parameter: inc1 filter one (Hz)
# 10th parameter: inc2 filter two (Hz)
# 11th parameter: overlap per each BVP segment (sec or percentage)
# 12th parameter: windows size (WIN) BVP segment (sec or number of windows)
# 13th parameter: spline type selector 1-> Cubic Spline, 0-> Moving Average Filter
# 14th parameter: name of the nohup output file that will be save in the /log folder
# 15th parameter: MATLAB_BINARY "matlab binary execution path defined by the user (i.e, in Matlab /usr/local/../bin/matlab and in mac /Applications/Matlab.../bin/matlab)"

ppwd=$(pwd)
pid=$(echo $$)
echo $ppwd

if [ $# -eq 4 ]; then
   echo "The number of parameters is not enough for the execution, please refer to readme to define the parameters accordingly."
fi

if [ $# -eq 14 ]; then
   echo "Defining default input and output folders for HR execution..."
   ## Default input and output folders
   input="$ppwd/example_data/"
   output="$ppwd/example_data/out_dir/"
else  ## this part of the code define the output and input
   input="$1"
   output="$1"
fi

if [ $# -eq 15 ]; then
   #rm "$ppwd/log/$pid-${14}"
   echo "Selecting user-defined input and output folders for HR execution..."
   input="$1"
   output="$1"
fi

if [ $# -eq 16 ]; then
    #rm "$ppwd/log/$pid-${15}"
    echo "Selecting user-defined input and output folders for HR execution..."
    input="$1"
    output="$2"
fi
if [ -f "$output/ibiIBIselHEAL.txt" ] && [ "$2" == "1" ]; then
    echo "The current folder is already updated, Do you want to run it anyways? (yes/no)"
    read indicator
    echo $indicator
    if [ "$indicator" == "no" ]; then
        echo "Quitting..."
        exit 1
    fi
fi
if [ -f "$output/iibibase.txt" ] && [ "$2" == "0" ]; then
    echo "The current folder is already updated, Do you want to run it anyways? (yes/no)"
    read  indicator
    if [ "$indicator" == "no" ]; then
        echo "Quitting..."
	exit 1
    fi
fi
echo "Running..."
echo "$input  $output"
if [ $# -eq 15 ]; then
    nohup "${15}" -r "heal_t_call('$input','$output',$2,$3,$4,{[$5 $6],[$7 $8]},[$9 ${10}],${11},${12},${13}); exit;" &> "$ppwd/log/$pid-${14}" &
fi
if [ $# -eq 16 ]; then
    nohup "${16}" -r "heal_t_call('$input','$output',$3,$4,$5,{[$6 $7],[$8 $9]},[${10} ${11}],${12},${13},${14}); exit;" &> "$ppwd/log/$pid-${15}" &
fi
if [ $# -eq 14 ]; then
    nohup "${14}" -r "heal_t_call('$input','$output',$1,$2,$3,{[$4 $5],[$6 $7]},[$8 $9],${10},${11},${12}); exit;" &> "$ppwd/log/$pid-${13}" &
fi
#echo "'$1' $2 $3 $4 $5"
