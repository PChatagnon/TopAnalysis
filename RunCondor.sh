#!/bin/zsh
outpath=/eos/user/p/pchatagn/OutputDirFPMC/MINIAOD
N=300

extime="tomorrow" #"testmatch" #testmatch tomorrow workday
condor="condor_generator.sub"
echo "executable  = run.sh" > $condor
echo "output      = ${condor}.out" >> $condor
echo "error       = ${condor}.err" >> $condor
echo "log         = ${condor}.log" >> $condor
echo "+JobFlavour =\"${extime}\"">> $condor
echo "requirements = (OpSysAndVer =?= \"CentOS7\")" >> $condor  # SLCern6 CentOS7

TYPINC='QCD' # QED QCD
NFLUX='9'
folder='exclu_ttbar'
command='11706'
for mode in lephad; do #hadlep leplep; do
for i in {1..${N}}; do
    if [ "$mode" = 'lephad' ]; then
      ii=$(printf "%02d" `expr ${i} + 0`)
    elif [ "$mode" = 'hadlep' ]; then
      ii=$(printf "%02d" `expr ${i} + ${N}`)
    elif [ "$mode" = 'leplep' ]; then
      ii=$(printf "%02d" `expr ${i} + ${N} + ${N}`)
    fi
    if test ! -f "${outpath}/${folder}_${mode}_${TYPINC}/${folder}_${mode}_${TYPINC}_${ii}.root"; then
      echo submit ${outpath}/${folder}_${mode}_${TYPINC}/${folder}_${mode}_${TYPINC}_${ii}.root
      echo "arguments   = ${folder} ${command} ${ii} ${TYPINC} ${NFLUX} ${mode}" >> $condor
      echo "queue 1" >> $condor
    fi
done
done

echo "Submitting $condor"
echo condor_submit $condor
condor_submit $condor
