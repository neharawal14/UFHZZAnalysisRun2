echo "================= PSet.py file =================="
cat PSet.py

# fix MCFM issue
echo "================= Fix MCFM issue ================"
echo "$PWD"
rm *.dat
rm *.DAT
rm br.sm*
rm -rf Pdfdata
cp -r ${CMSSW_BASE}/src/JHUGenMELA/MELA/data/Pdfdata ${CMSSW_BASE}/src/
MCFM_LIBS_PATH=${CMSSW_BASE}/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCFM_LIBS_PATH}
echo ${LD_LIBRARY_PATH}
echo "================= DONE =========================="


# Copy Input Files
echo "================= Copy Input Files ================"
echo "${PWD}"

#rm -rf ${CMSSW_BASE}/src/ZZMatrixElement/MEKD
#mv ${CMSSW_BASE}/MEKD ${CMSSW_BASE}/src/ZZMatrixElement/

#rm -rf ${CMSSW_BASE}/src/ZZMatrixElement/MELA
#mv ${CMSSW_BASE}/MELA ${CMSSW_BASE}/src/ZZMatrixElement/

rm -rf ${CMSSW_BASE}/src/JHUGenMELA/MELA 
mv ${CMSSW_BASE}/MELA ${CMSSW_BASE}/src/JHUGenMELA/
cd ${CMSSW_BASE}/src/JHUGenMELA/MELA
eval $(./setup.sh env)
cd -

rm -rf ${CMSSW_BASE}/src/KinZfitter/KinZfitter/ParamZ1/
mv ${CMSSW_BASE}/ParamZ1 ${CMSSW_BASE}/src/KinZfitter/KinZfitter/

rm -rf ${CMSSW_BASE}/src/KinZfitter/HelperFunction/hists/
mv ${CMSSW_BASE}/hists ${CMSSW_BASE}/src/KinZfitter/HelperFunction/

echo "================= DONE =========================="

cmsRun -j FrameworkJobReport.xml -p PSet.py 
