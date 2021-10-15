SCRAM_ARCH=slc7_amd64_gcc700; export SCRAM_ARCH
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src/
cmsenv
git cms-init
git clone -b UL_19 https://ferrico@github.com/ferrico/UFHZZAnalysisRun2.git
git cms-merge-topic asculac:Electron_XGBoost_MVA_16UL_17UL
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS
git cms-addpkg RecoEgamma/EgammaTools
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/
git cms-addpkg RecoEgamma/PhotonIdentification
git cms-addpkg RecoEgamma/ElectronIdentification
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools
git cms-addpkg RecoJets/JetProducers
git cms-addpkg PhysicsTools/PatAlgos/
git clone -b v2.3.5 https://github.com/JHUGen/JHUGenMELA
sh JHUGenMELA/MELA/setup.sh -j 8
git clone -b UL_19 https://github.com/ferrico/KinZfitter.git
scramv1 b -j 8





----------------------
OLD STUFF
----------------------
git cms-init
git cms-merge-topic mkovac:Electron_XGBoost_MVA_2016_and_2018_CMSSW_10_2_15
git cms-merge-topic cms-egamma:EgammaPostRecoTools 
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029
#git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd -
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
#scram b -j 8
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS
cp /raid/raid9/qguo/Run2/after/Run2_2/new/CMSSW_10_2_18/src/SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h 
cp /raid/raid9/qguo/Run2/after/Run2_2/new/CMSSW_10_2_18/src/GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc 
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa
#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v19 v1.9)

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v223 v2.2.3)
# replace ZZMatrixElement/MELA/setup.sh -j 8
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/COLLIER/
  pkgname="collier-1.2.0"
  pkgdir="COLLIER-1.2"
  tarname=$pkgname".tar.gz"
  tarweb="https://www.hepforge.org/archive/collier/"$tarname
  libname="libcollier.so"
  tmpdir="colliertmp"
  wget $tarweb
  mkdir $tmpdir
  tar -xvzf $tarname -C $tmpdir
  rm $tarname
  mv $tmpdir"/"$pkgdir"/src/"* ./
  rm -rf $tmpdir
  make
  mv $libname "../data/"$SCRAM_ARCH"/"$libname
popd
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
  make all
  mv libjhugenmela.so ../data/${SCRAM_ARCH}/
popd
popd
git clone -b UL_19 https://github.com/ferrico/KinZfitter.git #https://github.com/VBF-HZZ/KinZfitter.git
scram b -j 8

