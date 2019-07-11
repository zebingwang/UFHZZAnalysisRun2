git cms-init
git cms-merge-topic mkovac:Electron_XGBoost_MVA_2016_and_2018_CMSSW_10_3_1
git cms-merge-topic cms-egamma:EgammaPostRecoTools 
#git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd -
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
#scram b -j 8
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS
cp /raid/raid7/dsperka/Run2/HZZ4l/CMSSW_10_2_5/src/SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h SimDataFormats/HTXS/interface/
cp /raid/raid7/dsperka/Run2/HZZ4l/CMSSW_10_2_5/src/GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc GeneratorInterface/RivetInterface/src/
cp /raid/raid7/dsperka/Run2/HZZ4l/CMSSW_10_2_5/src/GeneratorInterface/RivetInterface/plugins/HTXSRivetProducer.cc GeneratorInterface/RivetInterface/plugins/
git clone https://github.com/mkovac/MuonMVAReader.git MuonMVAReader
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa
git clone -b 2018_L https://github.com/qyguo/UFHZZAnalysisRun2.git
#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v14 v1.4)

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v220 v2.2.0)
# replace ZZMatrixElement/MELA/setup.sh -j 8
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/COLLIER/
  pkgname="collier-1.2"
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
git clone https://github.com/VBF-HZZ/KinZfitter.git
scram b -j 8

