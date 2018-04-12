git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa

git cms-merge-topic guitargeek:ElectronID_MVA2017_V2_HZZ_940pre3
rm -rf RecoEgamma/ElectronIdentification/data
git clone -b ElectronID_MVA2017_V2 https://github.com/guitargeek/RecoEgamma-ElectronIdentification RecoEgamma/ElectronIdentification/data/
git clone -b 94X https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git
git cms-merge-topic cms-egamma:EGM_94X_v1
cd EgammaAnalysis/ElectronTools/data
git clone https://github.com/ECALELFS/ScalesSmearings.git
cd ScalesSmearings/
git checkout Run2017_17Nov2017_v1
cd ../../../../

git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v215 v2.1.5)
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
scram b -j 8

