git cms-init
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS
cp /raid/raid7/dsperka/Run2/HZZ4l/CMSSW_10_2_5/src/SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h SimDataFormats/HTXS/interface/
cp /raid/raid7/dsperka/Run2/HZZ4l/CMSSW_10_2_5/src/GeneratorInterface/RivetInterface/src/HiggsTemplateCrossSections.cc GeneratorInterface/RivetInterface/src/
cp /raid/raid7/dsperka/Run2/HZZ4l/CMSSW_10_2_5/src/GeneratorInterface/RivetInterface/plugins/HTXSRivetProducer.cc GeneratorInterface/RivetInterface/plugins/
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa
git clone -b 102X https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v218 v2.1.8)
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
git clone https://github.com/VBF-HZZ/KinZfitter.git
scram b -j 8

