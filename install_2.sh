git cms-init
git clone -b FullRun_II https://github.com/VBF-HZZ/UFHZZAnalysisRun2.git
git clone https://github.com/VBF-HZZ/UFHZZAnalysisRun2-Accessary.git
mv UFHZZAnalysisRun2-Accessary/* ./
rm -rf UFHZZAnalysisRun2-Accessary
git clone https://github.com/VBF-HZZ/KinZfitter.git
scram b -j 8

