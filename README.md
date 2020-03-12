# FNT2
Fast Neutron Tomography sort code

It compiles on the command line (although that will be slowest):

 root -l src/Helper.cpp+ src/Gate.cpp+ src/Channel.cpp+ src/FNT.cpp+ src/analysis.cpp+ fntsort.cpp+

If you see a floating point exception, make sure to check the files in line 10-15 of FNT.cpp actually exist, and the same for the root files in the file called by fileFiles.

If you std::out_of_range error, check the histograms all exist in histos.txt and remember that non-channel histograms end in 0 (= 'all channels') after the name listed



Example output:

root [0]

Processing src/Helper.cpp+...

Info in <TUnixSystem::ACLiC>: creating shared library /mnt/d/fnt/macros/git/FNT2/./src/Helper_cpp.so

Processing src/Gate.cpp+...

Info in <TUnixSystem::ACLiC>: creating shared library /mnt/d/fnt/macros/git/FNT2/./src/Gate_cpp.so

Processing src/Channel.cpp+...

Info in <TUnixSystem::ACLiC>: creating shared library /mnt/d/fnt/macros/git/FNT2/./src/Channel_cpp.so

Processing src/FNT.cpp+...

Info in <TUnixSystem::ACLiC>: creating shared library /mnt/d/fnt/macros/git/FNT2/./src/FNT_cpp.so

Processing src/analysis.cpp+...

Info in <TUnixSystem::ACLiC>: creating shared library /mnt/d/fnt/macros/git/FNT2/./src/analysis_cpp.so

Processing fntsort.cpp++...

Info in <TUnixSystem::ACLiC>: creating shared library /mnt/d/fnt/macros/git/FNT2/./fntsort_cpp.so

Tree output is in out.root and histograms will be saved in histograms.root

No tree exists!  Creating...

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    100%

There are 240834975 entries, with a minimum timestamp of 1188244930856500 (1e16) and a maximum timestamp of 2309649679480218 (2e16)

Analysis object 0x7fffc6059420 has compiled successfully, continuing bespoke analysis from analysis.cpp...

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    100%

Moving histograms to timeCfolder ...

Moving histograms to nrjCfolder ...

Moving histograms to nrj2Cfolder ...

Moving histograms to time_adjCfolder ...

Moving histograms to nrj_adjCfolder ...

Moving histograms to nrj2_adjCfolder ...

Analysis complete!

root [5]

 
