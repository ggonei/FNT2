# FNT2
Fast Neutron Tomography sort code

It compiles on the command line (although that will be slowest):

 root -l src/Helper.cpp+ src/Gate.cpp+ src/Channel.cpp+ src/FNT.cpp+ src/analysis.cpp+ fntsort.cpp+

# Troubleshooting
If you see a floating point exception, make sure to check the files in line 10-15 of FNT.cpp actually exist, and the same for the root files in the file called by fileFiles

If you std::out_of_range error, check the histograms all exist in histos.txt and remember that non-channel histograms end in 0 (= 'all channels') after the name listed

If you get a <TBufferFile::WriteByteCount> error, or a segmentation fault when newHists is written, then one of your 2-D histograms is probably too big

