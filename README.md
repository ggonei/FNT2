# FNT2
Fast Neutron Tomography sort code

It compiles on the command line (although that will be slowest):

 root -l src/Helper.cpp+ src/Gate.cpp+ src/Channel.cpp+ src/FNT.cpp+ src/analysis.cpp+ fntsort.cpp+

# Speed
unordered_map is a hash table and faster than map here but does require C11, so switch to map on old systems

No compression is quicker by ~20% but the file size is 10x larger, so is NOT recommended

# Troubleshooting
If you see a floating point exception, make sure to check the files in line 10-15 of FNT.cpp actually exist, and the same for the root files in the file called by fileFiles

If you std::out_of_range error, check the histograms all exist in histos.txt and remember that non-channel histograms end in 0 (= 'all channels') after the name listed, remove the template code to find the issue

If you get a <TBufferFile::WriteByteCount> error, or a segmentation fault when newHists is written, then one of your 2-D histograms is probably too big

If your entire environment crashes with absolutely no error or useful messages even in linked debuggers you probably ran out of memory - decrease number of bins in histograms especially multi-D and per channel histograms and those where peak finder is used

If after the above you still have errors, did you delete out.root?  It might have been made before the crash occurred, and so the program is now using that corrupt file.  Delete out.root and try again.
