# Electron Analysis Tree Maker

### Written by Bingchu Huang and Zach Miller

This code was designed to work within the STAR Collaboration Framework of
provided libraries, and is run within the ROOT Analysis program provided by
CERN. The goal of this code is to trim the data recorded by the STAR detector
into a manageable amount of data without losing any of the data that is
valuable for our particular analysis. A few notes before details:

* This builds on top of the work by the STAR collaborators at Lawerence
  Berkeley National Laboratory (LBNL). The StPicoDstMaker Class is written and
  provided by
    them, but is contained here because all of our classes reference it.

* The initial design of this code was done by Bingchu Huang with input from Zach Miller (and the UIC Group). The current version has been refactored many times by both of us and is maintained by both of us for different applications.

## Purpose

The amount of data that comes out of the STAR detector is HUGE. Even once it's
been encoded into Tuples, there are 100s of terrabytes to deal with. In this
form, the files are known as "MicroDSTs (muDSTs)." The muDSTs contain all of
the information from each event in the STAR detector... even the information
that is not relvant to our analysis. The LBNL group trimmed down this file
format even more to a smaller distribution known as a "PicoDST (pDST)."
This was done by removing data that is not necessary to any analyses, but that
was available in the muDST. The PicoDSTs are useful for analysis but still
contain much information that is not necessary for our analyses, and at a
precision that is prohibitive due to disk space (each pDST is 3+ GB). Thus,
the analysis tree maker sorts through the picoDST and only keeps events and
tracks that are useful, and only keeps variable for each of these classes
(event, track, etc) that are directly needed for analyses. The anaTreeMaker
also reduces the precision of the variables to only the needed values. For
example if we only need to be able to read 5 digits of a float value that runs
between 0-5... we multiply by 10000 and store it as a short integer, which
takes less memory than a pure float. The encoding and decoding of these values
is built into to the constructor and the "getVariable()" functions for each
class.

## How to Run (Z. Miller - July 2016)

### Reading Analysis Tree
If you are only going to read the analysis tree, an example analysis program
can be found in StMyAnaTreeMaker. This is invoked by macro/readAnaTree.C. In
general, the user will want to make many edits to the StMyAnaTreeMaker in
order to make the histograms of interest. However, there are many options that
can be set in readAnaTree.C. The most important of which is to tell the
AnaTreeReader what triggers you which to analyze... and what they are being
called. So if trigger 460201 is a BHT0 trigger, the user will call:
anaTreeMaker->addTrigger(460201,0). This syntax is further explained in the
comments of readAnaTree.C. A script as been provided to submit the tree
reading jobs to condor (the computing cluster) and is called
script/readAnaTreeSubmit.xml . This saves the relevant files and then runs
macro/readAnaTree.C, assuming that you point it to the correct file lists and
outputs a file of the histograms defined in StMyAnaTreeMaker.

Also useful is StMyAnaTreeMaker/histMaker.C (.h),
this de-couple histogram plotting and manipulation from the tree reading,
which allows quick plotting changes without needing to re-read the tree.
histMaker.C reads in the output of readAnaTree.C and grabs histograms, then
puts together canvases and plots the histograms on the canvases. The output of
this is set to be a single PDF file of all the canvases. This program should
be fairly easy to follow based on function names... but the main things are:

1. Load the output anaTree File
2. Get the Histograms (getHistograms() Function)
3. Manipulate the histograms (doProjections() function)
4. Prepare the various canvases for plots (prepareCanvas() function)
5. Put the plots on the canvas (various drawXYZ() functions)
6. make the PDF (makePDF() function)

Many other functions are provided to clean up the manipulation of the
histograms. histMaker.h contains many pre-defined const values, such as the pT
binning and the color palette selectors. Usage of histMaker.C is:

root -l -b -q 'histMaker.C("path/to/anaTree/Output.root")'

and will create a file called "path/to/anaTree/Output.root.pdf".


### Making Analysis Tree
The main portion if this lies in StPicoAnaTreeMaker/. The various classes
control what is saved for each branch (StEventHeader, StElectronTrack,
etc...), which is (usually) all done in the class constructor. The class
constructor uses the picoDst to get the information and store it in the
variable names, after properly converting it (if applicable) to a smaller data
size. This also makes it necessary to use the class functions to get the data
for analysis, as it re-converts it to the proper values (e.g. store 1.007 as 1007 to only use a short int, so when accessing you need to have it converted back to 1.007). 

To run this, there is a pre-made macro: macros/makeAnaTree.C. Various options
are set here using each classes "setting functions." The current work flow is
to read muDST, make the picoDST, then use the picoDST to make the anaTree...
then only save the anaTree. Future updates will allow use to read picoDst
directly if they become the collaboration standard. The key things to check
are that you are setting the Vz cuts correctly and have input the triggers you
wish to store (using the addTrigger functions). You must also make sure the
correct list of runs is fed to the analyzer so it knows which runs it should
process. The current runlist is called "runNumberList_run15pp_transverse", and is saved in the directory above StRoot. After you have this done, you
can submit using script/pp200.xml. If using the Purity Maker (option in the
makeAnaTree.C), you also need to make sure there is runList in
StPicoElecPurityMaker/ that contains the runs you want analyze.

## pico anaTree (B.Huang - March 2016)
main script and macro for submitting jobs starting from MuDst: script/AuAu200.xml and macros/makeAnaTree.C
Some options: 
prodMod: 0 for mb, 1 for ht, 2 for mtd (it's different in picoDstMaker, there is no "2" definedor different definition in picoDstMaker). It requires a run list as input. For Au+Au, a recenter file is required to correct event plane.
To save the recenter correction histograms: use macros/makeRecenter.C
All setters can be found in StPicoAnaTreeMaker.h 

submit jobs: script/submitAll.sh

hadd histograms and trees for each run: script/haddHist.sh and script/haddTree.sh
