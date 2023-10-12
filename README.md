# LocalDiffusion
Code with JBM

The files within this repository were the state of the art for inferenceMap that we used when performing the analysis in Obara et al. Nature 2023. 
The directory "LocalDiffusion_MatlabScripts" contains everything needed to reproduce all the data in the paper (and this is the code we used to perform everything). 
The directory "CompiledStandAloneVersions" contains code to allow the same analysis by someone without a matlab license. These are dated as of the time we used them, and may not be forward compatible with all operating systems. If you have errors, try the matlab source code or the updated versions (see below). 
The directory "ChrisWrapperCode" is matlab scripts to ease the conversion of the outputs from either of the above packages into the structures and analysis pipelines used in the remainder of the paper. Email Chris at cjobara@gmail.com if you are trying to do this and he'll try to help you...it's a bit confusng.

Note: There are more updated versions of local diffusion analysis with inferenceMap that have been updated since these. They will have less bugs and are supported on modern operating systems, we provide this code here for archival purposes. You can download the updated versions and appropriate documentation for free at the following links:
https://pypi.org/project/tramway/
https://github.com/DecBayComp/TRamWAy

