# MajProj_ALightweightMedicalImageEncryptionSchemeUsingChaoticMapsAndImageScrambling_SujithaSudevan
This repository contains the code implemented for the project "A Lightweight Medical Image Encryption Scheme Using Chaotic Maps and Image Scrambling"

## A Lightweight Medical Image Encryption Scheme Using Chaotic Maps and Image Scrambling

This repository is divided into 2 parts:
1. **/src** folder: Contains 3 folders /algorithms, /dataset, and /security analysis and the proposed algorithm in MATLAB.
2. **Architectural Block Diagram**: Contains the pictorial representation of the flow of the code implementation.

**/src** folder contains these sub-folders:
1. **/algorithms**: These are the recent papers used for comparison in performance and security analysis with the proposed algorithm.
	
	A. Implementation of the paper "**A Novel Bit-level Image Encryption Algorithm Based on Chaotic Maps**" in Matlab.
	
	Link for paper: https://www.sciencedirect.com/science/article/abs/pii/S0143816615002109
	
	**pwlcm_algorithm.m**
	
	B. Implementation of the paper "**A Novel Chaotic Image Encryption Algorithm Based on Improved Baker Map and Logistic Map**" in Matlab.
	
	Link for paper: https://link.springer.com/article/10.1007/s11042-019-7453-3
	
	**bakers2d_1dlogistic_algorithm.m**
	
	C. Implementation of the paper "**Medical Image Encryption Scheme Using Multiple Chaotic Maps**" in Matlab.
	
	Link for paper: https://www.sciencedirect.com/science/article/abs/pii/S0167865521003913
	
	**arnoldscat_2dlscm_algorithm.m**

2. **/dataset**: Contains different images of different pixel sizes. Also contains Image Description (.txt file) which contains details about all the images used for implementation.

3. **/security analysis**: Contains code for implementation of security analysis in Matlab
	
	A. Pseudo Code for Performance analysis (.txt files)
		
		a. end-to-end runtime analysis.txt
		
		b. end-to-end memory analysis.txt
	
	B. Implementation of security analysis in Matlab
		
		a. key_sensitivity_analysis.m
		
		b. histogram_analysis.m
		
		c. shannons_entropy_analysis.m
		
		d. differential_cryptanalysis.m
		
		e. contrast_analysis.m

4. Implementation of the proposed algorithm "**A Lightweight Medical Image Encryption Scheme Using Chaotic Maps and Image Scrambling**" in Matlab.
	
	**proposed_algorithm.m**
