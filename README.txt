Understanding how multiple inputs may combine to create the output of a given target is a fundamental challenge in many fields, in particular in neuroscience. Shannon's Information Theory is the most suitable frame to cope with this problem: the interaction information decomposition (IID) indeed measures the balance between redundant and synergetic interaction within the classical multivariate entropy-based approach; the recent extension of Shannon's Information Theory leading to the Partial Information Decomposition (PID) led to provide specific measures for the information that several variables convey individually (unique information), redundantly (shared information) or only jointly (synergistic information) about the output.

The contribution of our work [1] is the proposal of an analytical frame where both IID and PID can be exactly evaluated across multiple timescales, for multivariate Gaussian processes, on the basis of simple vector autoregressive (VAR) identification. In doing this, our work opens the way for both the theoretical analysis and the practical implementation of information modification in processes that exhibit multiscale dynamical structures.

The Matlab toolbox presented here allows to compute analytically the parameters of a VAR model represented at multiple time scales, and to obtain from these parameters exact multiscale values for the information measures composing the IID and PID involving two drivers and a target inside a multivariate stochastic process. This is done exploiting the theory of state-space models, elaborating recent results provided in [2] and [3]. The reader is referred to the main paper [1] for details about the computation.

[1] Faes, L.; Marinazzo, D.; Stramaglia, S. Multiscale Information Decomposition: Exact Computation for Multivariate Gaussian Processes. Entropy 2017, in press.

[2]  Barnett, L.; Seth, A.K. Granger causality for state-space models. Phys. Rev. E 2015, 91, 040101.

[3] Solo, V. State-space analysis of Granger-Geweke causality measures with application to fMRI. Neural
Computation 2016, 28, 914–949.

---------------------------

The code is provided free of charge. It is neither exhaustively tested nor particularly well documented. The authors accept no liability for its use.
Use, modification and redistribution of the code is allowed in any way users see fit. Authors ask only that authorship is acknowledged and ref. [1] is cited upon utilization of the code in integral or partial form. 



To get started, we recommend that you run and work through the two demonstration scripts.

Luca Faes, Daniele Marinazzo, Sebastiano Stramaglia, August 2017.


Demonstration scripts
--------------------

test_simulation - Computes multiscale IID and PID measures for a simulated Gaussian 4-variate process, according to Sect. 4 in [1]

test_application - Performs VAR identification and computes multiscale IID and PID measures for an example of the public data used as application in [1]. (loads the data in example_data.mat)

Main computational functions
---------------------------

iss_PCOV    - Calculate partial variances from the innovations form state space parameters 

iss_ds    - Calculate innovations form state space parameters for a downsampled state space model

ss2iss    - Compute innovations form parameters for a general state space model by solution of a discrete algebraic Riccati equation (DARE)

varma2iss    - Compute parameters of an innovations form state space model from the parameters of the equivalent vector ARMA model

mos_idMVAR - Model Order Selection for identification of strictly causal VAR model

idMVAR - identification of strictly causal VAR model




NOTE: the iss_ds and ss2iss functions are taken from the State-Space Granger Causality Matlab® Toolbox  - http://users.sussex.ac.uk/~lionelb/downloads/ssgc.zip



Contacts
--------

Luca Faes
Healthcare Research and Innovation Program, FBK, Trento
BIOTech, Dept. of Industrial, Engineering, University of Trento
via delle Regole 101, 38123 Mattarello, Trento, Italy
Tel +39 0461 282773       Fax +39 0461 283659
faes.luca@gmail.com
www.lucafaes.net

Daniele Marinazzo
Department of Data Analysis
Faculty of Psychological and Educational Sciences
1, Henri Dunantlaan
B-9000 Ghent, Belgium
Phone: +32 (0) 9 264 6375
email: daniele.marinazzo@ugent.be
Twitter: dan_marinazzo

Sebastiano Stramaglia
Dipartimento di Fisica, Università degli Studi Aldo Moro, Bari, Italy
INFN, Sezione di Bari, Italy
email: sebastiano.stramaglia@ba.infn.it
