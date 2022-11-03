# SeqMem
 Object oriented approach to processing statMatrix data for the Sequence Memory task
The objective here is to employ an object oriented approach to analyzing the statMatrix data. Though a very inefficient storage method, the statMatrix format is an intuitive organization scheme. In the future it will likely be phased out or reworked but for now it's what we're working with. This object oriented approach employs a three tiered class hierarchy. The main superclass (SeqMem here) compiles spiking and LFP data for a session into a single object and contains methods for organizing/extracting trial data as well as some generic processing methods (e.g. band-pass filtering & instantaneous phase calculation, manual artifact rejection UI coming in future updates). Intermediate classes (MLB_SM here) correspond to common analysis techniques and/or methods that can be applied in various ways with different approaches. Base subclasses (two at present I'm developing are PFC_TrialEvent_MLB_SM and PFC_UniSum_MLB_SM) represent specific analyses created to address specific questions about the data. For example, PFC_TrialEvent_MLB_SM and PFC_UniSum_MLB_SM are both implementations of the memoryless bayesian (MLB) algorithm but ..._TrialEvent_... is an ensemble analysis whereas ..._UniSum_... is a individual unit analysis.
