# conditionalColoc

conditionalColoc is a software, written in Matlab, for analyzing the colocalization relationships between objects in three-channel microscopy images. 

Conceptually, it divides the channels into target, reference and condition, and addresses the question whether the extent of target colocalization with reference is influenced, positively or negatively, with target and/or reference colocalization with condition.

It handles two types of three-channel images:

1. Images where all three channels contain punctate objects. In this case, the analysis can be performed with all permutations of target, reference and condition.
2. Images where two of the three channels contain punctate objects, while one contains non-punctate objects. In this case, the non-punctate channel may be the condition or reference, but not the target.

Being object-based, conditionalColoc assesses the statistical significance of any observed colocalization or influence on colocalization through randomization and other statistical procedures.
