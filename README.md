# conditionalColoc

conditionalColoc is a software, written in Matlab, for analyzing the colocalization relationships between objects in three-channel microscopy images. 

Conceptually, it divides the channels into target, reference and condition, and addresses the question whether the extent of target colocalization with reference is influenced, positively or negatively, by target and/or reference colocalization with condition.

There are three main analysis functions:

1. A function that performs conditional colocalization analysis for the case of all channels containing punctate objects. In this case, the analysis can be performed with all permutations of target, reference and condition between the three channels. If punctate objects are detected with sub-pixel localization, the sub-pixel positions are used by the function for analysis.
2. A function that performs conditional colocalization analysis for the general case of punctate or non-punctate channels. Also in this case, the analysis can be performed with all permutations of target, reference and condition between the three channels. However, any sub-pixel positions of punctate objects will be rounded to integer pixel values.
3. A function that performs conditional colocalization analysis for the case where two channels contain punctate objects and only one channel contains non-punctate objects. In this case, the non-punctate channel may be the condition or reference, but not the target. The sub-pixel positions of punctate objects are retained in this function (if provided).

Being object-based, conditionalColoc rigorously assesses the statistical significance of any observed colocalization relationship through randomization and other statistical procedures.
