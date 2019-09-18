# Implementation: *Molecular heterogeneity drives reconfigurable nematic liquid crystal drops*

These are custom Python codes that are associated with the free energy modeling in our paper, "Molecular heterogeneity induces reconfigurable nematic liquid crystal drops".

These codes can be easily followed and run, and will help interested readers verify our results and explore other combinations of parameters. 

The first program “ShapeTransitionCondition.py” calculates the required saddle-splay elastic constant, K24, that permits shape transitions based on user-defined parameters. User-defined parameters include the composite elastic constant K (or the splay K11 and bend K33 elastic constants, i.e., if the user prefers not to use one-constant approximation), the anchoring energy coefficient Wa, and the interfacial tension γ. By running the code with different parameters, the user can, for example, appreciate how reduction of γ helps push K24 into reasonable range for the shape transition. 

The second program “SegregationReducedEnergy.py” demonstrates how chain-length-dependent spatial segregation (i.e., into a long-chain-rich shell and a short-chain-rich core) lowers overall system elastic energy. A bidisperse demixing model was employed for this calculation, and the user can choose different weight ratios of monomer to macromer in order to interrogate different system mean chain lengths. Other tunable parameters include the volume ratio between shell and core, and the mean chain length near the drop surface.
