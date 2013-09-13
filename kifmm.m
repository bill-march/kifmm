% Main function for kernel-independent fmm
function [] = kifmm(Sources, Targets, Kernel)



    Souce_Tree = BuildTree(Sources);
    Target_Tree = BuildTree(Targets);
    
    
