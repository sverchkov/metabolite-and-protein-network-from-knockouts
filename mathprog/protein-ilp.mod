# A GNU MathProg formulation of the ILP of inferring the protein network
#
# Copyright Yuriy Sverchkov 2018

set P;
/* proteins */

param loglikP{ko in P, effect in P};
/* log-likelihoods (obrained from data) of correlated effects subject to knockouts */

param loglikM{ko in P, effect in P};
/* log-likelihoods (obtained from data) of anticorrelated effects subject to knockouts */

#var edgeP{pS in P, pT in P}, binary;
#/* the binary matrix of positive pS -> pT edges */

#var edgeM{pS in P, pT in P}, binary;
#/* the binary matrix of negative pS -> pT edges */

var pathPP{pS in P, pM in P, pT in P}, binary;
/* the binary tensor for the path pS -(+)-> pM -(+)-> pT */

var pathPM{pS in P, pM in P, pT in P}, binary;
/* the binary tensor for the path pS -(+)-> pM -(-)-> pT */

var pathMP{pS in P, pM in P, pT in P}, binary;
/* the binary tensor for the path pS -(-)-> pM -(+)-> pT */

var pathMM{pS in P, pM in P, pT in P}, binary;
/* the binary tensor for the path pS -(-)-> pM -(-)-> pT */

maximize loglikelihood: sum{ko in P, effect in P} ( loglikP[ko, effect] * pathPP[ko, ko, effect] + loglikM[ko, effect] * pathPM[ko, ko, effect] );

s.t. path_exclusivity{pS in P, pM in P, pT in P}: pathPP[pS,pM,pT] + pathPM[pS,pM,pT] + pathMP[pS,pM,pT] + pathMM[pS,pM,pT] <= 1;
s.t. self_reachable{p in P}: pathPP[p,p,p] = 1;
s.t. self_m_self_p_other_conflict{pS in P, pT in P}: pathMP[pS,pS,pT] = 0;
s.t. self_m_self_m_other_conflict{p1 in P, p2 in P}: pathMM[p1,p1,p2] = 0;
s.t. m_edge_collapse{p1 in P, p2 in P}: pathPM[p1,p1,p2] = pathMP[p1,p2,p2];
s.t. p_edge_collapse{p1 in P, p2 in P}: pathPP[p1,p1,p2] = pathPP[p1,p2,p2];

# Path logic for plus-plus paths
s.t. apbpc_needs_ab{p1 in P, p2 in P, p3 in P}: pathPP[p1,p2,p3] <= pathPP[p1,p2,p2];
s.t. apbpc_needs_bc{p1 in P, p2 in P, p3 in P}: pathPP[p1,p2,p3] <= pathPP[p2,p3,p3];
s.t. apbpc_implies_ab_bc{p1 in P, p2 in P, p3 in P}: pathPP[p1,p2,p3] >= pathPP[p1,p2,p2] + pathPP[p2,p3,p3] - 1;

# Path logic for plus-minus paths
s.t. apbmc_needs_ab{p1 in P, p2 in P, p3 in P}: pathPM[p1,p2,p3] <= pathPP[p1,p2,p2];
s.t. apbmc_needs_bc{p1 in P, p2 in P, p3 in P}: pathPM[p1,p2,p3] <= pathMP[p2,p3,p3];
s.t. apbmc_implies_ab_bc{p1 in P, p2 in P, p3 in P}: pathPM[p1,p2,p3] >= pathPP[p1,p2,p2] + pathMP[p2,p3,p3] - 1;

# Path logic for minus-plus paths
s.t. ambpc_needs_ab{p1 in P, p2 in P, p3 in P}: pathMP[p1,p2,p3] <= pathMP[p1,p2,p2];
s.t. ambpc_needs_bc{p1 in P, p2 in P, p3 in P}: pathMP[p1,p2,p3] <= pathPP[p2,p3,p3];
s.t. ambpc_implies_ab_bc{p1 in P, p2 in P, p3 in P}: pathMP[p1,p2,p3] >= pathMP[p1,p2,p2] + pathPP[p2,p3,p3] - 1;

# Path logic for minus-minus paths
s.t. ambmc_needs_ab{p1 in P, p2 in P, p3 in P}: pathMM[p1,p2,p3] <= pathMP[p1,p2,p2];
s.t. ambmc_needs_bc{p1 in P, p2 in P, p3 in P}: pathMM[p1,p2,p3] <= pathMP[p2,p3,p3];
s.t. ambmc_implies_ab_bc{p1 in P, p2 in P, p3 in P}: pathMM[p1,p2,p3] >= pathMP[p1,p2,p2] + pathMP[p2,p3,p3] - 1;

# Transitivity constraints
s.t. transitivity_pp{p1 in P, p2 in P, p3 in P}: pathPP[p1,p2,p3] <= pathPP[p1,p3,p3];
s.t. transitivity_pm{p1 in P, p2 in P, p3 in P}: pathPM[p1,p2,p3] <= pathMP[p1,p3,p3];
s.t. transitivity_mp{p1 in P, p2 in P, p3 in P}: pathMP[p1,p2,p3] <= pathMP[p1,p3,p3];
s.t. transitivity_mm{p1 in P, p2 in P, p3 in P}: pathMM[p1,p2,p3] <= pathPP[p1,p3,p3];

solve;

data;

# Load from file

set P := a b c d;

param : P : P : loglikP : loglikM :=
              a            a          1        -1
              a            b          1        -1
              a            c          0         0
              a            d          0         0
;

end;
