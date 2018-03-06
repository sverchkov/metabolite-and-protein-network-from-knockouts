# A GNU MathProg formulation of the ILP of inferring the protein network
#
# Copyright Yuriy Sverchkov 2018

set P;
/* proteins */

set E := plus minus;
/* edge types */

param loglik{ko in P, resp in E, effect in P};
/* log-likelihoods (obrained from data) of response types given effects subject to knockouts */

var t{pS in P, e1 in E, pM in P, e2 in E, pT in P}, binary;
/* the binary tensor for the path pS -(e1)-> pM -(e2)-> pT */

maximize loglikelihood: sum{ko in P, resp in E, effect in P} loglik[ko, resp, effect] * t[ko, resp, effect, plus, effect];

s.t. self-reachable{p in P}: t[p,plus,p,plus,p] = 1;
s.t. self-conflict-start{p1 in P, e in E, p2 in P}: t[p1,minus,p1,e,p2] = 0;
s.t. self-conflict-end{p1 in P, e in E, p2 in P}: t[p1,e,p2,minus,p2] = 0;
s.t. edge-collapse{p1 in P, e in E, p2 in P}: t[p1,plus,p1,e,p2] = t[p1,e,p2,plus,p2];
s.t. transitivity-first-edge{pS in P, e1 in E, pM in P, e2 in E, pT in P}: t[pS,e1,pM,e2,pT] <= t[pS,e1,pM,plus,pM];
s.t. transitivity-last-edge{pS in P, e1 in E, pM in P, e2 in E, pT in P}: t[pS,e1,pM,e2,pT] <= t[pM,plus,pM,e2,pT];
s.t. transitivity-inference{pS in P, e1 in R, pM in P, e2 in R, pT in P}: t[pS,e1,pM,e2,pT] >= t[pS,e1,pM,plus,pM] + t[pM,plus,pM,e2,pT] - 1;

