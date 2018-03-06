* A gams formalization of the protein matrix inference ILP

Set
   protein 'proteins'
      / a
        b
        c
        d /
   edge 'edge types'
      / pos 'positive edge'
        neg 'negative edge' /;

Alias(protein,pm,pl)
Alias(edge,edgep)

Parameter loglik(protein,protein) "ko table data";
Table loglik(protein,protein)
  a b c d
a 1 1 1 1
b 1 1 1 1

#TODO

