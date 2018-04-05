##############################
## Load reactions from KEGG ##
##############################

library( "KEGGREST" )

# hsa is the Human (Homo sapiens) gene database

genes.by.enzymes <- keggLink( "hsa", "enzyme" )
enzymes.by.reaction <- keggLink( "enzyme", "reaction" )
compound.by.reaction <- keggLink( "compound", "reaction" )


