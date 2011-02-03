####################################################################
####################################################################
## Methods/Class Definitions
####################################################################
####################################################################


setClassUnion("AorN", c("array", "NULL"))

setClass("msSurv",representation(tree="graphNEL",
                                 ns="numeric",
                                 et="numeric",
                                 pos.trans="character",
                                 nt.states="numeric",
                                 dNs="array",
                                 Ys="array",
                                 ps="array",
                                 all.ajs="array",
                                 Fs="array",
                                 Gs="array",
                                 out="array",
                                 cov.p="array",
                                 sum.dNs="array",
                                 dNs.K="array",
                                 Ys.K="array",
                                 sum.dNs.K="array",
				 all.I_dA="array",
				 cov.dA="array",
				 Fs.var="AorN",
				 Gs.var="AorN"))
