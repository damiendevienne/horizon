require(ape)
require(reshape2)
rm(list=ls())
###
CreateSpeciesGroup<-function(tree, nbgroup=100) {
	sps_per_group<-Ntip(tree)/nbgroup
	nbdesc<-Descendants(tree, (Ntip(tree)+1):(Ntip(tree)+Nnode(tree)))
	size_nbdesc<-unlist(lapply(nbdesc, length))
	possible_nodes<-which(abs(size_nbdesc-sps_per_group)==min(abs(size_nbdesc-sps_per_group)))
	node2rm<-sample(possible_nodes, 1)
}

### read tables
hp<-read.table("../input/hp_links.txt", header=T)
hpl<-read.table("../input/host_and_plants_matrix.tab") ##hosts and plants
H<-read.table("../input/hosts.txt", header=T, row.names=1)
P<-read.table("../input/parasitoids.txt", header=T,row.names=1)
hh<-read.table("../input/host_links_by_plant_sp.txt", header=T)

##read trees
host_tree<-read.tree("../input/host_tree_JUNE2019.tre")
par_tree<-read.tree("../input/par_tree_JUNE2019.tre")

##clean trees and tables
#clean parameters
S_col_events<-3
good_H <- H[which(H$n_col_event>=S_col_events),]
##NEW: remove from H those species absent in the tree:
good_H <- good_H[is.element(rownames(good_H),host_tree$tip.label),]
#good_H <- good_H[which(is.na(good_H$subfamily)==F),]
good_P <- P[which(P$n_col_event>=S_col_events),]
H<-good_H
P<-good_P
#good_HP <- hp[is.element(hp$host,row.names(good_H)) & is.element(hp$par,row.names(good_P)),]
host_tree<-keep.tip(host_tree, rownames(H))
par_tree<-keep.tip(par_tree, rownames(P))

##
## WRITE THESE TREES AND COMPUTE THE GROUPS FOR VARIOUS 
## VALUES OF N (nb of groups)
## ALL THIS BLOCK NEEDS TO BE DONE ONLY ONCE.
## SO I COMMENT IT.
##

# write.tree(host_tree, file="../output/host_tree_final_clean")
# write.tree(par_tree, file="../output/par_tree_final_clean")
# for (i in seq(50,150, by=10)) {
#     print(i)
#     system(paste("python3 MakeGroupsInTree.py ../output/host_tree_final_clean ",i," > ../output/host_groups_N_",i,sep=""))
#     system(paste("python3 MakeGroupsInTree.py ../output/par_tree_final_clean ",i," > ../output/par_groups_N_",i,sep=""))
# }

##
## END OF THE BLOCK
##


### Get DISTANCES in trees
Dh<-cophenetic(host_tree)
Dp<-cophenetic(par_tree)

### READ CLUSTERS FROM EXTERNAL FILES, GIVEN N
### AND COMPUTE THE EXHAUSTIVE LIST OF CUSTER MEMBERSHIP
##1. for hosts
N_host<-100 #nb de clusters hôtes
Host_Clust<-read.table(paste("../output/host_groups_N_",N_host,sep=""))
Clhost_tree_init<-Host_Clust$V2
names(Clhost_tree_init)<-Host_Clust$V1
##ADD missing species (not in a cluster) as members of their own clusters
missing_hosts<-host_tree$tip.label[!is.element(host_tree$tip.label, names(Clhost_tree_init))]
missing_groups_hosts<-(max(Clhost_tree_init)+1):(max(Clhost_tree_init)+length(missing_hosts))
names(missing_groups_hosts)<-missing_hosts
Clhost_tree<-c(Clhost_tree_init, missing_groups_hosts) ##Correct list of groups

##2. for parasites
N_par<-100 #nb de clusters parasites
Par_Clust<-read.table(paste("../output/par_groups_N_",N_host,sep=""))
Clpar_tree_init<-Par_Clust$V2
names(Clpar_tree_init)<-Par_Clust$V1
##ADD missing species (not in a cluster) as members of their own clusters
missing_pars<-par_tree$tip.label[!is.element(par_tree$tip.label, names(Clpar_tree_init))]
missing_groups_par<-(max(Clpar_tree_init)+1):(max(Clpar_tree_init)+length(missing_pars))
names(missing_groups_par)<-missing_pars
Clpar_tree<-c(Clpar_tree_init, missing_groups_par) ##Correct list of groups




GetTriplets<-function(S_lien=5,S_nl=10, D_seuil_min=0.1) {	
	hp_all<-table(hp[,1:2])
	hp_l<-hp[(hp$n_links>S_lien),]
	
	hp_lien<-hp_all
	hp_lien[cbind(as.character(hp_l$host), as.character(hp_l$par))]<-hp_lien[cbind(as.character(hp_l$host), as.character(hp_l$par))]+1
	hp_lien<-hp_lien-1
	hp_lien[hp_lien<0]<-0
	
	P_nl<-rownames(P[P$n_obs_tot>S_nl,])
	H_nl<-rownames(H[H$n_obs_tot>S_nl,])
	hp_nl1<-hp_all
	hp_nl1[hp_nl1>0]<-0 #initialize
	hp_nl1[H_nl,P_nl]<-1 # Cells with a one are those for which we have strong confidence 
	hp_nl2<-hp_all
	hp_nl2<-1-hp_nl2 #0 means that there is a link (even soft); 1 that there is not.
	hp_nonlien<-hp_nl1*hp_nl2 #because we want those cells that have a 1 in nl1 (strong confidence) and a 1 in nl2 (this strong non link is absent )
	
	###Then we subsample those host and parasites that DO have strong links AND strong non-links
	eligible_H<-rownames(hp_all)[(apply(hp_lien,1,sum)>0)&(apply(hp_nonlien,1,sum)>0)]
	eligible_P<-colnames(hp_all)[(apply(hp_lien,2,sum)>0)&(apply(hp_nonlien,2,sum)>0)]
	all_H<-rownames(hp_all)
	all_P<-colnames(hp_all)
	
	print ("Create H1H2, P1P2 pairs...")
	####Create a list of H1H1 pairs with ok taxonomy (same subfamily but d>D_seuil_min)
	Dhsmall<-Dh[eligible_H,eligible_H]
	Dhsmall[upper.tri(Dhsmall)] <- NA
	diag(Dhsmall)<-NA
	H1H2<-melt(Dhsmall, na.rm=T)	# remove links involving species with not enough collecting events (and host species with not subfamily)
	# add cluster of each partner (TOREDO)
	H1H2<-cbind(H1H2, gr1=Clhost_tree[as.character(H1H2$Var1)], gr2=Clhost_tree[as.character(H1H2$Var2)])
	H1H2<-H1H2[(H1H2$value>=D_seuil_min)&(H1H2$gr1==H1H2$gr2),]


	Dpsmall<-Dp[eligible_P,eligible_P]
	Dpsmall[upper.tri(Dpsmall)] <- NA
	diag(Dpsmall)<-NA
	P1P2<-melt(Dpsmall, na.rm=T)
	# add cluster of each partner (TOREDO)
	P1P2<-cbind(P1P2, gr1=Clpar_tree[as.character(P1P2$Var1)], gr2=Clpar_tree[as.character(P1P2$Var2)])
	P1P2<-P1P2[(P1P2$value>=D_seuil_min)&(P1P2$gr1==P1P2$gr2),]
    ##Now we will test all of the possible pair of pairs
    ##For All Pairs of H1H2 we test against all eligible_P
    TheTripletTest_H<-function(arr) {
      ##i and j are the references of the rows of H1H2 and P1P2 respectively
      H1<-arr[1]
      H2<-arr[2]
      sp_forming_triplets<-names(which(apply(hp_nonlien[c(H1,H2),all_P]-hp_lien[c(H1,H2),all_P],2,prod)=="-1"))
      if (length(sp_forming_triplets)>0) return(cbind(H1,H2,sp_forming_triplets))
      else return(NA)
    }
    TheTripletTest_P<-function(arr) {
      ##i and j are the references of the rows of H1H2 and P1P2 respectively
      P1<-arr[1]
      P2<-arr[2]
      sp_forming_triplets<-names(which(apply(hp_nonlien[all_H, c(P1,P2)]-hp_lien[all_H, c(P1,P2)],1,prod)=="-1"))
      
      if (length(sp_forming_triplets)>0) return(cbind(P1,P2,sp_forming_triplets))
      else return(NA)
    }
    print("Get Host Triplets...")
    Trip_H<-apply(H1H2,1, TheTripletTest_H)
    print("Get Parasite Triplets...")
    Trip_P<-apply(P1P2,1, TheTripletTest_P)
    
    Trip_H<-Trip_H[!is.na(Trip_H)]
    Trip_P<-Trip_P[!is.na(Trip_P)]
    TRIPH<-do.call(rbind, Trip_H)
    TRIPP<-do.call(rbind, Trip_P)
    rownames(TRIPH)<-NULL
    rownames(TRIPP)<-NULL
    return(list(TRIPH=TRIPH,TRIPP=TRIPP))
}

GetTripletPlants<-function(S_lien=5,S_nl=10, D_seuil_minH1H2=0.1) {
  print("Get Plant Triplets...")
  # for each pair of hosts, we compute their link and their non-link
  # according to plant sharing. We will then use that to make triplets of interest
  Hnames<-rownames(hpl)
  
  X<-t(combn(Hnames,2))
  # remove combinations including at least one species with not enough collecting event => doesn't work on april 16
  # good_H <- H[which(H$n_col_event>=S_col_events),]
  # X <- X[is.element(X[,1],row.names(good_H)) & is.element(X[,2],row.names(good_H)),]
  
  nbtot<-nrow(X)
  X2<-cbind(X, 1:nrow(X))
  
  # had to be done once to compute plant triplets
  {
    ##########################
    ##WHAT FOLLOWS WAS DONE ONCE. NO NEED TO REDO IT (stored)
    ##########################
    # getLinkScoreAndNotLinkScore<-function(arr) {
    # 	setTxtProgressBar(pb, (as.numeric(arr[3])/nbtot))
    # 	Link<-sum(apply(hpl[c(arr[1],arr[2]),],2,min))
    # 	minEach<-min(apply(hpl[c(arr[1],arr[2]),],1,sum))
    # 	return(c(Link,minEach))
    # }
    # fun<-function() {
    # 	pb<-txtProgressBar(style=3)
    # 	res<-apply(X2,1,getLinkScoreAndNotLinkScore) ##VERY LONG! Do not regenerate ! 
    # 	close(pb)
    # 	return(res)
    # }
    # MinNumberSharedPlants<-fun()
    # ###WE SAVE ALL THE SESSION FOR LATER: 
    # #save(MinNumberSharedPlants, file="THE-LONG-HostPlantsData.Rdata")
    ##########################
    ##########################
    ##########################
  }
  
  load("../input/THE-LONG-HostPlantsData.Rdata")
  ##from this data it is easy to get those strong pairs and strong non pairs
  MinNumberSharedPlants<-t(MinNumberSharedPlants)
  MinNumberSharedPlants<-MinNumberSharedPlants[which(is.element(X[,1], host_tree$tip.label)&is.element(X[,2],host_tree$tip.label)),] ##KEEP those cases where both hosts are in the tree
  #idem for Hnames
  X<-X[which(is.element(X[,1], host_tree$tip.label)&is.element(X[,2],host_tree$tip.label)),]
  Hnames<-Hnames[is.element(Hnames, host_tree$tip.label)]
  
  X_lien<-X[which(MinNumberSharedPlants[,1]>=S_lien),]
  X_nl<-X[which((MinNumberSharedPlants[,1]==0)&(MinNumberSharedPlants[,2]>=S_nl)),]
  
  MAT<-matrix(0, nrow=length(Hnames), ncol=length(Hnames))
  colnames(MAT)<-Hnames
  rownames(MAT)<-Hnames
  MAT_nl<-MAT
  for (i in 1:nrow(X_nl)) {
    w<-X_nl[i,]
    MAT_nl[w[1],w[2]]<-MAT_nl[w[2],w[1]]<-1
  }
  ##For each link in X_lien (there are few), that we call H1-H3, 
  ## we search in MAT_nl all H2 so that 
  # A/ H1-H2 & H3-H2 are non links
  # B/ H3 is outside of the clade of H1 and H2
  # C/ H1 and H2 belong to the same clade
  # D/ d(H1,H2)>=D_seuil_minH1H2 
  H2forH1H3<-list()
  for (i in 1:nrow(X_lien)) {
    H3isH1validatedH2<-NULL
    H1isH1validatedH2<-NULL
    dh1h2<-NULL
    allnl<-MAT_nl[X_lien[i,],]
    possibleH2<-names(which(apply(allnl,2,sum)==2)) # condition A/
    grH2<-Clhost_tree[possibleH2]
    grH1<-Clhost_tree[X_lien[i,1]]
    grH3<-Clhost_tree[X_lien[i,2]]
    H3isH1validatedH2<-names(which((grH2==grH3)&(grH1!=grH2)))
    H1isH1validatedH2<-names(which((grH2==grH1)&(grH3!=grH2)))
    
    r1<-NULL
    r2<-NULL
    if (length(H3isH1validatedH2)>0) r1<-cbind(X_lien[i,2], H3isH1validatedH2, X_lien[i,1])
    if (length(H1isH1validatedH2)>0) r2<-cbind(X_lien[i,1], H1isH1validatedH2, X_lien[i,2])
    
    r1r2<-rbind(r1,r2)
    
    if (is.matrix(r1r2)==F){ # for cases where r1r2 has only one row and is not a matrix after last step
      r1r2 <- cbind(r1r2[1],r1r2[2],r1r2[3])
      dh1h2<-Dh[r1r2[,1],r1r2[,2]]
    }
    else{
      dh1h2<-diag(Dh[r1r2[,1],r1r2[,2]])
    }
    H2forH1H3[[i]]<-r1r2[(dh1h2>=D_seuil_minH1H2),] #condition D/ 
  }
  MAT<-do.call(rbind, H2forH1H3)
  return(MAT)

}


# a function to get all SGP duplets: pairs of parasites that are closely related and specialist vs generalist based on taxonomy (genus)
GetSGPDuplets<-function(threshold_n_obs=10, D_seuil_min_SG=0.1, low_par_spec_index = 0.34, high_par_spec_index = 0.98) {
  
  par_spec_tab <- P
  good_par_spec_tab <- par_spec_tab[which(par_spec_tab$n_obs_with_host>=threshold_n_obs),]
  
  tab_dist_P <- melt(Dp)

  gr1<-Clpar_tree[as.character(tab_dist_P$Var1)]
  gr2<-Clpar_tree[as.character(tab_dist_P$Var2)]

  
  good_tab_dist_P <- tab_dist_P[(tab_dist_P$value>=D_seuil_min_SG)&(gr1==gr2),]
  good_tab_dist_P$score1<-good_par_spec_tab[as.character(good_tab_dist_P$Var1),]$spec_index
  good_tab_dist_P$score2<-good_par_spec_tab[as.character(good_tab_dist_P$Var2),]$spec_index
    
  
  test_P_duplet <- (good_tab_dist_P$score1 <= low_par_spec_index & good_tab_dist_P$score2 >= high_par_spec_index & !is.na(good_tab_dist_P$score1+good_tab_dist_P$score2)) # +     (good_tab_dist_P$score1 >= high_par_spec_index & good_tab_dist_P$score2 <= low_par_spec_index) => je mets la fin en commentaires car toutes les paires d'espÃ¨ces sont listÃ© dans les 2 ordres...
  
  SGPDup_tab <- good_tab_dist_P[test_P_duplet,]
  SGPDup_tab$group<-Clpar_tree[as.character(SGPDup_tab$Var1)]
  
  SGPDup_mat <- array(0, c(length(row.names(P)), length(row.names(P))))
  dimnames(SGPDup_mat)<-list(row.names(P),row.names(P))
  
  for (i in 1:nrow(SGPDup_tab)) {
    SGPDup_mat[SGPDup_tab[i,1],SGPDup_tab[i,2]]<-SGPDup_tab$group[i]
    SGPDup_mat[SGPDup_tab[i,2],SGPDup_tab[i,1]]<-SGPDup_tab$group[i]
    # SGPDup_mat[SGPDup_tab[i,1],SGPDup_tab[i,2]]<-1
    # SGPDup_mat[SGPDup_tab[i,2],SGPDup_tab[i,1]]<-1
  }
  return(SGPDup_mat)
}
# geometric mean
geomean<-function(x,eps=1e-6) {
	exp(mean(log(x + eps)))
}
# Get list of species from a list of triplets in string format (with & separators)
GetAllMono<-function(triplets) {
	return(unique(unlist(strsplit(as.character(triplets),"&"))))
}

# a function to get triplet scores from a list of triplets strings + one new potential triplet string
GetTripletScores<-function(newtriplet, selectedtriplets) {
	triplets<-c(newtriplet, selectedtriplets)
	allmono <-  GetAllMono(triplets)
	
	# separer H et P
	all_mono_H <- unique(row.names(H)[which(is.element(row.names(H),allmono))])
	all_mono_P <- unique(row.names(P)[which(is.element(row.names(P),allmono))])

  # the -1 in the score is because we don't want to count the value of 0 in the cubes as a cluster.
	sc_HHP<-length(unique(array(CUBE_HHP[all_mono_H,all_mono_H,all_mono_P])))-1
	sc_PPH<-length(unique(array(CUBE_PPH[all_mono_P,all_mono_P,all_mono_H])))-1
	sc_HHH<-length(unique(array(CUBE_HHH[all_mono_H,all_mono_H,all_mono_H])))-1
  # we add the duplet score:
  sc_SGP<-length(unique(array(SGPDup_mat[all_mono_P,all_mono_P])))-1
	
	return(c(sc_HHP, sc_PPH, sc_HHH, sc_SGP))
}

# a function to get triplet scores from a list of triplets strings + one new potential triplet string
GetTripletScores2<-function(availabletriplets, selectedtriplets) {
	allmono1 <-  GetAllMono(selectedtriplets)
	all_mono1_H <- unique(row.names(H)[which(is.element(row.names(H),allmono1))])
	all_mono1_P <- unique(row.names(P)[which(is.element(row.names(P),allmono1))])

	allmono2 <-  str2trip(availabletriplets)
	splitHandP<-function(arr, h_s, p_s) {
		is_H <- unique(row.names(H)[which(is.element(row.names(H),arr))])
		is_P <- unique(row.names(P)[which(is.element(row.names(P),arr))])
		is_H<-unique(c(is_H,h_s))
		is_P<-unique(c(is_P,p_s))
		res<-NULL
		res$is_H<-is_H
		res$is_P<-is_P
		res
	}
	splitted<-apply(allmono2,1,splitHandP, h_s=all_mono1_H, p_s=all_mono1_P)
	getscoresFromListOfList<-function(L) {
		sc_HHP<-length(unique(array(CUBE_HHP[L$is_H,L$is_H,L$is_P])))-1
		sc_PPH<-length(unique(array(CUBE_PPH[L$is_P,L$is_P,L$is_H])))-1
		sc_HHH<-length(unique(array(CUBE_HHH[L$is_H,L$is_H,L$is_H])))-1
		# we add the duplet score:
		sc_SGP<-length(unique(array(SGPDup_mat[L$is_P,L$is_P])))-1	
		return(c(sc_HHP, sc_PPH, sc_HHH, sc_SGP))
	}
	AllScores<-lapply(splitted, getscoresFromListOfList)
	RES<-do.call(cbind, AllScores)
	colnames(RES)<-availabletriplets
	RES
}




# a function to get triplet scores when removing teh first element in the second list. 
GetTripletScores_rm<-function(triplettoremove, selectedtriplets) {
	#for each triplet removed, the species composing it are remioved from all the others.
	allmono<-GetAllMono(selectedtriplets)
	allmono<-setdiff(allmono, array(str2trip(triplettoremove)))

	# separer H et P
	all_mono_H <- unique(row.names(H)[which(is.element(row.names(H),allmono))])
	all_mono_P <- unique(row.names(P)[which(is.element(row.names(P),allmono))])

  # the -1 in the score is because we don't want to count the value of 0 in the cubes as a cluster.
	sc_HHP<-length(unique(array(CUBE_HHP[all_mono_H,all_mono_H,all_mono_P])))-1
	sc_PPH<-length(unique(array(CUBE_PPH[all_mono_P,all_mono_P,all_mono_H])))-1
	sc_HHH<-length(unique(array(CUBE_HHH[all_mono_H,all_mono_H,all_mono_H])))-1
  # we add the duplet score:
    sc_SGP<-length(unique(array(SGPDup_mat[all_mono_P,all_mono_P])))-1
	
	return(c(sc_HHP, sc_PPH, sc_HHH, sc_SGP))
}
# a function to put a triplet list (tab format) into string format (not removing duplicates, if there are any...)
trip2str<-function(alltriplets){
	alltrip_string<-apply(alltriplets[,1:3],1,function(x) paste(sort(x), collapse="&"))
	return(alltrip_string)
}
# the opposite function
str2trip<-function(alltriplets) {
	return(do.call(rbind, strsplit(alltriplets,"&")))
}

#function to convert a list of species to a list of selected triplets, and a list of availbalble triplets, given a list of total triplets.
Sp2Triplets<-function(species, allTRIP) {
	allTRIPstr<-trip2str(allTRIP)
	trip_tab_na<-matrix(match(allTRIP, species), ncol=3,byrow=F)
	is_selected<-!is.na(apply(trip_tab_na,1,sum))
	selected_triplets_real<-allTRIPstr[is_selected]
	available_triplets_real<-allTRIPstr[!is_selected]
	return(list(selected=selected_triplets_real, available=available_triplets_real))
}
printit<-function(scs, scg, sps) {
	cat("scores: ",paste(scs, collapse="-"),"\n")
	cat("geomean score: ",scg,"\n")
	cat("species: ",sps,"\n")
}


TRP<-GetTriplets(S_lien=20, S_nl=40, D_seuil_min=0.1) ##HHP and PPH triplets; uncomment if recomputing
TRPl<-GetTripletPlants(S_lien=20,S_nl=40, D_seuil_minH1H2=0.1) #HHH triplets (plant triplets); uncomment if recomputing
SGPDup_mat<-GetSGPDuplets(threshold_n_obs=5, D_seuil_min_SG=0.1, low_par_spec_index = 0.34, high_par_spec_index = 0.98)


## Put triplets in cube format. NEW: we put in the cube the 
## NEW: we put in the cube the number of the cluster the H1/H2 or P1/P2 belong to.
All_H<-unique(c(TRP$TRIPP[,3], TRP$TRIPH[,1], TRP$TRIPH[,2], array(TRPl)))
All_P<-unique(c(TRP$TRIPH[,3], TRP$TRIPP[,1], TRP$TRIPP[,2]))
number_of_H<-length(All_H)
number_of_P<-length(All_P)

CUBE_HHP<- array(0, c(number_of_H,number_of_H,number_of_P))
dimnames(CUBE_HHP)<-list(All_H,All_H, All_P)
for (i in 1:nrow(TRP$TRIPH)) {
  CUBE_HHP[TRP$TRIPH[i,1],TRP$TRIPH[i,2],TRP$TRIPH[i,3]]<-Clhost_tree[TRP$TRIPH[i,1]]
}

CUBE_PPH<- array(0, c(number_of_P,number_of_P,number_of_H))
dimnames(CUBE_PPH)<-list(All_P,All_P, All_H)
for (i in 1:nrow(TRP$TRIPP)) {
  CUBE_PPH[TRP$TRIPP[i,1],TRP$TRIPP[i,2],TRP$TRIPP[i,3]]<-Clpar_tree[TRP$TRIPP[i,1]]
}

CUBE_HHH<- array(0, c(number_of_H,number_of_H,number_of_H))
dimnames(CUBE_HHH)<-list(All_H,All_H, All_H)
for (i in 1:nrow(TRPl)) {
  CUBE_HHH[TRPl[i,1],TRPl[i,2],TRPl[i,3]]<-Clhost_tree[TRPl[i,1]]
}
  



trip_HHP <- TRP$TRIPH
trip_PPH <- TRP$TRIPP
trip_HHH <- TRPl
ALLTRIPLETS_TOTAL <- rbind(trip_HHP,trip_PPH,trip_HHH)
AllTriplets<-trip2str(ALLTRIPLETS_TOTAL)

########## PARAMETERS FOR OPTIMIZATION
# Number of triplets per category for start.
Freq.remove<-10 #every 10rounds, instead of adding best triplet, we remove the one with less impact when removed.
Freq.reload_all_TRIPLETS<-30
cpt<-0 #count the number of loops of the optimization process
maxSP<-125 #max number of species to select (stoping criteria)
n_rand_trip_HHP<-5 # nb of randomly chosen HHP triplets for starting
n_rand_trip_PPH<-5 # nb of randomly chosen PPH triplets for starting
n_rand_trip_HHH<-5 # nb of randomly chosen HHH triplets for starting
SampleSizeAvailable<-40
SAVESCORE<-array() #TO SAVE SCORES.

######### INITIATION: 
##### the code is run 100 times, picking the nb of random triplets of each type
##### and the best selection of triplets encounterd (highets score) is kept. 
#####
startscore<-0
for (i in 1:100) {
	if (i==1) cat(paste("init ","#"," - ",sep=" "), paste("number of species","scores (hhp-pph-hhh-dupl)", "geomean_score", sep="    "), "\n")
	# Initiation
	rand_HHP_trip <-trip_HHP[sample(nrow(trip_HHP), n_rand_trip_HHP), ] # pick randomly some HHP triplets
	rand_PPH_trip <-trip_PPH[sample(nrow(trip_PPH), n_rand_trip_PPH), ] # pick randomly some PPH triplets
	rand_HHH_trip <-trip_HHH[sample(nrow(trip_HHH), n_rand_trip_HHH), ] # pick randomly some HHH triplets
	SpeciesSelected<-unique(array(rbind(rand_HHP_trip,rand_PPH_trip, rand_HHH_trip)))
	UPDATED_TRIPLETS<-Sp2Triplets(SpeciesSelected, ALLTRIPLETS_TOTAL)
	SelectedTrip<-UPDATED_TRIPLETS$selected ##REAL selected
	AvailableTrip<-UPDATED_TRIPLETS$available ##REAL available
	n_rand_trip <- length(SelectedTrip)
	# calculate the score from the current list of selected
	TripScore_s_Init <- GetTripletScores(newtriplet=SelectedTrip[n_rand_trip], selectedtriplets=SelectedTrip[-n_rand_trip])
	# Initial score:
	TripScoreInit <- geomean(TripScore_s_Init)
	if (TripScoreInit>startscore) {
		startscore<-TripScoreInit
		SelectedTripOK<-SelectedTrip
		AvailableTripOK<-AvailableTrip
	}
	cat(paste("init ",i," - ",sep=" "), paste(length(SpeciesSelected),paste(TripScore_s_Init, collapse="-"), geomean(TripScore_s_Init), sep="    "), "\n")
}

########### STARTING CONDITIONS
SelectedTrip<-SelectedTripOK
AvailableTrip<-AvailableTripOK
n_rand_trip <- length(SelectedTrip)
TripScore_s_Init <- GetTripletScores(newtriplet=SelectedTrip[n_rand_trip], selectedtriplets=SelectedTrip[-n_rand_trip])
TripScoreInit <- geomean(TripScore_s_Init)
printit(TripScore_s_Init,TripScoreInit, length(GetAllMono(SelectedTrip)))

########### OPTIMIZATION.
while (length(GetAllMono(SelectedTrip))<maxSP) {
	cpt <- cpt + 1
	cat(paste("cpt: ",cpt,"\n"))
	cat(paste(length(GetAllMono(SelectedTrip))," species -> ",sep=""))
	### CORRIGER ça
	if (cpt%%Freq.remove==0) {
		GetScoresDetails<-sapply(SelectedTrip, GetTripletScores_rm, selectedtriplets=SelectedTrip)
		GetScores<-apply(GetScoresDetails,2,geomean)
		TripletOfInterest<-sample(names(which(GetScores==max(GetScores))),1)
		# get list of species after removing that triplet
		NewListOfSpecies<-setdiff(GetAllMono(SelectedTrip), array(str2trip(TripletOfInterest)))
		cat(paste("removed",TripletOfInterest,"\n"))
	}

	else {
		# GetScoresDetails<-sapply(AvailableTrip, GetTripletScores, selectedtriplets=SelectedTrip)
		GetScoresDetails<-GetTripletScores2(AvailableTrip, SelectedTrip)
		GetScores<-apply(GetScoresDetails,2,geomean)
		TripletOfInterest<-sample(names(which(GetScores==max(GetScores))),1)
		# get new list of species after adding those in the new triplet.
		NewListOfSpecies <- unique(c(GetAllMono(SelectedTrip), array(str2trip(TripletOfInterest))))
		cat(paste("added",TripletOfInterest,"\n"))
	}
	if ((cpt == 1) || ((cpt-1)%%Freq.reload_all_TRIPLETS==0)) { ##we just did the first run OR we are just after a turn where we retook everyone 
		cat("\n### RESAMPLING BEST TRIPLETS ###\n")
		# we order all Triplets by there scores for each category of score
		OrderedMatricScores<-apply(GetScoresDetails,1,order)
		# we count how many different triplets are selected given where we cut in the list of ordered triplets
		NbTripletsConsidered<-sapply(1:nrow(OrderedMatricScores),function(x, D) length(unique(array(D[x:nrow(D),]))), D=OrderedMatricScores)
		# we only consider best triplets so that the number of triplets considered equals SampleSizeAvailable
		FromeWhereToGet<-which(NbTripletsConsidered<=SampleSizeAvailable)[1]
		# get all triplets from this index:
		TripletsToKeep<-unique(array(OrderedMatricScores[FromeWhereToGet:nrow(OrderedMatricScores),]))
		TripletsToKeepName<-AvailableTrip[TripletsToKeep]
		##we create the list of 'ALLTRIPLETS' that we will keep for some Time.
		ALLSPECIES<-GetAllMono(c(SelectedTrip,TripletsToKeepName))
		ALLTRIPLETS<-str2trip(Sp2Triplets(ALLSPECIES, ALLTRIPLETS_TOTAL)$selected)
	}
	if (cpt%%Freq.reload_all_TRIPLETS==0) {
		cat("\n### REINJECTING ALL TRIPLETS... ###\n")
		ALLTRIPLETS<-ALLTRIPLETS_TOTAL
	}

	UPDATED_TRIPLETS <- Sp2Triplets(NewListOfSpecies, ALLTRIPLETS)
	SelectedTrip<-UPDATED_TRIPLETS$selected ##REAL selected
	AvailableTrip<-UPDATED_TRIPLETS$available ##REAL available
	# print results

	cat("---\nplaying with ",nrow(ALLTRIPLETS), " triplets in total\n")
	printit(GetScoresDetails[,TripletOfInterest], geomean(GetScoresDetails[,TripletOfInterest]), length(GetAllMono(SelectedTrip)))
	cat("--\n")
#	cat(paste(c("sc_HHP:", "sc_PPH:", "sc_HHH:", "sc_SGP:"),GetScoresDetails[,TripletOfInterest]))
#	cat("\n")

}


