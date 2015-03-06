SUBROUTINE mean_shift(numproc,nblimit,nbideal,dataw,bandWidth)
	USE module_structure

   !INCLUDE 'mpif.h'
    !IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
  
    INTEGER :: nbideal
    INTEGER :: nblimit
    INTEGER :: numproc
    DOUBLE PRECISION :: bandWidth !bandwidth parameter

    !=== IN/OUT ===
    TYPE(type_data) :: dataw
    
    !###########################################
    ! DECLARATIONS
    !###########################################      
    INTEGER ::point_num												!number of points
    INTEGER ::dim_num													!number of dimensions
    INTEGER ::cluster_num												!number of clusters
    
    !#### Variables  ####
    INTEGER :: numClust 												!the cluster number  
    DOUBLE PRECISION :: bandSq										!square of bandWidth
    DOUBLE PRECISION :: stopThresh									!when mean has converged
    INTEGER :: beenVisitedFlag(dataw%nb)								!track if a point has been seen already
    INTEGER :: numInitPts												!number of points to possibly use as initialization points
    INTEGER :: thisClusterVotes(dataw%nb)								!used to resolve conflicts on cluster membership
    INTEGER :: stInd													!start point of mean
    DOUBLE PRECISION :: myMean(dataw%dim)								!mean of this cluster
    DOUBLE PRECISION :: myOldMean(dataw%dim)							!old mean computed for this cluster
    INTEGER :: myMembers(dataw%nb)										!1 if the point belongs to the cluster, else 0
    INTEGER :: mergeWith												!used to merge clusters
    DOUBLE PRECISION :: clustCent(dataw%dim,dataw%nbclusters)			!centers of each cluster
	INTEGER :: clusterVotes(dataw%nbclusters,dataw%nb)					!number of votes for each point for each cluster
    INTEGER :: i
    INTEGER :: j
    INTEGER :: num
    DOUBLE PRECISION :: sqDist
    
    
    
    !INITIALIZE STUFF    
    point_num=dataw%nb
    dim_num=dataw%dim
    numClust = 1
    bandSq = bandWidth**2
    stopThresh = 1e-3*bandWidth
    beenVisitedFlag(:) = 0
    numInitPts = point_num
    clusterVotes(:,:) = 0		
    
    DO WHILE (numInitPts>0)
    
		!take the first point as start of mean
		DO i=1, point_num
			IF (beenVisitedFlag(i)==0) THEN
				stInd = i
				EXIT
			ENDIF
		ENDDO
		myMean = dataw%point(stInd)%coord	!initialize mean to this points location
		DO j=1, dim_num
			myMean(j) = dataw%point(i)%coord(j)
		ENDDO
		myMembers(:) = 0
		thisClusterVotes(:) = 0	!used to resolve conflicts on cluster membership
		
		DO
		
			DO i=1, point_num
				!dist squared from mean to all points still active
				sqDist = 0
				DO j=1, dim_num
					sqDist = sqDist + (dataw%point(i)%coord(j) - myMean(j))**2
				ENDDO
				IF (sqDist < bandSq) THEN
					thisClusterVotes(i) = thisClusterVotes(i) + 1	!add a vote for all the in points belonging to this cluster
					myMembers(i) = 1								!add any point within bandWidth to the cluster
					beenVisitedFlag(i) = 1							!mark that these points have been visited
				ENDIF
			ENDDO
			
			myOldMean = myMean
			
			!compute the new mean
			DO i=1, point_num
				num = 0
				IF (myMembers(i)==1) THEN
					DO j=1, dim_num
						myMean(j) = myMean(j) + dataw%point(i)%coord(j)
					ENDDO
					num = num + 1
				ENDIF
			ENDDO
			myMean = myMean/num
			
			!compute the distance from myMean to myOldMean
			sqDist = 0
			DO j=1, dim_num
				sqDist = sqDist + (myOldMean(j) - myMean(j))**2
			ENDDO
			
			!if mean doesn't move much stop this cluster
			IF (sqDist < stopThresh**2) THEN
			
				!check for merge posibilities
				mergeWith = 0
				DO cN=1, numclust-1
					!compute the distance from possible new clust max to old clust max
					sqDist = 0
					DO j=1, dim_num
						sqDist = sqDist + (clustCent(j,cN) - myMean(j))**2
					ENDDO
					IF (sqDist < (bandWidth/2)**2) THEN
						mergeWith = cN
						EXIT
					ENDIF
				ENDDO
				
				IF (mergeWith > 0) THEN		!something to merge
				
					clustCent(:,mergeWith) = (myMean+clustCent(:,mergeWith))/2					!mean of centers
					clusterVotes(mergeWith,:) = clusterVotes(mergeWith,:) + thisClusterVotes;	!add these votes to the merged cluster
					
				ELSE
				
					numClust = numClust + 1
					clustCent(:,numClust) = myMean
					clusterVotes(numClust,:) = thisClusterVotes
				ENDIF
				EXIT
				
			ENDIF
		
		ENDDO
			
		DO i=1, point_num
			numInitPts = 0
			IF (beenVisitedFlag(i)==0) THEN
				numInitPts = numInitPts + 1
			ENDIF
		ENDDO
			
	ENDDO
	
	DO i=1, point_num
		dataw%point(i)%cluster = MAXLOC(clusterVotes(:,i), DIM=1)
	ENDDO
	
END SUBROUTINE mean_shift
						
			
			
			
			
			
