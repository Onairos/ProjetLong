PROGRAM teste_clusters
  USE module_teste_clusters

  IMPLICIT NONE
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  TYPE(type_test),DIMENSION(:),POINTER :: test
  CHARACTER*30 :: files
  INTEGER :: i
  INTEGER :: nbtests
  
  !###########################################
  ! INSTRUCTIONS
  !###########################################
  PRINT *
  PRINT *,'-------------------------------------------------------'
  PRINT *,'Programme de test de CLUSTERS et VISUCLUSTERS '
  PRINT *,'-------------------------------------------------------'
  PRINT *
  PRINT *,'teste pour plusieurs configurations differentes :'
  PRINT *,'  - le programme CLUSTERS'
  PRINT *,'  - le programme VISUCLUSTERS pour des sorties paraview'
  PRINT *,'  - le programme VISUCLUSTERS pour des sorties gmsh'
  PRINT *
  PRINT *,'le programme verifie l''existence d''un fichier genere'
  PRINT *,'  en fin d''execution de ces programmes.'
  PRINT *,'il retourne "T" ou "F" selon que ce programme s''est'
  PRINT *,'  bien execute ou non.'
  PRINT *
  PRINT *,'-------------------------------------------------------'
  
  !creation des datas
  PRINT *
  PRINT *,'creation des datas...'
  CALL create_data

  !declaration des tests
  nbtests=12
  ALLOCATE(test(nbtests))

  !test coord mono
  test(1)%dir='test_coord_mono'
  test(1)%nbproc=1
  test(1)%output='verif.coord_mono.txt'
  test(1)%visup='paraview.coord_mono.txt'
  test(1)%visug='gmsh.coord_mono.txt'
  test(1)%fichier='cible'
  test(1)%datatype='COORD'
  test(1)%decoupetype='INTERFACE'
  ALLOCATE(test(1)%decoupe(2))
  test(1)%decoupe(:)=1
  test(1)%epaisseur=0.

  !test coord multi interface
  test(2)%dir='test_coord_multi_interface'
  test(2)%nbproc=5
  test(2)%output='verif.coord_multi_interface.txt'
  test(2)%visup='paraview.coord_multi_interface.txt'
  test(2)%visug='gmsh.coord_multi_interface.txt'
  test(2)%fichier='cible'
  test(2)%datatype='COORD'
  test(2)%decoupetype='INTERFACE'
  ALLOCATE(test(2)%decoupe(2))
  test(2)%decoupe(1)=2
  test(2)%decoupe(2)=2
  test(2)%epaisseur=0.3

  !test coord multi recouvrement
  test(3)%dir='test_coord_multi_recouvre'
  test(3)%nbproc=4
  test(3)%output='verif.coord_multi_recouvre.txt'
  test(3)%visup='paraview.coord_multi_recouvre.txt'
  test(3)%visug='gmsh.coord_multi_recouvre.txt'
  test(3)%fichier='cible'
  test(3)%datatype='COORD'
  test(3)%decoupetype='RECOUVREMENT'
  ALLOCATE(test(3)%decoupe(2))
  test(3)%decoupe(1)=2
  test(3)%decoupe(2)=2
  test(3)%epaisseur=0.3

  !test image mono
  test(4)%dir='test_image_mono'
  test(4)%nbproc=1
  test(4)%output='verif.image_mono.txt'
  test(4)%visup='paraview.image_mono.txt'
  test(4)%visug='gmsh.image_mono.txt'
  test(4)%fichier='image1d'
  test(4)%datatype='IMAGE'
  test(4)%decoupetype='INTERFACE'
  ALLOCATE(test(4)%decoupe(2))
  test(4)%decoupe(:)=1
  test(4)%epaisseur=0.

  !test image multi interface
  test(5)%dir='test_image_multi_interface'
  test(5)%nbproc=7
  test(5)%output='verif.image_multi_interface.txt'
  test(5)%visup='paraview.image_multi_interface.txt'
  test(5)%visug='gmsh.image_multi_interface.txt'
  test(5)%fichier='image1d'
  test(5)%datatype='IMAGE'
  test(5)%decoupetype='INTERFACE'
  ALLOCATE(test(5)%decoupe(2))
  test(5)%decoupe(1)=3
  test(5)%decoupe(2)=2
  test(5)%epaisseur=1.01

  !test image multi recouvrement
  test(6)%dir='test_image_multi_recouvre'
  test(6)%nbproc=6
  test(6)%output='verif.image_multi_recouvre.txt'
  test(6)%visup='paraview.image_multi_recouvre.txt'
  test(6)%visug='gmsh.image_multi_recouvre.txt'
  test(6)%fichier='image1d'
  test(6)%datatype='IMAGE'
  test(6)%decoupetype='RECOUVREMENT'
  ALLOCATE(test(6)%decoupe(2))
  test(6)%decoupe(1)=3
  test(6)%decoupe(2)=2
  test(6)%epaisseur=1.01

  !test seuil mono
  test(7)%dir='test_seuil_mono'
  test(7)%nbproc=1
  test(7)%output='verif.seuil_mono.txt'
  test(7)%visup='paraview.seuil_mono.txt'
  test(7)%visug='gmsh.seuil_mono.txt'
  test(7)%fichier='image1d'
  test(7)%datatype='SEUIL'
  test(7)%decoupetype='INTERFACE'
  ALLOCATE(test(7)%decoupe(1))
  test(7)%decoupe(:)=1
  test(7)%epaisseur=0.

  !test seuil multi interface
  test(8)%dir='test_seuil_multi_interface'
  test(8)%nbproc=9
  test(8)%output='verif.seuil_multi_interface.txt'
  test(8)%visup='paraview.seuil_multi_interface.txt'
  test(8)%visug='gmsh.seuil_multi_interface.txt'
  test(8)%fichier='image1d'
  test(8)%datatype='SEUIL'
  test(8)%decoupetype='INTERFACE'
  ALLOCATE(test(8)%decoupe(1))
  test(8)%decoupe(:)=8
  test(8)%epaisseur=0.01

  !test seuil multi recouvrement
  test(9)%dir='test_seuil_multi_recouvre'
  test(9)%nbproc=8
  test(9)%output='verif.seuil_multi_recouvre.txt'
  test(9)%visup='paraview.seuil_multi_recouvre.txt'
  test(9)%visug='gmsh.seuil_multi_recouvre.txt'
  test(9)%fichier='image1d'
  test(9)%datatype='SEUIL'
  test(9)%decoupetype='RECOUVREMENT'
  ALLOCATE(test(9)%decoupe(1))
  test(9)%decoupe(:)=8
  test(9)%epaisseur=0.01

  !test geom mono
  test(10)%dir='test_geom_mono'
  test(10)%nbproc=1
  test(10)%output='verif.geom_mono.txt'
  test(10)%visup='paraview.geom_mono.txt'
  test(10)%visug='gmsh.geom_mono.txt'
  test(10)%fichier='image1d'
  test(10)%datatype='GEOM'
  test(10)%decoupetype='INTERFACE'
  ALLOCATE(test(10)%decoupe(3))
  test(10)%decoupe(:)=1
  test(10)%epaisseur=0.

  !test geom multi interface
  test(11)%dir='test_geom_multi_interface'
  test(11)%nbproc=28
  test(11)%output='verif.geom_multi_interface.txt'
  test(11)%visup='paraview.geom_multi_interface.txt'
  test(11)%visug='gmsh.geom_multi_interface.txt'
  test(11)%fichier='image1d'
  test(11)%datatype='GEOM'
  test(11)%decoupetype='INTERFACE'
  ALLOCATE(test(11)%decoupe(3))
  test(11)%decoupe(:)=3
  test(11)%epaisseur=0.01

  !test geom multi recouvrement
  test(12)%dir='test_geom_multi_recouvre'
  test(12)%nbproc=27
  test(12)%output='verif.geom_multi_recouvre.txt'
  test(12)%visup='paraview.geom_multi_recouvre.txt'
  test(12)%visug='gmsh.geom_multi_recouvre.txt'
  test(12)%fichier='image1d'
  test(12)%datatype='GEOM'
  test(12)%decoupetype='RECOUVREMENT'
  ALLOCATE(test(12)%decoupe(3))
  test(12)%decoupe(:)=3
  test(12)%epaisseur=0.01

  !lancement
  files=''
  DO i=1,nbtests
     PRINT *
     PRINT *,'TEST : '//test(i)%dir
     IF (files/='t') THEN
        PRINT *,'  > lancer le test ? [o,n,t]'
        READ *,files
     ENDIF
     IF (files/='n') THEN
        CALL create_executable(test(i))
        CALL create_test(test(i))  
        CALL execute_test(test(i))
     ENDIF
  ENDDO
  PRINT *
  PRINT *,'-------------------------------------------------------'

END PROGRAM teste_clusters
