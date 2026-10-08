#if !defined(LFRIC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

PROGRAM ukca_main

USE ukca_shell_mod, ONLY: ukca_shell

IMPLICIT NONE

! A simple top level routine for the UKCA box model.

CALL ukca_shell()

END PROGRAM ukca_main
#endif
