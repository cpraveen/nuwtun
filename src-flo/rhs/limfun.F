C
C     This file is part of NuWTun, see <http://nuwtun.berlios.de>, and was
C     originally taken from ISAAC Version 4.2, release date October 2001. 
C     This file may have been modified; for a list of changes, see the 
C     changes.txt file in the docs directory and the subversion log.
C
C     Portions Copyright (C) 2001 Joseph H. Morrison
C
C     This code is part of ISAAC.
C
C     This program is distributed under the terms of the ISAAC Public Source
C     License. This program is distributed WITHOUT ANY WARRANTY; without
C     even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C     PURPOSE. 
C
C     You should have received a copy of the ISAAC Public Source License
C     with this program. If you did not, you may get a copy of the license
C     at <http://isaac-cfd.sourceforge.net>
C
      SUBROUTINE SMILIM(DQP1, DQP2, DQM1)
      SMI(D1, D2, EPSQ) = ( D1*(D2*D2+2.E0*EPSQ) 
     1                    + D2*(2.E0*D1*D1+EPSQ) )
     2                / ( 2.E0*D1*D1 - D1*D2 + 2.E0*D2*D2 + 3.E0*EPSQ )
      EPSQ = 1.0E-4
      XP   = SMI( DQP1, DQP2, EPSQ )
      XM   = SMI( DQP2, DQP1, EPSQ ) 
      DQP1 = - 0.5E0 * XP
      DQM1 =   0.5E0 * XM
      RETURN
      END
C
C 2.  MinMod Limiter
C
      SUBROUTINE RMMLIM(DQP1, DQP2, DQM1, OMKAP, OPKAP, B)
      RMINMD(D1, D2, B) = 0.5E0 * ( SIGN(1.E0,D1) + SIGN(1.E0,D2) )
     1                          * MIN ( ABS(D1), ABS(B*D2) )
      XM   = RMINMD( DQP1, DQP2, B )
      XP   = RMINMD( DQP2, DQP1, B ) 
      DQP1 = - 0.25E0 * ( OMKAP*XP + OPKAP*XM )
      DQM1 =   0.25E0 * ( OMKAP*XM + OPKAP*XP )
      RETURN
      END
C
C 3.  Smooth Limiter
C
      SUBROUTINE SMTLIM(DQP1, DQP2, DQM1, RKAP)
      SMTH(D1, D2, EP) = ( 2.E0*D1*D2 + EP ) 
     1                 / ( D1*D1 + D2*D2 + EP )
      EPS   = 1.E-06
      XP    = DQP2
      XM    = DQP1
      SM    = SMTH (XM, XP, EPS)
      DQP1  = - 0.25E0 * SM * 
     1                ( (1.E0 - SM*RKAP)*XP
     2              + (1.E0 + SM*RKAP)*XM )
      DQM1  =   0.25E0 * SM * 
     1                ( (1.E0 - SM*RKAP)*XM
     2              + (1.E0 + SM*RKAP)*XP )
      RETURN
      END
C
C 4.  Unlimited Scheme
C
      SUBROUTINE RNOLIM(DQP1, DQP2, DQM1, OMKAP, OPKAP)
      XP   = DQP2
      XM   = DQP1
      DQP1 = - 0.25E0*(OMKAP*XP + OPKAP*XM)
      DQM1 =   0.25E0*(OMKAP*XM + OPKAP*XP)
      RETURN
      END
