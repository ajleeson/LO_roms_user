      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: fennel_inp.h 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in Fennel et al. (2006) ecosystem model input    !
!  parameters. They are specified in input script "fennel.in".         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status
      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc
      real(r8), dimension(NBT,Ngrids) :: Rbio
      real(dp), dimension(nRval) :: Rval
      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in Fennel et al. (2006) biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lbiology')
              Npts=load_l(Nval, Cval, Ngrids, Lbiology)
            CASE ('BioIter')
              Npts=load_i(Nval, Rval, Ngrids, BioIter)
            CASE ('AttSW')
              Npts=load_r(Nval, Rval, Ngrids, AttSW)
            CASE ('AttChl')
              Npts=load_r(Nval, Rval, Ngrids, AttChl)
            CASE ('PARfrac')
              Npts=load_r(Nval, Rval, Ngrids, PARfrac)
            CASE ('Vp0')
              Npts=load_r(Nval, Rval, Ngrids, Vp0)
            CASE ('I_thNH4')
              Npts=load_r(Nval, Rval, Ngrids, I_thNH4)
            CASE ('D_p5NH4')
              Npts=load_r(Nval, Rval, Ngrids, D_p5NH4)
            CASE ('NitriR')
              Npts=load_r(Nval, Rval, Ngrids, NitriR)
            CASE ('K_NO3')
              Npts=load_r(Nval, Rval, Ngrids, K_NO3)
            CASE ('K_NH4')
              Npts=load_r(Nval, Rval, Ngrids, K_NH4)
            CASE ('K_PO4')
              Npts=load_r(Nval, Rval, Ngrids, K_PO4)
            CASE ('K_Phy')
              Npts=load_r(Nval, Rval, Ngrids, K_Phy)
            CASE ('Chl2C_m')
              Npts=load_r(Nval, Rval, Ngrids, Chl2C_m)
            CASE ('ChlMin')
              Npts=load_r(Nval, Rval, Ngrids, ChlMin)
            CASE ('PhyCN')
              Npts=load_r(Nval, Rval, Ngrids, PhyCN)
            CASE ('R_P2N')
              Npts=load_r(Nval, Rval, Ngrids, R_P2N)
            CASE ('PhyIP')
              Npts=load_r(Nval, Rval, Ngrids, PhyIP)
            CASE ('PhyIS')
              Npts=load_r(Nval, Rval, Ngrids, PhyIS)
            CASE ('PhyMin')
              Npts=load_r(Nval, Rval, Ngrids, PhyMin)
            CASE ('PhyMR')
              Npts=load_r(Nval, Rval, Ngrids, PhyMR)
            CASE ('ZooAE_N')
              Npts=load_r(Nval, Rval, Ngrids, ZooAE_N)
            CASE ('ZooBM')
              Npts=load_r(Nval, Rval, Ngrids, ZooBM)
            CASE ('ZooCN')
              Npts=load_r(Nval, Rval, Ngrids, ZooCN)
            CASE ('ZooER')
              Npts=load_r(Nval, Rval, Ngrids, ZooER)
            CASE ('ZooGR')
              Npts=load_r(Nval, Rval, Ngrids, ZooGR)
            CASE ('ZooMin')
              Npts=load_r(Nval, Rval, Ngrids, ZooMin)
            CASE ('ZooMR')
              Npts=load_r(Nval, Rval, Ngrids, ZooMR)
            CASE ('LDeRRN')
              Npts=load_r(Nval, Rval, Ngrids, LDeRRN)
            CASE ('LDeRRC')
              Npts=load_r(Nval, Rval, Ngrids, LDeRRC)
            CASE ('CoagR')
              Npts=load_r(Nval, Rval, Ngrids, CoagR)
            CASE ('SDeRRN')
              Npts=load_r(Nval, Rval, Ngrids, SDeRRN)
            CASE ('SDeRRC')
              Npts=load_r(Nval, Rval, Ngrids, SDeRRC)
            CASE ('RDeRRN')
              Npts=load_r(Nval, Rval, Ngrids, RDeRRN)
            CASE ('RDeRRC')
              Npts=load_r(Nval, Rval, Ngrids, RDeRRC)
            CASE ('wPhy')
              Npts=load_r(Nval, Rval, Ngrids, wPhy)
            CASE ('wLDet')
              Npts=load_r(Nval, Rval, Ngrids, wLDet)
            CASE ('wSDet')
              Npts=load_r(Nval, Rval, Ngrids, wSDet)
            CASE ('pCO2air')
              Npts=load_r(Nval, Rval, Ngrids, pCO2air)
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  nl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  nl_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu2(i,ng)=Rbio(itrc,ng)
                  tl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu4(i,ng)=Rbio(itrc,ng)
                  tl_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('LtracerSponge')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSponge(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Akt_bak(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_AKT_fac')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_Akt_fac(i,ng)=Rbio(itrc,ng)
                  tl_Akt_fac(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, NBT, Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Tnudg(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('Hadvection')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idbio(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idbio(iTrcStr), idbio(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Hadvection)
            CASE ('Vadvection')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idbio(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idbio(iTrcStr), idbio(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Vadvection)
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idbio(iTrcStr), idbio(iTrcEnd),             &
     &                      Vname(1,idTvar(idbio(itracer))), LBC)
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('LnudgeTCLM')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LnudgeTCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTvar)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idTvar(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTsur)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idTsur(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idTvar)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  Qout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idsurT)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idsurT(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idsurT(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Qout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idTsur)')
              Npts=load_l(Nval, Cval, NBT, Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  Qout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,60) ng
            WRITE (out,70) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,80) AttSW(ng), 'AttSW',                          &
     &            'Light attenuation of seawater (m-1).'
            WRITE (out,80) AttChl(ng), 'AttChl',                        &
     &            'Light attenuation by chlorophyll (1/(mg_Chl m-2)).'
            WRITE (out,90) PARfrac(ng), 'PARfrac',                      &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,90) Vp0(ng), 'Vp0',                              &
     &            'Eppley temperature-limited growth parameter',        &
     &            '(nondimensional).'
            WRITE (out,80) I_thNH4(ng), 'I_thNH4',                      &
     &            'Radiation threshold for nitrification (W/m2).'
            WRITE (out,80) D_p5NH4(ng), 'D_p5NH4',                      &
     &            'Half-saturation radiation for nitrification (W/m2).'
            WRITE (out,80) NitriR(ng), 'NitriR',                        &
     &            'Nitrification rate (day-1).'
            WRITE (out,90) K_NO3(ng), 'K_NO3',                          &
     &            'Inverse half-saturation for phytoplankton NO3',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,90) K_NH4(ng), 'K_NH4',                          &
     &            'Inverse half-saturation for phytoplankton NH4',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,90) K_PO4(ng), 'K_PO4',                          &
     &            'Inverse half-saturation for phytoplankton PO4',      &
     &            'uptake (1/(mmol_P m-3)).'
            WRITE (out,90) K_Phy(ng), 'K_Phy',                          &
     &            'Zooplankton half-saturation constant for ingestion', &
     &            '(mmol_N m-3)^2.'
            WRITE (out,80) Chl2C_m(ng), 'Chl2C_m',                      &
     &            'Maximum chlorophyll to carbon ratio (mg_Chl/mg_C).'
            WRITE (out,80) ChlMin(ng), 'ChlMin',                        &
     &            'Chlorophyll minimum threshold (mg_Chl/m3).'
            WRITE (out,80) PhyCN(ng), 'PhyCN',                          &
     &            'Phytoplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,80) R_P2N(ng), 'R_P2N',                          &
     &            'Phytoplankton P:N ratio (mol_P/mol_N).'
            WRITE (out,80) PhyIP(ng), 'PhyIP',                          &
     &            'Phytoplankton NH4 inhibition parameter (1/mmol_N).'
            WRITE (out,90) PhyIS(ng), 'PhyIS',                          &
     &            'Phytoplankton growth, initial slope of P-I curve',   &
     &            '(mg_C/(mg_Chl Watts m-2 day)).'
            WRITE (out,80) PhyMin(ng), 'PhyMin',                        &
     &            'Phytoplankton minimum threshold (mmol_N/m3).'
            WRITE (out,80) PhyMR(ng), 'PhyMR',                          &
     &            'Phytoplankton mortality rate (day-1).'
            WRITE (out,90) ZooAE_N(ng), 'ZooAE_N',                      &
     &            'Zooplankton nitrogen assimilation efficiency',       &
     &            '(nondimensional).'
            WRITE (out,80) ZooBM(ng), 'ZooBM',                          &
     &            'Rate for zooplankton basal metabolism (1/day).'
            WRITE (out,80) ZooCN(ng), 'ZooCN',                          &
     &            'Zooplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,80) ZooER(ng), 'ZooER',                          &
     &            'Zooplankton specific excretion rate (day-1).'
            WRITE (out,80) ZooGR(ng), 'ZooGR',                          &
     &            'Zooplankton maximum growth rate (day-1).'
            WRITE (out,80) ZooMin(ng), 'ZooMin',                        &
     &            'Zooplankton minimum threshold (mmol_N/m3).'
            WRITE (out,80) ZooMR(ng), 'ZooMR',                          &
     &            'Zooplankton mortality rate (day-1).'
            WRITE (out,80) LDeRRN(ng), 'LDeRRN',                        &
     &            'Large detritus N re-mineralization rate (day-1).'
            WRITE (out,80) LDeRRC(ng), 'LDeRRC',                        &
     &            'Large detritus C re-mineralization rate (day-1).'
            WRITE (out,80) CoagR(ng), 'CoagR',                          &
     &            'Coagulation rate (day-1).'
            WRITE (out,80) SDeRRN(ng), 'SDeRRN',                        &
     &            'Remineralization rate for small detritus N (day-1).'
            WRITE (out,80) SDeRRC(ng), 'SDeRRC',                        &
     &            'Remineralization rate for small detritus C (day-1).'
            WRITE (out,80) RDeRRN(ng), 'RDeRRN',                        &
     &            'Remineralization rate for river detritus N (day-1).'
            WRITE (out,80) RDeRRC(ng), 'RDeRRC',                        &
     &            'Remineralization rate for river detritus C (day-1).'
            WRITE (out,80) wPhy(ng), 'wPhy',                            &
     &            'Phytoplankton sinking velocity (m/day).'
            WRITE (out,80) wLDet(ng), 'wLDet',                          &
     &            'Large detritus sinking velocity (m/day).'
            WRITE (out,80) wSDet(ng), 'wSDet',                          &
     &            'Small detritus sinking velocity (m/day).'
            WRITE (out,80) pCO2air(ng), 'pCO2air',                      &
     &            'CO2 partial pressure in air (ppm by volume).'
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,110) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,100) Akt_bak(i,ng), 'Akt_bak', i,               &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) Tnudg(i,ng), 'Tnudg', i,                  &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,110) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning ON  point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning OFF point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,110) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,110) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Hout(idTvar(i),ng)) WRITE (out,120)                 &
     &              Hout(idTvar(i),ng), 'Hout(idTvar)',                 &
     &              'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Hout(idTsur(i),ng)) WRITE (out,120)                 &
     &              Hout(idTsur(i),ng), 'Hout(idTsur)',                 &
     &              'Write out tracer flux ', i,                        &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
            END IF
            IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Qout(idTvar(i),ng)) WRITE (out,120)                 &
     &              Qout(idTvar(i),ng), 'Qout(idTvar)',                 &
     &              'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Qout(idsurT(i),ng)) WRITE (out,120)                 &
     &              Qout(idsurT(i),ng), 'Qout(idsurT)',                 &
     &              'Write out surface tracer ', i,                     &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NBT
                i=idbio(itrc)
                IF (Qout(idTsur(i),ng)) WRITE (out,120)                 &
     &              Qout(idTsur(i),ng), 'Qout(idTsur)',                 &
     &              'Write out tracer flux ', i,                        &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters.
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          nl_tnu4(i,ng)=SQRT(ABS(nl_tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO
  30  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)
  40  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
  50  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Fennel Model Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  70  FORMAT (1x,i10,2x,a,t32,a)
  80  FORMAT (1p,e11.4,2x,a,t32,a)
  90  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
 100  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 120  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 130  FORMAT (10x,l1,2x,a,t32,a,1x,a)
      RETURN
      END SUBROUTINE read_BioPar
