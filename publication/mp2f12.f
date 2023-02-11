      program mp2f12
#include "MRCCCOMMON"
      integer nbfc,ifocka,ifockb,iexcha,iexchb,icmoa,icmob,ieva,ievb
      integer dblalloc,verbosity,ihpja,ihpjb
      character(len=4) c4
      character(len=5) scftype
      character(len=6) core
C
      call mrccini
      write(iout,*) 'MP2-F12 calculation'
      write(iout,*)
      call memalloc
C
      call getvar('nbf       ',nbf)
      call getvar('nal       ',nal)
      call getvar('nbe       ',nbe)
      call getvar('ncore     ',ncore)
      call getkey('core',4,core,6)
      call getkey('scftype',7,scftype,5)
      call getkey('verbosity',9,c4,4)
      read(c4,"(i4)") verbosity
       if(core.eq.'corr  ') ncore=0
      nocc=nal+nbe
      nvirtal=nbf(5)-nal
      nvirtbe=nbf(5)-nbe
      nvirt=nvirtal+nvirtbe
      nbasis=nbf(5)
      nbfc=nbf(7)
      dfnbasis=nbf(3)
C
      ifocka=dblalloc(nbfc*nbfc)
      iexcha=dblalloc(nbfc*nbfc)
      ihpja =dblalloc(nbfc*nbfc)
      icmoa =dblalloc(nbfc*nbfc)
      ieva  =dblalloc(nbfc)
       if(scftype.ne.'rhf  ') then
        ifockb=dblalloc(nbfc*nbfc)
        iexchb=dblalloc(nbfc*nbfc)
        ihpjb =dblalloc(nbfc*nbfc)
        icmob =dblalloc(nbfc*nbfc)
        ievb  =dblalloc(nbfc)
       else
        ifockb=ifocka
        iexchb=iexcha
        ihpjb =ihpja
        icmob =icmoa
        ievb  =ieva
       endif
C
      call mp2f12core(nbasis,nbf(6),nbfc,scrfile1,nocc,nvirt,nal,nbe,
     $nvirtal,nvirtbe,ncore,scftype,iout,dcore(ifocka),dcore(ifockb),
     $dcore(iexcha),dcore(iexchb),dcore(icmoa),dcore(icmob),
     $dcore(ieva),dcore(ievb),dcore,imem,verbosity,diisfile,errfile,
     $ifltln,dcore(ihpja),dcore(ihpjb),ifcfile,varsfile,icore,dfnbasis,
     $nbf,maxcor,scrfile2)
C
      call mrccend(0)
      end
C
************************************************************************
      subroutine mp2f12core(nbasis,nbfa,nbfc,scrfile1,nocc,nvirt,nal,
     $nbe,nvirtal,nvirtbe,ncore,scftype,iout,focka,fockb,excha,exchb,
     $cmoa,cmob,eva,evb,dcore,imem,verbosity,diisfile,errfile,ifltln,
     $hpja,hpjb,ifcfile,varsfile,icore,dfnbasis,nbf,maxcor,scrfile2)
************************************************************************
* Driver for MP2-F12 calculations
************************************************************************
      implicit none
      integer nbfc,scrfile1,nbasis,nbfa,i,j,nbfan,nbfcn,iout,nalc,nbec
      integer nocc,nvirt,nal,nbe,nvirtal,nvirtbe,ncore,imem,dblalloc,nn
      integer nbf(*),iepaa,iepbb,iepab,ifcfile,if12ij
      integer verbosity,diisfile,errfile,ifltln,varsfile,imem1
      integer icf12aa,icf12bb,icf12ab,icf12ba,nb,itscalea,itscaleb
      integer icore(*),isinv,is,is12,iconv,intalloc,dfnbasis,maxcor
      integer icf122aa,icf122bb,icf122ab,icf122ba,scrfile2
      integer ivaa,ivbb,ivab
      integer ivifaa2,icf12aa2,nel,nlm
      integer ivpia,ivpib,icaia,icaib,iuaia,iuaib
      integer icf12r12aa,icf12r12bb,icf12r12ab,icf12r12ba
      real*8 focka(nbfc*nbfc),fockb(nbfc*nbfc),eva(nbfc),evb(nbfc),eref
      real*8 excha(nbfc*nbfc),exchb(nbfc*nbfc),dcore(*),emp2f12,ecabs
      real*8 hpja(nbfc*nbfc),hpjb(nbfc*nbfc),cctol,tfact,ecoup
      real*8 cmoa(nbfc,nbfc),cmob(nbfc,nbfc)
      real*8 ovltol
      character(len=4) c4
      character(len=5) scftype
      character(len=16) c1
      character(len=32) calc
C Variables for integral calculation
      integer natoms,nangmax,ncontrmax,nprimmax,ncartmax,nsphermax
      integer nmboys,nbfshmax,inang,incontr,inprim,igexp,igcoef,icoord
      integer ictostr,icf,iboysval,iindarr,igcn,ipre,inshrange,inangmin
      integer iatnum,inzipr,inzjpr,inzkpr,inzlpr,inzint,ithad,ithcf2
      integer iscoord,irqqij,irqqkl,ihrec,nbset,inbf,ispctostr,iatchg
      integer ncent,iebf,necpatoms,nbfatmax,ispre,idfipre,itcdfint
      integer idfrqq,inatrange,idfipra,ilogkc,igck,ikp,icprea
      integer ivifaa,ivifbb,ivifab,ivifba,itcf12int,icpreb,iint1
      integer freememspace
      real*8 itol
      logical cartg,lcc,lblock
      common/memcom/ imem1
C
      call getvar('eref      ',eref)
      call getkey('calc',4,calc,32)
      call getkey('ovltol',6,c1,16)
      read(c1,*) ovltol
      nalc=nal-ncore
      nbec=nbe-ncore
      lcc=trim(calc).ne.'mp2-f12'
C nbasis -size of AO basis
C nbfa - size of CABS basis
C nbfc - size of merged AO + CABS basis
C nbfan - size of non-redundant CABS basis
C nbfcn - size of non-redundant AO + CABS basis
C Construct CABS MOs
      isinv=dblalloc(nbasis*nbasis)
      is=dblalloc(nbfc*nbfc)
      is12=dblalloc(nbasis*nbfa)
      iconv=intalloc(nbfc)
      call getvar('conv      ',icore(iconv))
      open(scrfile1,file='S12MATCABS',form='unformatted')
      call rtdmx(dcore(imem),dcore(imem),dcore(is12),scrfile1,nbasis,
     $nbfa)
      call roeint(dcore(imem),dcore(imem),dcore(is),scrfile1,nbfc)
      close(scrfile1)
      open(scrfile1,file='SROOT_AO',form='unformatted')
      read(scrfile1)
      call roeint(dcore(imem),dcore(imem),dcore(isinv),scrfile1,nbasis)
      close(scrfile1)
      call dsymm('l','l',nbasis,nbasis,1.d0,dcore(isinv),nbasis,
     $dcore(isinv),nbasis,0.d0,dcore(imem),nbasis)
      call dgemm('n','n',nbasis,nbfa,nbasis,1.d0,dcore(imem),nbasis,
     $dcore(is12),nbasis,0.d0,dcore(imem+nbasis**2),nbasis)
      call dcopy(nbasis*nbfa,dcore(imem+nbasis**2),1,dcore(is12),1)
      call aoproj(nbfc,nbfa,nbasis,icore(iconv),cmoa,dcore(is12))
      call dsymm('l','u',nbfc,nbfa,1.d0,dcore(is),nbfc,cmoa(1,nbasis+1),
     $nbfc,0.d0,dcore(imem+nbfc*nbfc),nbfc)
      call dgemm('t','n',nbfa,nbfa,nbfc,1.d0,cmoa(1,nbasis+1),nbfc,
     $dcore(imem+nbfc*nbfc),nbfc,0.d0,dcore(imem),nbfa)
      call dsyev('V','U',nbfa,dcore(imem),nbfa,eva,
     $dcore(imem+nbfa*nbfa),3*nbfa**2,i)
      nbfan=0
       do while(eva(nbfan+1).lt.ovltol.and.(nbfan+1).le.nbfa)
        nbfan=nbfan+1
       enddo
       do j=nbfan+1,nbfa
        call dscal(nbfa,1.d0/dsqrt(eva(j)),dcore(imem+(j-1)*nbfa),1)
       enddo
      call dgemm('n','n',nbfc,nbfa-nbfan,nbfa,1.d0,cmoa(1,nbasis+1),
     $nbfc,dcore(imem+nbfan*nbfa),nbfa,0.d0,dcore(imem+nbfa*nbfa),nbfc)
      nbfan=nbfa-nbfan
      nbfcn=nbfc-(nbfa-nbfan)
      call dcopy(nbfc*nbfan,dcore(imem+nbfa*nbfa),1,cmoa(1,nbasis+1),1)
       if(scftype.ne.'rhf  ') cmob=cmoa
C Read HF MO coefficients
      open(scrfile1,file='MOCOEF',form='unformatted')
      call moread(nbfc,nbasis,scrfile1,icore(iconv),cmoa,dcore(imem),
     $dcore(imem+nbasis**2))
       if(scftype.ne.'rhf  ') call moread(nbfc,nbasis,scrfile1,
     $icore(iconv),cmob,dcore(imem),dcore(imem+nbasis**2))
      close(scrfile1)
      call dbldealloc(isinv)
C Read Fock-matrix
      open(scrfile1,file='FOCK',form='unformatted')
      call roeint(dcore(imem),dcore(imem),focka,scrfile1,nbfc)
       if(scftype.ne.'rhf  ') then
        call roeint(dcore(imem),dcore(imem),fockb,scrfile1,nbfc)
        read(scrfile1)
       endif
      read(scrfile1)
      call roeint(dcore(imem),dcore(imem),excha,scrfile1,nbfc)
       if(scftype.ne.'rhf  ')
     $call roeint(dcore(imem),dcore(imem),exchb,scrfile1,nbfc)
      close(scrfile1)
C Canonicalize CABS virtuals
      call cancabs(nbfc,nbfan,nbasis,focka,cmoa,dcore(imem+nbfan*nbfan),
     $eva,dcore(imem),iout)
       if(scftype.ne.'rhf  ')
     $call cancabs(nbfc,nbfan,nbasis,fockb,cmob,dcore(imem+nbfan*nbfan),
     $evb,dcore(imem),iout)
C Transform Fock- and exchange matrix to CABS MO basis
      call transf(nbfc,nbfcn,focka,cmoa,dcore(imem),eva,.true.)
      call transf(nbfc,nbfcn,excha,cmoa,dcore(imem),eva,.false.)
      call dcopy(nbfcn*nbfcn,focka,1,hpja,1)
      call daxpy(nbfcn*nbfcn,-1.d0,excha,1,hpja,1)
       if(scftype.ne.'rhf  ') then
        call transf(nbfc,nbfcn,fockb,cmob,dcore(imem),evb,.true.)
        call transf(nbfc,nbfcn,exchb,cmob,dcore(imem),evb,.false.)
        call dcopy(nbfcn*nbfcn,fockb,1,hpjb,1)
        call daxpy(nbfcn*nbfcn,-1.d0,exchb,1,hpjb,1)
       endif
       if(verbosity.ge.3) then
        write(iout,*)
        write(iout,"(' Number of alpha electrons:    ',14x,i6)") nal
        write(iout,"(' Number of beta  electrons:    ',14x,i6)") nbe
        write(iout,"(' Number of core orbitals:      ',14x,i6)") ncore
        write(iout,"(' Number of AOs:                ',14x,i6)") nbasis
        write(iout,"(' Number of auxiliary functions:',14x,i6)") nbfa
        write(iout,
     $"(' Number of non-redundant auxiliary functions:',i6)") nbfan
        write(iout,
     $"(' Number of non-redundant CABS functions:     ',i6)") nbfcn
       endif
C Calculate CABS singles correction
      write(iout,*)
      write(iout,*) 'Calculating CABS singles correction...'
      call getkey('cctol',5,c4,4)
      read(c4,*) i
      cctol=10.d0**(-i)
      ecabs=0d0
      call cabssing(nbfcn,nbfan,nal,nvirtal,ecabs,iout,focka,eva,
     $eva(nal+1),diisfile,errfile,ifltln,cctol,dcore(imem),
     $dcore(imem+nbfan*nal),dcore(imem+2*nbfan*nal),verbosity)
       if(scftype.eq.'rhf  ') then
        ecabs=2.d0*ecabs
       else
        call cabssing(nbfcn,nbfan,nbe,nvirtbe,ecabs,iout,fockb,evb,
     $evb(nbe+1),diisfile,errfile,ifltln,cctol,dcore(imem),
     $dcore(imem+nbfan*nbe),dcore(imem+2*nbfan*nbe),verbosity)
       endif
      call timer
c     write(iout,*)
c     write(iout,"(' Reference energy [au]:           'f24.12)") eref
c     write(iout,"(' CABS singles correction [au]:    'f24.12)") ecabs
c     write(iout,"(' Corrected reference energy [au]: 'f24.12)")
c    $eref+ecabs
c     call mrccend(0)
C Allocate memory for integrals, fitting coefficients, pair energies
      nb=nbasis-ncore
      nn=nalc
       if(lcc) nn=nb

C TODO calculate this
      nlm = nalc
      if(nlm.lt.nalc)then
       write(iout,*)'nlm.lt.nalc in mp2f12core... fix this!'
      stop
      endif
      

      freememspace = maxcor-imem+imem1
      write(iout,"('M.left (Mb):',14x,i6)")(8*freememspace)/(1024*1024)

      itscalea=dblalloc(nalc*2)
      freememspace = maxcor-imem+imem1

      write(iout,"('Bottleneck G/C size (Mb):    ',14x,i6)")
     $(8*nbfcn*dfnbasis*nn)/(1024*1024)

      if(freememspace.gt.nbfcn*dfnbasis*nn*3)then
      ivifaa  =dblalloc(nn  *nbfcn*dfnbasis)
      icf12aa =dblalloc(nalc*nbfcn*dfnbasis)
      lblock = .False.
      else
       write(iout,*)'Only a block of G/C is allocated...'
       lblock = .True.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write(iout,"('Memalloc attempt for G block (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*(nlm+1))/(1024*1024)
      ivifaa  =dblalloc((nlm+1)*nbfcn*dfnbasis)
      icf12aa =dblalloc((nlm+1)*nbfcn*dfnbasis)
      endif

       lblock = .True.


      iepaa=dblalloc(nalc*(nalc+1)/2)
       if(scftype.ne.'rhf  ') then
        itscaleb=dblalloc(nbec*2)
        ivifbb  =dblalloc(nn  *nbfcn*dfnbasis)
        icf12bb =dblalloc(nalc*nbfcn*dfnbasis)
        ivifab  =dblalloc(nn  *nbfcn*dfnbasis)
        icf12ab =dblalloc(nalc*nbfcn*dfnbasis)
        ivifba  =dblalloc(nn  *nbfcn*dfnbasis)
        icf12ba =dblalloc(nalc*nbfcn*dfnbasis)
        iepbb=dblalloc(nbec*(nbec+1)/2)
        iepab=dblalloc(nalc*nbec)
       else
        itscaleb=itscalea
        ivifbb=  ivifaa
        icf12bb= icf12aa
        ivifab=  ivifaa
        icf12ab= icf12aa
        ivifba=  ivifaa
        icf12ba= icf12aa
        iepbb=iepaa
        iepab=iepaa
       endif
       if(lcc) then

       write(iout,*)'Memory allocation for icf12r12aa: '
       write(iout,*)'dfnbasis: ',dfnbasis
       write(iout,*)'nalc: ',nalc
       write(iout,*)'nlm: ',nlm
       write(iout,*)'nb: ',nb
       write(iout,*)'nbec: ',nbec
       write(iout,*)'nb*nalc*dfnbasis: ',nb*nalc*dfnbasis
       write(iout,*)'nbfan*nvirtbe*nalc: ',nbfan*nvirtbe*nalc

        icf12r12aa=dblalloc(max(nb*nalc*dfnbasis,nbfan*nvirtbe*nalc))
         if(scftype.ne.'rhf  ') then
          icf12r12bb=dblalloc(max(nb*nbec*dfnbasis,nbfan*nvirtbe*nalc))
          icf12r12ab=dblalloc(max(nb*nalc*dfnbasis,nbfan*nvirtbe*nalc))
          icf12r12ba=dblalloc(max(nb*nbec*dfnbasis,nbfan*nvirtbe*nalc))
         else
          icf12r12bb=icf12r12aa
          icf12r12ab=icf12r12aa
          icf12r12ba=icf12r12aa
         endif
       else
        icf12r12aa=imem
        icf12r12bb=imem
        icf12r12ab=imem
        icf12r12ba=imem
       endif
C Read data for integral calculation
      open(varsfile,file='VARS',form='UNFORMATTED')
      call tedatr(varsfile,natoms,nangmax,ncontrmax,nprimmax,
     $ncartmax,cartg,nsphermax,nmboys,i,itol,nbfshmax,
     $inang,incontr,inprim,igexp,igcoef,icoord,ictostr,icf,iboysval,
     $iindarr,igcn,ipre,inshrange,inangmin,iatnum,icore,dcore,inzipr,
     $inzjpr,inzkpr,inzlpr,inzint,ithad,ithcf2,iscoord,irqqij,irqqkl,
     $ihrec,iout,0,nbset,inbf,ispctostr,iatchg,ncent,iebf,necpatoms,
     $i,i,i,i,i,i,.false.,nbfatmax,ispre)
      close(varsfile)
C Allocate memory for integral calculation
      idfipre=dblalloc((nangmax+1)*natoms)
      itcdfint=dblalloc((dfnbasis+1)*dfnbasis/2)
      itcf12int=dblalloc(dfnbasis*dfnbasis)
      idfrqq=dblalloc((4*nangmax+1)*nprimmax*(nangmax+1)*natoms)
      inatrange=intalloc(2*natoms*nbset)
      call getvar('natrange  ',icore(inatrange))
      idfipra=dblalloc(natoms)
      ilogkc=intalloc((nangmax+1)*natoms)
      igck=dblalloc(nprimmax*ncontrmax*(nangmax+1)*natoms)
      ikp=intalloc(nprimmax*(nangmax+1)*natoms)
      icprea=dblalloc(natoms*(nangmax+1))
       if(scftype.ne.'rhf  ') then
        icpreb=dblalloc(natoms*(nangmax+1))
       else
        icpreb=icprea
       endif
C Calculate and transform three-center integrals for MP2-F12
      icf122aa=dblalloc(nalc*nbfcn*dfnbasis)
       if(scftype.ne.'rhf  ') then
        icf122bb=dblalloc(nbec*nbfcn*dfnbasis)
        icf122ab=dblalloc(nalc*nbfcn*dfnbasis)
        icf122ba=dblalloc(nbec*nbfcn*dfnbasis)
       else
        icf122bb=icf122aa
        icf122ab=icf122aa
        icf122ba=icf122aa
       endif


      if(lblock)then
       write(iout,*)"Executing block mp2int..."
      call mp2int_bl(nbf(1),nalc,nbec,ncore,nb,dfnbasis,nbfc,nbfcn,cmoa,
     $cmob,dcore(iepaa),dcore(iepbb),dcore(iepab),focka,fockb,hpja,hpjb,
     $dcore,imem,natoms,iout,imem1,maxcor,nangmax,ncontrmax,nprimmax,
     $ncartmax,cartg,nsphermax,nmboys,itol,nbfshmax,icore(inang),
     $icore(incontr),icore(inprim),dcore(igexp),dcore(igcoef),
     $dcore(icoord),dcore(ictostr),dcore(icf),dcore(iboysval),
     $icore(iindarr),icore(igcn),dcore(ipre),icore(inshrange),
     $icore(ithad),dcore(ithcf2),dcore(iscoord),dcore(irqqij),
     $dcore(irqqkl),dcore(ihrec),nbset,dcore(ispctostr),dcore(idfipre),
     $icore,dcore(itcdfint),dcore(idfrqq),dcore(ispre),icore(inatrange),
     $dcore(idfipra),icore(ilogkc),icore(ikp),dcore(igck),dcore(icprea),
     $dcore(icpreb),dcore(itcf12int),dcore(ivifaa),dcore(ivifbb),
     $dcore(ivifab),dcore(icf12aa),dcore(icf12bb),dcore(icf12ab),
     $dcore(icf122aa),dcore(icf122bb),dcore(icf122ab),dcore(ivifba),
     $dcore(icf12ba),dcore(icf122ba),scrfile2,scftype,lcc,
     $dcore(icf12r12aa),dcore(icf12r12bb),dcore(icf12r12ab),
     $dcore(icf12r12ba),verbosity,ifltln)
      else
      write(iout,*)"Executing in-memory mp2int..."
      call mp2int(nbf(1),nalc,nbec,ncore,nb,dfnbasis,nbfc,nbfcn,cmoa,
     $cmob,dcore(iepaa),dcore(iepbb),dcore(iepab),focka,fockb,hpja,hpjb,
     $dcore,imem,natoms,iout,imem1,maxcor,nangmax,ncontrmax,nprimmax,
     $ncartmax,cartg,nsphermax,nmboys,itol,nbfshmax,icore(inang),
     $icore(incontr),icore(inprim),dcore(igexp),dcore(igcoef),
     $dcore(icoord),dcore(ictostr),dcore(icf),dcore(iboysval),
     $icore(iindarr),icore(igcn),dcore(ipre),icore(inshrange),
     $icore(ithad),dcore(ithcf2),dcore(iscoord),dcore(irqqij),
     $dcore(irqqkl),dcore(ihrec),nbset,dcore(ispctostr),dcore(idfipre),
     $icore,dcore(itcdfint),dcore(idfrqq),dcore(ispre),icore(inatrange),
     $dcore(idfipra),icore(ilogkc),icore(ikp),dcore(igck),dcore(icprea),
     $dcore(icpreb),dcore(itcf12int),dcore(ivifaa),dcore(ivifbb),
     $dcore(ivifab),dcore(icf12aa),dcore(icf12bb),dcore(icf12ab),
     $dcore(icf122aa),dcore(icf122bb),dcore(icf122ab),dcore(ivifba),
     $dcore(icf12ba),dcore(icf122ba),scrfile2,scftype,lcc,
     $dcore(icf12r12aa),dcore(icf12r12bb),dcore(icf12r12ab),
     $dcore(icf12r12ba),verbosity,ifltln)
      endif

      call dbldealloc(inang)
C Calculate MP2-F12 pair energies

C legyen ez most itt 1 és akkor mindenből csak 1 fér be a memóriába 

      nlm = 15 
      call mp2pair(nbfcn,ncore,nal,nbe,nalc,nbec,nbasis,dfnbasis,nbfan,
     $iout,scftype,verbosity,eref,emp2f12,ecabs,ifcfile,dcore,imem,
     $focka,fockb,eva,evb,eva(nal+1),evb(nbe+1),nvirtal,nvirtbe,excha,
     $exchb,dcore(ivifaa),dcore(ivifbb),dcore(ivifab),dcore(ivifba),
     $dcore(icf12aa),dcore(icf12bb),dcore(icf12ab),dcore(icf12ba),
     $dcore(iepaa),dcore(iepbb),dcore(iepab),dcore(itscalea),
     $dcore(itscaleb),tfact,cctol,ecoup,maxcor,imem1,ifltln,lblock,nlm)
       if(trim(calc).eq.'mp2-f12') return
C Calculate CC intermediates
      if12ij=dblalloc(max(nbfcn*nbfcn,2*nbfan*nvirtbe))
c     write(iout,*)'alloc1: ',max(nbfcn*nbfcn,2*nbfan*nvirtbe)
c     write(iout,*)'alloc2: ',nbfcn*nbfcn
c     write(iout,*)'alloc3: ',2*nbfan*nvirtbe
c     write(iout,*)'start add: ',loc(dcore(if12ij))

      ivpia=dblalloc(nb*nalc)
      icaia=dblalloc(nvirtal*nalc)
      iuaia=dblalloc(nvirtal*nalc)
       if(scftype.eq.'rhf  ') then
        ivaa=dblalloc(nb*nb*(nalc+1)*nalc/2)
        ivbb=ivaa
        ivab=ivaa
        ivpib=ivpia
        icaib=icaia
        iuaib=iuaia
       else
        ivaa=dblalloc((nb-1)*nb/2*(nalc-1)*nalc/2)
        ivbb=dblalloc((nb-1)*nb/2*(nbec-1)*nbec/2)
        ivab=dblalloc(nb*nb*nalc*nbec)
        ivpib=dblalloc(nb*nbec)
        icaib=dblalloc(nvirtbe*nbec)
        iuaib=dblalloc(nvirtbe*nbec)
       endif

C Call ccimed  block version if vifaa too large
      if(.true.)then
C legyen ez most itt 1 és akkor mindenből csak 1 fér be a memóriába 
C ha nbfcn*dfnbasis*nalc befér a memóriába akkor csak a cabijcalc_bl belső ciklusát kell blokkosítani
      nlm = 15
      ivifaa2  =dblalloc((nlm+1)*nbfcn*dfnbasis)
      icf12aa2 =dblalloc((nlm+1)*nbfcn*dfnbasis)

      call ccimed_bl(nb,dfnbasis,nalc,nbec,nvirtal,nvirtbe,ncore,focka,
     $fockb,nbfcn,nbfan,emp2f12,ecabs,scrfile1,scftype,dcore(ivaa),
     $dcore(ivaa),dcore(ivbb),dcore(ivab),dcore(ivpia),dcore(ivpib),
     $dcore(icaia),dcore(icaib),dcore(iuaia),dcore(iuaib),dcore(ivifaa),
     $dcore(ivifbb),dcore(ivifab),dcore(ivifba),dcore(icf12aa),
     $dcore(icf12bb),dcore(icf12ab),dcore(icf12ba),iout,dcore,imem,
     $imem1,maxcor,dcore(if12ij),dcore(if12ij),dcore(if12ij),
     $dcore(if12ij+nbfan*nvirtbe),verbosity,dcore(icf12r12aa),
     $dcore(icf12r12bb),dcore(icf12r12ab),dcore(icf12r12ba),
     $dcore(itscalea),dcore(itscaleb),tfact,dcore(iepaa),dcore(iepbb),
     $dcore(iepab),ecoup,ifltln,dcore(ivifaa2),dcore(icf12aa2),nlm)

      else
      call ccimed(nb,dfnbasis,nalc,nbec,nvirtal,nvirtbe,ncore,focka,
     $fockb,nbfcn,nbfan,emp2f12,ecabs,scrfile1,scftype,dcore(ivaa),
     $dcore(ivaa),dcore(ivbb),dcore(ivab),dcore(ivpia),dcore(ivpib),
     $dcore(icaia),dcore(icaib),dcore(iuaia),dcore(iuaib),dcore(ivifaa),
     $dcore(ivifbb),dcore(ivifab),dcore(ivifba),dcore(icf12aa),
     $dcore(icf12bb),dcore(icf12ab),dcore(icf12ba),iout,dcore,imem,
     $imem1,maxcor,dcore(if12ij),dcore(if12ij),dcore(if12ij),
     $dcore(if12ij+nbfan*nvirtbe),verbosity,dcore(icf12r12aa),
     $dcore(icf12r12bb),dcore(icf12r12ab),dcore(icf12r12ba),
     $dcore(itscalea),dcore(itscaleb),tfact,dcore(iepaa),dcore(iepbb),
     $dcore(iepab),ecoup)
      endif
C
      return
      end
C
************************************************************************
      subroutine mp2pair(nbfcn,ncore,nal,nbe,nalc,nbec,nbasis,dfnbasis,
     $nbfan,iout,scftype,verbosity,eref,ef12,ecabs,ifcfile,dcore,imem,
     $focka,fockb,eoa,eob,eva,evb,nvirtal,nvirtbe,excha,exchb,vifaa,
     $vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,epaa,epbb,epab,
     $tscalea,tscaleb,tfact,cctol,ecoup,maxcor,imem1,ifltln,lblock,nlm)
************************************************************************
* Calculate MP2-F12 energy
************************************************************************
      implicit none
      integer nbfcn,nocc,nvirt,i,ncore,nal,nbe,nbasis,dfnbasis,nbfan
      integer verbosity,nvirtal,nvirtbe,ifcfile,iout,nalc,nbec,ivifaaf
      integer iscr,irsa,if12ij,if12sy,ivintpq,ivintao,ivintoa,icf12aaf
      integer dblalloc,irsb,icf12bbf,icf12abf,icf12baf,ivifbbf
      integer ivifabf,ivifbaf,imp2imed
      integer bsize,ifltln,nblocks,nlm
      integer ivifaab, icf12aab, icf12aafb, ivifaafb

      integer maxcor,imem1,imem,freememspace
      logical lblock

      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn),eoa(*),eob(*),eva(*)
      real*8 excha(nbfcn,nbfcn),exchb(nbfcn,nbfcn),dcore(*),evb(*),tfact
      real*8 vifaa,vifbb,vifab,cf12aa,cf12bb,cf12ab,epab(*),ecoup,cctol
      real*8 emp2,ecabs,ef12,et1,eref,vifba,cf12ba,epaa(*),epbb(*)
      real*8 tscalea(nalc,2),tscaleb(nbec,2)
      character(len=5) scftype
C Calculate MP2 singles contribution
      et1=0.d0
      tscalea=0.d0
      tscaleb=0.d0
      call mp2sing(nbfcn,ncore,nal,nvirtal,focka,eoa,eva,et1,
     $tscalea(1,2))
       if(scftype.eq.'rhf  ') then
        et1=2.d0*et1
        call dscal(nalc,2.d0,tscalea(1,2),1)
       else
        call mp2sing(nbfcn,ncore,nbe,nvirtbe,fockb,eob,evb,et1,
     $tscaleb(1,2))
       endif
      emp2=et1
C Calculate F12 contribution
      write(iout,*)
      write(iout,*) 'Calculating pair energies...'

c     write(iout,*) 'ncore before file: ',ncore
      open(imp2imed,status='unknown',file='MP2IMED',form='unformatted')
      write(imp2imed)ncore
      close(imp2imed)

      open(imp2imed,status='unknown',file='MP2IMED',form='unformatted')
      read(imp2imed)ncore
c     write(iout,*) 'ncore from file: ',ncore
      close(imp2imed)


      iscr    =dblalloc(nbfcn*nbfcn)
      irsa    =dblalloc(nbfcn*nbfcn)
      if12ij  =dblalloc(nbfcn*nbfcn)
      if12sy  =dblalloc(nbfcn*nbfcn)
      ivintpq =dblalloc(nbasis*nbasis)
      ivintao =dblalloc(nbfan*nal)
      ivintoa =dblalloc(nal*nbfan)
      icf12aaf=dblalloc(nbfcn*dfnbasis*nalc)
c     freememspace = 751775
      freememspace = maxcor-imem+imem1
c     write(iout,"(' nbfcn*dfnbasis:    ',14x,i8)") nbfcn*dfnbasis
c     write(iout,"(' nbfcn:    ',14x,i6)") nbfcn
c     write(iout,"(' nalc:    ',14x,i6)") nalc

      if (lblock)then
       write(iout,*)
     $'Warning, mp2pair requires more mem., disk will be used!'
       bsize = nlm
       if(bsize.eq.0)then
       write(iout,*)
     $  'Not enough mem. even for a single block...'
        call mrccend(1)
       endif
       write(iout,"(' nalc:    ',14x,i6)") nalc
       write(iout,"(' bsize:    ',14x,i6)") bsize
       do while(mod(nalc,bsize).ne.0)
        bsize = bsize - 1
       enddo
       nblocks = nalc / bsize
       write(iout,"(' block size:    ',14x,i6)") bsize
       write(iout,"(' nblocks:    ',14x,i6)") nblocks

       icf12aaf=dblalloc(nbfcn*dfnbasis*(bsize+1))
       ivifaaf =dblalloc(nbfcn*dfnbasis*(bsize+1))
      else
       icf12aaf=dblalloc(nbfcn*dfnbasis*nalc)
       ivifaaf =dblalloc(nbfcn*dfnbasis*nalc)
      endif
      if(.not.lblock) write(iout,*)
     $'mp2pair executes in-memory...'


      write(iout,"(' G/C size (Mb):    ',14x,i6)") 
     $(8*nbfcn*dfnbasis*nalc)/(1024*1024)
c     icf12aafb = dblalloc(nbfcn*dfnbasis*(bsize+1))
c     ivifaafb  = dblalloc(nbfcn*dfnbasis*(bsize+1))

c     write(iout,"(' iscr:    ',14x,i6)") nbfcn*nbfcn
c     write(iout,"(' irsa:    ',14x,i6)") nbfcn*nbfcn
c     write(iout,"(' if12ij:    ',14x,i6)") nbfcn*nbfcn
c     write(iout,"(' if12sy:    ',14x,i6)") nbfcn*nbfcn
c     write(iout,"(' ivintpq:    ',14x,i6)") nbasis*nbasis
c     write(iout,"(' ivintao:    ',14x,i6)") nbfan*nal
c     write(iout,"(' ivintoa:    ',14x,i6)") nal*nbfan
c     write(iout,"(' nalc:    ',14x,i6)") nalc

      ef12=0.d0
      ecoup=0.d0
       if(verbosity.gt.2) write(iout,*)
     $'    MOs     MP2 pair energy [au]  F12 pair energy [au]'
       if(scftype.eq.'rhf  ') then
C Closed-shell
      write(iout,"(' 1. F12 contribution [au]: 'f24.12)") ef12

C ezt át kell gondolni hogy itt mit fogok átadni dcore(ivifaab),dcore(icf12aab) helyeken: semmit?
      if(lblock)then

C Bence módosítások
C intermediers for C_abij are calculated from density-fitted G and C (mem. heavy) matrices as
C Q^{ab}_{ij} = G^{P}_{ia}*C^{P}_{jb} + C^{P}_{ia}*G^{P}_{jb}
C allocate mem. for the blocks of G and C
C bsize nbfcn*dfnbasis blocks fit in the mem. 
C one additional element is allocated for index "j"
      ivifaab   = dblalloc(nbfcn*dfnbasis*(bsize+1))
      icf12aab  = dblalloc(nbfcn*dfnbasis*(bsize+1))

       call epair_st_bl(nal,nalc,nbfcn,nbfan,nvirtal,ncore,nbasis,
     $dfnbasis,focka,excha,eoa,eva,iout,verbosity,ef12,emp2,ecoup,
     $dcore(ivifaab),dcore(icf12aab),epaa,dcore(irsa),
     $dcore(if12ij),dcore(if12sy),
     $dcore(iscr),dcore(ivintpq),dcore(ivintao),dcore(ivintoa),
     $dcore(icf12aaf),dcore(ivifaaf),tscalea,cctol,bsize,ifltln,nblocks)

      else
       call epair_st(nal,nalc,nbfcn,nbfan,nvirtal,ncore,nbasis,
     $dfnbasis,focka,excha,eoa,eva,iout,verbosity,ef12,emp2,ecoup,
     $vifaa,cf12aa,epaa,dcore(irsa),dcore(if12ij),dcore(if12sy),
     $dcore(iscr),dcore(ivintpq),dcore(ivintao),dcore(ivintoa),
     $dcore(icf12aaf),dcore(ivifaaf),tscalea,cctol)
      endif

      write(iout,"(' 2. F12 contribution [au]: 'f24.12)") ef12
       else
C alpha-alpha
        call epair_ss(nal,nalc,nbfcn,nbfan,nvirtal,ncore,nbasis,
     $dfnbasis,focka,excha,eoa,eva,iout,verbosity,ef12,emp2,ecoup,
     $vifaa,cf12aa,epaa,dcore(irsa),dcore(if12ij),dcore(if12sy),
     $dcore(iscr),dcore(ivintpq),dcore(ivintao),dcore(ivintoa),
     $dcore(icf12aaf),dcore(ivifaaf),'a',tscalea,cctol)
C beta-beta
        call epair_ss(nbe,nbec,nbfcn,nbfan,nvirtbe,ncore,nbasis,
     $dfnbasis,fockb,exchb,eob,evb,iout,verbosity,ef12,emp2,ecoup,
     $vifbb,cf12bb,epbb,dcore(irsa),dcore(if12ij),dcore(if12sy),
     $dcore(iscr),dcore(ivintpq),dcore(ivintao),dcore(ivintoa),
     $dcore(icf12aaf),dcore(ivifaaf),'b',tscaleb,cctol)
C alpha-beta
      irsb    =dblalloc(nbfcn*nbfcn)
      icf12bbf=dblalloc(nbfcn*dfnbasis*nbec)
      icf12abf=dblalloc(nbfcn*dfnbasis*nalc)
      icf12baf=dblalloc(nbfcn*dfnbasis*nbec)
      ivifbbf =dblalloc(nbfcn*dfnbasis*nbec)
      ivifabf =dblalloc(nbfcn*dfnbasis*nalc)
      ivifbaf =dblalloc(nbfcn*dfnbasis*nbec)
        call epair_os(nal,nbe,nalc,nbec,nbfcn,nbfan,nvirtal,nvirtbe,
     $ncore,nbasis,dfnbasis,focka,fockb,excha,exchb,eoa,eva,eob,evb,
     $iout,verbosity,ef12,emp2,ecoup,vifaa,vifbb,vifab,vifba,cf12aa,
     $cf12bb,cf12ab,cf12ba,epab,dcore(irsa),dcore(irsb),dcore(if12ij),
     $dcore(if12sy),dcore(iscr),dcore(ivintpq),dcore(ivintao),
     $dcore(ivintoa),dcore(icf12aaf),dcore(icf12bbf),dcore(icf12abf),
     $dcore(icf12baf),dcore(ivifaaf),dcore(ivifbbf),dcore(ivifabf),
     $dcore(ivifbaf),tscalea,tscaleb,cctol)
       endif
      call timer
C Calculating scaling parameters for (T+)
       if(verbosity.gt.2) then
        write(iout,*)
        write(iout,*) 'Contribution of orbitals to correlation energy:'
        write(iout,*)
     $' MO          MP2 [au]              F12 [au]           ' //
     $'(MP2+F12)/MP2'
         if(scftype.eq.'rhf  ') then
           do i=1,nalc
             if(dabs(tscalea(i,2)).gt.cctol)
     $write(iout,"(i4,1x,f20.12,2x,f20.12,2x,f20.12)") ncore+i,
     $tscalea(i,2),tscalea(i,1),(tscalea(i,1)+tscalea(i,2))/tscalea(i,2)
           enddo
         else
           do i=1,nalc
             if(dabs(tscalea(i,2)).gt.cctol)
     $write(iout,"(i4,a1,f20.12,2x,f20.12,2x,f20.12)") ncore+i,'a',
     $tscalea(i,2),tscalea(i,1),(tscalea(i,1)+tscalea(i,2))/tscalea(i,2)
           enddo
           do i=1,nbec
             if(dabs(tscaleb(i,2)).gt.cctol)
     $write(iout,"(i4,a1,f20.12,2x,f20.12,2x,f20.12)") ncore+i,'b',
     $tscaleb(i,2),tscaleb(i,1),(tscaleb(i,1)+tscaleb(i,2))/tscaleb(i,2)
           enddo
         endif
       endif
       do i=1,nalc
         if(dabs(tscalea(i,2)).gt.cctol) then
          tscalea(i,1)=(tscalea(i,1)+tscalea(i,2))/tscalea(i,2)
         else
          tscalea(i,1)=0.d0
         endif
       enddo
       if(scftype.ne.'rhf  ') then
         do i=1,nbec
           if(dabs(tscaleb(i,2)).gt.cctol) then
            tscaleb(i,1)=(tscaleb(i,1)+tscaleb(i,2))/tscaleb(i,2)
           else
            tscaleb(i,1)=0.d0
           endif
         enddo
       endif
C Print results
      tfact=(emp2+ef12)/emp2
      write(iout,*)
      write(iout,"(' Reference energy [au]:           'f24.12)") eref
      write(iout,"(' CABS singles correction [au]:    'f24.12)") ecabs
      write(iout,"(' Corrected reference energy [au]: 'f24.12)")
     $eref+ecabs
       if(scftype.eq.'rohf ')
     $write(iout,"(' MP2 singles contribution [au]:   'f24.12)") et1
      write(iout,"(' MP2 correlation energy [au]:     'f24.12)") emp2
      write(iout,"(' MP2 energy [au]:                 'f24.12)")
     $eref+emp2
      write(iout,"(' F12 contribution [au]:           'f24.12)") ef12
      write(iout,"(' MP2-F12 correlation energy [au]: 'f24.12)")
     $emp2+ef12
      write(iout,"(' MP2-F12 energy [au]:             'f24.12)")
     $eref+ecabs+emp2+ef12
C
      write(iout,*)
      open(ifcfile,status='unknown',file='iface',position='append')
      call prtenergc('MP2-F12',eref+ecabs+emp2+ef12,eref+ecabs,1)
      close(ifcfile)
      ef12=ef12-ecoup
      call dbldealloc(iscr)
C
      return
      end
C
************************************************************************
      subroutine epair_st(nal,nalc,nbfcn,nbfan,nvirtal,ncore,nbasis,
     $dfnbasis,fa,exa,eoa,eva,iout,verbosity,ef12,emp2,ecoup,vifaa,
     $cf12aa,epaa,rsij,f12ij,f12sy,scr,vintpq,vintao,vintoa,cf12aaf,
     $vifaaf,tscalea,cctol)
************************************************************************
* Calculate MP2-F12 pair energies for RHF
************************************************************************
      implicit none
      integer nbfcn,nbfan,nal,nalc,nvirtal,i,a,j,b,ncore,ij
      integer o,p,q,nbasis,dfnbasis,iout,verbosity
      real*8 eoa(*),eva(*),fa(nbfcn,nbfcn),eij,eijb,emp2,tmp,ef12
      real*8 exa(nbfcn,nbfcn),e12,ecoup,eco,epaa(*),tscalea(nalc,2)
      real*8 vifaa(nbfcn,dfnbasis,nalc),cf12aa(nbfcn,dfnbasis,nalc)
      real*8 rsij(nbfcn,nbfcn),f12ij(nbfcn,nbfcn),f12sy(nbfcn,nbfcn)
      real*8 scr(nbfcn,nbfcn),vintpq(nbasis,nbasis),vintao(nbfan,nal)
      real*8 vintoa(nal,nbfan),cf12aaf(nbfcn,dfnbasis,nalc),cctol
      real*8 vifaaf(nbfcn,dfnbasis,nalc)
C Transform fitting coefficients with Fock
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,cf12aa,nbfcn*dfnbasis,0d0,cf12aaf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,vifaa,nbfcn*dfnbasis,0d0,vifaaf,nbfcn*dfnbasis)
C Loop over pairs
      ij=0
       do i=1,nalc
         do j=1,i

C Call epair_st_inner (block version also calls this subroutine)

         write(iout,*)'ij,i,j: ',ij,i,j

          call epair_st_inner(ij,i,j,epaa,nbfcn,dfnbasis,f12ij,exa,
     $ f12sy,scr,nal,fa,nbfan,nbasis,nvirtal,rsij,vintpq,vintao,eoa,eva,
     $ ecoup,ef12,nalc,ncore,iout,verbosity,emp2,tscalea,vintoa,cctol,
     $ vifaa(1,1,i),vifaa(1,1,j),cf12aa(1,1,i),cf12aa(1,1,j),
     $ vifaaf(1,1,i),vifaaf(1,1,j),cf12aaf(1,1,i),cf12aaf(1,1,j))
c
C C Intermediate B1, B2; V, F12/r12 contribution; X, F12^2 contribution
C           ij=ij+1
C           e12=epaa(ij)
C C Intermediate B, term 3
C           call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0, vifaa(1,1,i),
C      $nbfcn,cf12aa(1,1,j),nbfcn,0d0,f12ij,nbfcn)
C           call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,cf12aa(1,1,i),
C      $nbfcn, vifaa(1,1,j),nbfcn,1d0,f12ij,nbfcn)
C           f12sy=transpose(f12ij)
C           call daxpy(nbfcn*nbfcn,7d0,f12ij,1,f12sy,1)
C           call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,exa,nbfcn,f12sy,
C      $nbfcn,0d0,scr,nbfcn)
C           call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,f12sy,nbfcn,exa,
C      $nbfcn,1d0,scr,nbfcn)
C C Intermediate B, term 4
C           call dgemm('n','n',nbfcn,nal,nbfcn,-1d0,fa,nbfcn,f12sy,
C      $nbfcn,1d0,scr,nbfcn)
C           call dgemm('n','n',nal,nbfcn,nbfcn,-1d0,f12sy,nbfcn,fa,
C      $nbfcn,1d0,scr,nbfcn)
C C Intermediate B, term 5
C           call dgemm('n','n',nal,nbfan,nal,1d0,fa,nbfcn,
C      $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
C           call dgemm('n','n',nbfan,nal,nal,1d0,f12sy(nbasis+1,1),
C      $nbfcn,fa,nbfcn,1d0,scr(nbasis+1,1),nbfcn)
C C Intermediate B, term 6
C           call dgemm('n','n',nbasis,nvirtal,nbasis,-1d0,fa,nbfcn,
C      $f12sy(1,nal+1),nbfcn,1d0,scr(1,nal+1),nbfcn)
C           call dgemm('n','n',nvirtal,nbasis,nbasis,-1d0,
C      $f12sy(nal+1,1),nbfcn,fa,nbfcn,1d0,scr(nal+1,1),nbfcn)
C C Intermediate B, term 7
C           call dgemm('n','n',nal,nbfan,nbfcn,-2d0,fa,nbfcn,
C      $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
C           call dgemm('n','n',nbfan,nal,nbfcn,-2d0,f12sy(nbasis+1,1),
C      $nbfcn,fa,nbfcn,1d0,scr(nbasis+1,1),nbfcn)
C C Intermediate B, term 8
C           rsij=0
C           call dgemm('n','n',nvirtal,nbasis,nbfan,1d0,f12ij(nal+1,
C      $nbasis+1),nbfcn,fa(nbasis+1,1),nbfcn,0d0,rsij(nal+1,1),nbfcn)
C           call dgemm('n','n',nbasis,nvirtal,nbfan,1d0,fa(1,nbasis+1),
C      $nbfcn,f12ij(nbasis+1,nal+1),nbfcn,1d0,rsij(1,nal+1),nbfcn)
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a)
C            do a=nal+1,nbasis
C             scr(nal+1:nbasis,a)=scr(nal+1:nbasis,a)-
C      $14d0*rsij(nal+1:nbasis,a)-2d0*rsij(a,nal+1:nbasis)
C            enddo
C C$OMP END PARALLEL DO
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(o)
C            do o=1,nal
C             scr(nal+1:nbasis,o)=scr(nal+1:nbasis,o)-
C      $14d0*rsij(nal+1:nbasis,o)-2d0*rsij(o,nal+1:nbasis)
C            enddo
C C$OMP END PARALLEL DO
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o)
C            do a=nal+1,nbasis
C           scr(1:nal,a)=scr(1:nal,a)-14d0*rsij(1:nal,a)-2d0*rsij(a,1:nal)
C            enddo
C C$OMP END PARALLEL DO
C           e12=e12+ddot(nbfcn*nbfcn,f12ij,1,scr,1)
C C Intermediate X, F12 contribution
C           call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0, vifaaf(1,1,i),
C      $nbfcn,cf12aa (1,1,j),nbfcn,0d0,vintpq,nbasis)
C           call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,cf12aaf(1,1,i),
C      $nbfcn, vifaa (1,1,j),nbfcn,1d0,vintpq,nbasis)
C           call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0, vifaa (1,1,i),
C      $nbfcn,cf12aaf(1,1,j),nbfcn,1d0,vintpq,nbasis)
C           call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,cf12aa (1,1,i),
C      $nbfcn, vifaaf(1,1,j),nbfcn,1d0,vintpq,nbasis)
C           call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
C      $ vifaaf(nbasis+1,1,i),nbfcn,cf12aa (1,1,j),nbfcn,0d0,vintao,nbfan)
C           call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
C      $cf12aaf(nbasis+1,1,i),nbfcn, vifaa (1,1,j),nbfcn,1d0,vintao,nbfan)
C           call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
C      $ vifaa (nbasis+1,1,i),nbfcn,cf12aaf(1,1,j),nbfcn,1d0,vintao,nbfan)
C           call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
C      $cf12aa (nbasis+1,1,i),nbfcn, vifaaf(1,1,j),nbfcn,1d0,vintao,nbfan)
C           call dgemm('n','t',nal,nbfan,dfnbasis,1d0, vifaaf(1,1,i),
C      $nbfcn,cf12aa (nbasis+1,1,j),nbfcn,0d0,vintoa,nal)
C           call dgemm('n','t',nal,nbfan,dfnbasis,1d0,cf12aaf(1,1,i),
C      $nbfcn, vifaa (nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
C           call dgemm('n','t',nal,nbfan,dfnbasis,1d0, vifaa (1,1,i),
C      $nbfcn,cf12aaf(nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
C           call dgemm('n','t',nal,nbfan,dfnbasis,1d0,cf12aa (1,1,i),
C      $nbfcn, vifaaf(nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
C           tmp=0.d0
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q) REDUCTION(+:tmp)
C            do p=1,nbasis
C              do q=1,nbasis
C               tmp=tmp+f12sy(p,q)*vintpq(p,q)
C              enddo
C            enddo
C C$OMP END PARALLEL DO
C           e12=e12+tmp
C           tmp=0.d0
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
C            do o=1,nal
C              do a=nbasis+1,nbfcn
C               tmp=tmp+f12sy(o,a)*vintoa(o,a-nbasis)
C      $               +f12sy(a,o)*vintao(a-nbasis,o)
C              enddo
C            enddo
C C$OMP END PARALLEL DO
C           e12=e12+tmp
C C Intermediate V, F12 contribution
C           call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,vifaa(1,1,i),
C      $nbfcn,vifaa(1,1,j),nbfcn,0d0,vintpq,nbasis)
C           call dgemm('n','t',nbfan,nal,dfnbasis,1d0,vifaa(nbasis+1,1,i),
C      $nbfcn,vifaa(1,1,j),nbfcn,0d0,vintao,nbfan)
C           call dgemm('n','t',nal,nbfan,dfnbasis,1d0,vifaa(1,1,i),
C      $nbfcn,vifaa(nbasis+1,1,j),nbfcn,0d0,vintoa,nal)
C           tmp=0.d0
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q) REDUCTION(+:tmp)
C            do p=1,nbasis
C              do q=1,nbasis
C               tmp=tmp-(5d0*f12ij(p,q)-f12ij(q,p))*vintpq(p,q)
C              enddo
C            enddo
C C$OMP END PARALLEL DO
C           e12=e12+8d0*tmp
C           tmp=0.d0
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
C            do o=1,nal
C              do a=nbasis+1,nbfcn
C               tmp=tmp-(5d0*f12ij(a,o)-f12ij(o,a))*vintao(a-nbasis,o)
C      $               -(5d0*f12ij(o,a)-f12ij(a,o))*vintoa(o,a-nbasis)
C              enddo
C            enddo
C C$OMP END PARALLEL DO
C           e12=e12+8d0*tmp
C C Intermediate C
C           eij=eoa(ncore+i)+eoa(ncore+j)
C           eco=0.d0
C C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,b,eijb) REDUCTION(+:eco)
C            do b=nal+1,nbasis
C             eijb=eij-eva(b-nal)
C              do a=nal+1,nbasis
C               eco=eco+(40d0*vintpq(a,b)-8d0*vintpq(b,a)
C      $+7d0*rsij(a,b)+rsij(b,a))*rsij(a,b)/(eijb-eva(a-nal))
C              enddo
C            enddo
C C$OMP END PARALLEL DO
C           e12=e12+eco
C           eco=eco/16d0
C           e12=e12/16d0
C C
C            if(i.eq.j) then
C             eco=0.5d0*eco
C             e12=0.5d0*e12
C            endif
C           ecoup=ecoup+eco
C           ef12=ef12+e12
C C Calculate MP2 energy
C           tmp=0.d0
C C$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(a,b,eijb)
C C$OMP& REDUCTION(+:tmp)
C            do b=nal+1,nal+nvirtal
C             eijb=eij-eva(b-nal)
C             tmp=tmp+0.5d0*vintpq(b,b)**2/(eijb-eva(b-nal))
C              do a=nal+1,b-1
C               tmp=tmp+(
C      $(vintpq(a,b)-vintpq(b,a))**2+vintpq(a,b)*vintpq(b,a))/
C      $(eijb-eva(a-nal))
C              enddo
C            enddo
C C$OMP END PARALLEL DO
C           tmp=4d0*tmp
C            if(i.eq.j) tmp=0.5d0*tmp
C           emp2=emp2+tmp
C            if(verbosity.gt.2) write(iout,"(2(i4,a1),f20.12,2x,f20.12)")
C      $ncore+i,' ',ncore+j,' ',tmp,e12
C            if(e12.ge.cctol) write(iout,*)
C      $'Warning! Non-negative pair energy!'
C           tscalea(i,1)=tscalea(i,1)+0.5d0*e12
C           tscalea(i,2)=tscalea(i,2)+0.5d0*tmp
C           tscalea(j,1)=tscalea(j,1)+0.5d0*e12
C           tscalea(j,2)=tscalea(j,2)+0.5d0*tmp
C           epaa(ij)=0.d0
C            if(dabs(tmp).gt.cctol) epaa(ij)=(tmp+e12)/tmp
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine mfocktodisk(nblocks,nbfcn,dfnbasis,fa,ncore,bsize,iv,
     $ vifaab,vifaafb,iv2)
************************************************************************
* Multiply by F and write each block to disk
************************************************************************
      implicit none
      integer i,j,nbfcn,dfnbasis,ncore,k,sblock1,sblock2
      integer nblocks ! how many blocks
      integer bsize ! how many elements in a block
      integer iv,iv2 ! I/O file descriptors
      real*8 fa(nbfcn,nbfcn) ! Fock matrix
      real*8 vifaab(nbfcn,dfnbasis,bsize) ! Source matrix
      real*8 vifaafb(nbfcn,dfnbasis,bsize) ! Result
      do i=1,nblocks
       sblock1 = (i-1)*bsize ! vifaafb start block
       do j=1,nblocks
        sblock2 = (j-1)*bsize ! vifaab start block
        do k=sblock2 + 1,sblock2+bsize
         read(iv2,rec=k)vifaab(1:nbfcn,1:dfnbasis,k-sblock2)
        enddo
        if(j.eq.1)then
         call dgemm('n','n',nbfcn*dfnbasis,bsize,bsize,1d0,vifaab,
     $    nbfcn*dfnbasis,fa(ncore+sblock2+1,ncore+sblock1+1),
     $    nbfcn,0d0,vifaafb,nbfcn*dfnbasis)
        else
         call dgemm('n','n',nbfcn*dfnbasis,bsize,bsize,1d0,vifaab,
     $    nbfcn*dfnbasis,fa(ncore+sblock2+1,ncore+sblock1+1),
     $    nbfcn,1d0,vifaafb,nbfcn*dfnbasis)
        endif
       enddo
C Calculating blocks of vifaaf, write every block to the disk
       do j=1,bsize
        write(iv,rec=j+sblock1)vifaafb(1:nbfcn,1:dfnbasis,j)
       enddo
      enddo
      return
      end
C
************************************************************************
      subroutine epair_st_inner(ij,i,j,epaa,nbfcn,dfnbasis,f12ij,
     $ exa,f12sy,scr,nal,fa,nbfan,nbasis,nvirtal,rsij,vintpq,vintao,eoa,
     $ eva,ecoup,ef12,nalc,ncore,iout,verbosity,emp2,tscalea,vintoa,
     $ cctol,vifaai,vifaaj,cf12aai,cf12aaj,
     $ vifaafi,vifaafj,cf12aafi,cf12aafj)
************************************************************************
* Inner loop of epair_st
************************************************************************
      implicit none
C
      integer nbfcn,nbfan,nal,nalc,nvirtal,i,a,j,b,ncore,ij
      integer o,p,q,nbasis,dfnbasis,iout,verbosity
      real*8 eoa(*),eva(*),fa(nbfcn,nbfcn),ddot,eij,eijb,emp2,tmp,ef12
      real*8 exa(nbfcn,nbfcn),e12,ecoup,eco,epaa(*),tscalea(nalc,2)
      real*8 rsij(nbfcn,nbfcn),f12ij(nbfcn,nbfcn),f12sy(nbfcn,nbfcn)
      real*8 scr(nbfcn,nbfcn),vintpq(nbasis,nbasis),vintao(nbfan,nal)
      real*8 vintoa(nal,nbfan),cctol
C Separated indexes
      real*8 vifaai(nbfcn,dfnbasis)
      real*8 vifaaj(nbfcn,dfnbasis)

      real*8 cf12aai(nbfcn,dfnbasis)
      real*8 cf12aaj(nbfcn,dfnbasis)

      real*8 vifaafi(nbfcn,dfnbasis)
      real*8 vifaafj(nbfcn,dfnbasis)

      real*8 cf12aafi(nbfcn,dfnbasis)
      real*8 cf12aafj(nbfcn,dfnbasis)

C
C Intermediate B1, B2; V, F12/r12 contribution; X, F12^2 contribution
          ij=ij+1
          e12=epaa(ij)
C Intermediate B, term 3
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,vifaai(1,1),
     $nbfcn,cf12aaj(1,1),nbfcn,0d0,f12ij,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,cf12aai(1,1),
     $nbfcn, vifaaj(1,1),nbfcn,1d0,f12ij,nbfcn)
          f12sy=transpose(f12ij)
          call daxpy(nbfcn*nbfcn,7d0,f12ij,1,f12sy,1)
          call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,exa,nbfcn,f12sy,
     $nbfcn,0d0,scr,nbfcn)
          call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,f12sy,nbfcn,exa,
     $nbfcn,1d0,scr,nbfcn)
C Intermediate B, term 4
          call dgemm('n','n',nbfcn,nal,nbfcn,-1d0,fa,nbfcn,f12sy,
     $nbfcn,1d0,scr,nbfcn)
          call dgemm('n','n',nal,nbfcn,nbfcn,-1d0,f12sy,nbfcn,fa,
     $nbfcn,1d0,scr,nbfcn)
C Intermediate B, term 5
          call dgemm('n','n',nal,nbfan,nal,1d0,fa,nbfcn,
     $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
          call dgemm('n','n',nbfan,nal,nal,1d0,f12sy(nbasis+1,1),
     $nbfcn,fa,nbfcn,1d0,scr(nbasis+1,1),nbfcn)
C Intermediate B, term 6
          call dgemm('n','n',nbasis,nvirtal,nbasis,-1d0,fa,nbfcn,
     $f12sy(1,nal+1),nbfcn,1d0,scr(1,nal+1),nbfcn)
          call dgemm('n','n',nvirtal,nbasis,nbasis,-1d0,
     $f12sy(nal+1,1),nbfcn,fa,nbfcn,1d0,scr(nal+1,1),nbfcn)
C Intermediate B, term 7
          call dgemm('n','n',nal,nbfan,nbfcn,-2d0,fa,nbfcn,
     $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
          call dgemm('n','n',nbfan,nal,nbfcn,-2d0,f12sy(nbasis+1,1),
     $nbfcn,fa,nbfcn,1d0,scr(nbasis+1,1),nbfcn)
C Intermediate B, term 8
          rsij=0
          call dgemm('n','n',nvirtal,nbasis,nbfan,1d0,f12ij(nal+1,
     $nbasis+1),nbfcn,fa(nbasis+1,1),nbfcn,0d0,rsij(nal+1,1),nbfcn)
          call dgemm('n','n',nbasis,nvirtal,nbfan,1d0,fa(1,nbasis+1),
     $nbfcn,f12ij(nbasis+1,nal+1),nbfcn,1d0,rsij(1,nal+1),nbfcn)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a)
           do a=nal+1,nbasis
            scr(nal+1:nbasis,a)=scr(nal+1:nbasis,a)-
     $14d0*rsij(nal+1:nbasis,a)-2d0*rsij(a,nal+1:nbasis)
           enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(o)
           do o=1,nal
            scr(nal+1:nbasis,o)=scr(nal+1:nbasis,o)-
     $14d0*rsij(nal+1:nbasis,o)-2d0*rsij(o,nal+1:nbasis)
           enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o)
           do a=nal+1,nbasis
          scr(1:nal,a)=scr(1:nal,a)-14d0*rsij(1:nal,a)-2d0*rsij(a,1:nal)
           enddo
C$OMP END PARALLEL DO
          e12=e12+ddot(nbfcn*nbfcn,f12ij,1,scr,1)
C Intermediate X, F12 contribution
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0, vifaafi(1,1),
     $nbfcn,cf12aaj(1,1),nbfcn,0d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,cf12aafi(1,1),
     $nbfcn, vifaaj(1,1),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0, vifaai(1,1),
     $nbfcn,cf12aafj(1,1),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,cf12aai(1,1),
     $nbfcn, vifaafj(1,1),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $ vifaafi(nbasis+1,1),nbfcn,cf12aaj(1,1),nbfcn,0d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $cf12aafi(nbasis+1,1),nbfcn, vifaaj(1,1),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $ vifaai(nbasis+1,1),nbfcn,cf12aafj(1,1),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $cf12aai(nbasis+1,1),nbfcn, vifaafj(1,1),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0, vifaafi(1,1),
     $nbfcn,cf12aaj(nbasis+1,1),nbfcn,0d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0,cf12aafi(1,1),
     $nbfcn, vifaaj(nbasis+1,1),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0, vifaai(1,1),
     $nbfcn,cf12aafj(nbasis+1,1),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0,cf12aai(1,1),
     $nbfcn, vifaafj(nbasis+1,1),nbfcn,1d0,vintoa,nal)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q) REDUCTION(+:tmp)
           do p=1,nbasis
             do q=1,nbasis
              tmp=tmp+f12sy(p,q)*vintpq(p,q)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do o=1,nal
             do a=nbasis+1,nbfcn
              tmp=tmp+f12sy(o,a)*vintoa(o,a-nbasis)
     $               +f12sy(a,o)*vintao(a-nbasis,o)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
C Intermediate V, F12 contribution
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,vifaai(1,1),
     $nbfcn,vifaaj(1,1),nbfcn,0d0,vintpq,nbasis)
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,vifaai(nbasis+1,1),
     $nbfcn,vifaaj(1,1),nbfcn,0d0,vintao,nbfan)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0,vifaai(1,1),
     $nbfcn,vifaaj(nbasis+1,1),nbfcn,0d0,vintoa,nal)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q) REDUCTION(+:tmp)
           do p=1,nbasis
             do q=1,nbasis
              tmp=tmp-(5d0*f12ij(p,q)-f12ij(q,p))*vintpq(p,q)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+8d0*tmp
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do o=1,nal
             do a=nbasis+1,nbfcn
              tmp=tmp-(5d0*f12ij(a,o)-f12ij(o,a))*vintao(a-nbasis,o)
     $               -(5d0*f12ij(o,a)-f12ij(a,o))*vintoa(o,a-nbasis)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+8d0*tmp
C Intermediate C
          eij=eoa(ncore+i)+eoa(ncore+j)
          eco=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,b,eijb) REDUCTION(+:eco)
           do b=nal+1,nbasis
            eijb=eij-eva(b-nal)
             do a=nal+1,nbasis
              eco=eco+(40d0*vintpq(a,b)-8d0*vintpq(b,a)
     $+7d0*rsij(a,b)+rsij(b,a))*rsij(a,b)/(eijb-eva(a-nal))
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+eco
          eco=eco/16d0
          e12=e12/16d0
C
           if(i.eq.j) then
            eco=0.5d0*eco
            e12=0.5d0*e12
           endif
          ecoup=ecoup+eco
          ef12=ef12+e12
C Calculate MP2 energy
          tmp=0.d0
C$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(a,b,eijb)
C$OMP& REDUCTION(+:tmp)
           do b=nal+1,nal+nvirtal
            eijb=eij-eva(b-nal)
            tmp=tmp+0.5d0*vintpq(b,b)**2/(eijb-eva(b-nal))
             do a=nal+1,b-1
              tmp=tmp+(
     $(vintpq(a,b)-vintpq(b,a))**2+vintpq(a,b)*vintpq(b,a))/
     $(eijb-eva(a-nal))
             enddo
           enddo
C$OMP END PARALLEL DO
          tmp=4d0*tmp
           if(i.eq.j) tmp=0.5d0*tmp
          emp2=emp2+tmp
           if(verbosity.gt.2) write(iout,"(2(i4,a1),f20.12,2x,f20.12)")
     $ncore+i,' ',ncore+j,' ',tmp,e12
           if(e12.ge.cctol) write(iout,*)
     $'Warning! Non-negative pair energy!'
          tscalea(i,1)=tscalea(i,1)+0.5d0*e12
          tscalea(i,2)=tscalea(i,2)+0.5d0*tmp
          tscalea(j,1)=tscalea(j,1)+0.5d0*e12
          tscalea(j,2)=tscalea(j,2)+0.5d0*tmp
          epaa(ij)=0.d0
          if(dabs(tmp).gt.cctol) epaa(ij)=(tmp+e12)/tmp
C
      return
      end
C
************************************************************************
      subroutine epair_st_bl(nal,nalc,nbfcn,nbfan,nvirtal,ncore,nbasis,
     $dfnbasis,fa,exa,eoa,eva,iout,verbosity,ef12,emp2,ecoup,vifaa,
     $cf12aa,epaa,rsij,f12ij,f12sy,scr,vintpq,vintao,vintoa,cf12aaf,
     $vifaaf,tscalea,cctol,bsize,ifltln,nblocks)
************************************************************************
* Calculate MP2-F12 pair energies for RHF 
* BLOCK version
************************************************************************
      implicit none
      integer nbfcn,nbfan,nal,nalc,nvirtal,i,a,j,b,ncore,ij, ii, jj
      integer o,p,q,nbasis,dfnbasis,iout,verbosity,nalx,k,ivifaa, kk
      integer ii1,ii2,jj1,jjj1,ii3,ii4,jj2,bsize,bid,ivifaad,ifltln
      integer iblock,nblocks,jindex,kstart,mind
      real*8 eoa(*),eva(*),fa(nbfcn,nbfcn),eij,eijb,emp2,tmp,ef12
      real*8 exa(nbfcn,nbfcn),e12,ecoup,eco,epaa(*),tscalea(nalc,2)
      real*8 rsij(nbfcn,nbfcn),f12ij(nbfcn,nbfcn),f12sy(nbfcn,nbfcn)
      real*8 scr(nbfcn,nbfcn),vintpq(nbasis,nbasis),vintao(nbfan,nal)
      real*8 vintoa(nal,nbfan),cctol,multtmp
C Memory heavy matrices
      real*8 vifaa(nbfcn,dfnbasis,bsize+1)

      real*8 vifaaf(nbfcn,dfnbasis,bsize+1)

      real*8 cf12aa(nbfcn,dfnbasis,bsize+1)

      real*8 cf12aaf(nbfcn,dfnbasis,bsize+1)

C Transform fitting coefficients with Fock
c     call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
c    $nbfcn,cf12aa,nbfcn*dfnbasis,0d0,cf12aaf,nbfcn*dfnbasis)
c     call dsymm('r','l',nbfcn*dfnbasis,nalc/2,1d0,fa(ncore+1,ncore+1),
c    $nbfcn,vifaa,nbfcn*dfnbasis,0d0,vifaaf,nbfcn*dfnbasis)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(104,status='unknown',file='MPVIFAA',form='unformatted')
      read(104)nalx
      close(104)

      write(iout,*)'nalx ',nalx
      write(iout,*)'nalc ',nalc

      open(101,file='MPVIFAA2',access='direct',
     $recl=dfnbasis*nbfcn*ifltln,status='old')
c     do k=1,nalx 
c      read(101,rec=k)vifaa(1:nbfcn,1:dfnbasis,k)
c     enddo
      open(110,file='MPCF12AA',access='direct',
     $recl=dfnbasis*nbfcn*ifltln,status='old')
c     do k=1,nalc
c      read(110,rec=k)cf12aa(1:nbfcn,1:dfnbasis,k)
c     enddo

      open(102,status='unknown',file='MPVIFAAB',access='direct',
     $ recl=dfnbasis*nbfcn*ifltln)

      call mfocktodisk(nblocks,nbfcn,dfnbasis,fa,ncore,bsize,102,vifaa,
     $ vifaaf, 101)

      open(103,status='unknown',file='MPCF12AAB',access='direct',
     $ recl=dfnbasis*nbfcn*ifltln)

      call mfocktodisk(nblocks,nbfcn,dfnbasis,fa,ncore,bsize,103,cf12aa
     $ ,vifaaf, 110)

C Loop over pairs
      ij=0
      do iblock=1,nblocks ! loop over blocks

       kstart = (iblock-1)*bsize + 1
       do k=kstart,kstart+bsize-1
        mind = k-(iblock-1)*bsize
        read(101,rec=k)vifaa(1:nbfcn,1:dfnbasis,mind)
        read(102,rec=k)vifaaf(1:nbfcn,1:dfnbasis,mind)
        read(110,rec=k)cf12aa(1:nbfcn,1:dfnbasis,mind)
        read(103,rec=k)cf12aaf(1:nbfcn,1:dfnbasis,mind)
       enddo
C
       do i=1,bsize ! relative i index inside a block
        ii = i + (iblock - 1)*bsize ! absolute i index
        do j=1,ii
C Check whether j is in memory
         jj = j - (iblock - 1)*bsize
         if(j.lt.kstart)then ! j is not in memory if it is lower than the start of the block 
          jj = bsize+1
          read(101,rec=j)vifaa(1:nbfcn,1:dfnbasis,jj)
          read(102,rec=j)vifaaf(1:nbfcn,1:dfnbasis,jj)
          read(110,rec=j)cf12aa(1:nbfcn,1:dfnbasis,jj)
          read(103,rec=j)cf12aaf(1:nbfcn,1:dfnbasis,jj)
         endif
C Call epair_st_inner
         write(iout,*)'ij,ii,j: ',ij,ii,j
         call epair_st_inner(ij,ii,j,epaa,nbfcn,dfnbasis,f12ij,exa,
     $ f12sy,scr,nal,fa,nbfan,nbasis,nvirtal,rsij,vintpq,vintao,eoa,eva,
     $ ecoup,ef12,nalc,ncore,iout,verbosity,emp2,tscalea,vintoa,cctol,
     $ vifaa(1,1,i),vifaa(1,1,jj),cf12aa(1,1,i),cf12aa(1,1,jj),
     $ vifaaf(1,1,i),vifaaf(1,1,jj),cf12aaf(1,1,i),cf12aaf(1,1,jj))
         enddo
        enddo
       enddo
C
C Close block files
      close(101)
      close(102)
      close(103)
      close(110)
      return
C
      end
C
************************************************************************
      subroutine epair_os(nal,nbe,nalc,nbec,nbfcn,nbfan,nvirtal,nvirtbe,
     $ncore,nbasis,dfnbasis,fa,fb,exa,exb,eoa,eva,eob,evb,iout,
     $verbosity,ef12,emp2,ecoup,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,
     $cf12ab,cf12ba,epab,rsa,rsb,f12ij,f12sy,scr,vintpq,vintao,vintoa,
     $cf12aaf,cf12bbf,cf12abf,cf12baf,vifaaf,vifbbf,vifabf,vifbaf,
     $tscalea,tscaleb,cctol)
************************************************************************
* Calculate MP2-F12 opposite-spin pair energies
************************************************************************
      implicit none
      integer nbfcn,nbfan,nal,nbe,nalc,nbec,nvirtal,nvirtbe,ncore,k,l
      integer o,s,p,q,u,w,nbasis,c,m,n,iout,verbosity,i,a,j,b,dfnbasis
      real*8 eoa(*),eva(*),eob(*),evb(*),fa(nbfcn,nbfcn),fb(nbfcn,nbfcn)
      real*8 exa(nbfcn,nbfcn),exb(nbfcn,nbfcn),emp2,tmp,ef12,eij,eijb
      real*8 rsa(nbfcn,nbfcn),rsb(nbfcn,nbfcn),epab(nalc,nbec),cctol
      real*8 vifaa(nbfcn,dfnbasis,nalc),cf12aa(nbfcn,dfnbasis,nalc)
      real*8 vifbb(nbfcn,dfnbasis,nbec),cf12bb(nbfcn,dfnbasis,nbec)
      real*8 vifab(nbfcn,dfnbasis,nalc),cf12ab(nbfcn,dfnbasis,nalc)
      real*8 vifba(nbfcn,dfnbasis,nbec),cf12ba(nbfcn,dfnbasis,nbec)
      real*8 f12ij(nbfcn,nbfcn),f12sy(nbfcn,nbfcn),e12,ecoup,eco,ddot
      real*8 scr(nbfcn,nbfcn),vintpq(nbasis,nbasis),vintao(nbfan,nbe)
      real*8 vintoa(nal,nbfan),tscalea(nalc,2),tscaleb(nbec,2)
      real*8 vifaaf(nbfcn,dfnbasis,nalc),cf12aaf(nbfcn,dfnbasis,nalc)
      real*8 vifbbf(nbfcn,dfnbasis,nbec),cf12bbf(nbfcn,dfnbasis,nbec)
      real*8 vifabf(nbfcn,dfnbasis,nalc),cf12abf(nbfcn,dfnbasis,nalc)
      real*8 vifbaf(nbfcn,dfnbasis,nbec),cf12baf(nbfcn,dfnbasis,nbec)
C Transform fitting coefficients with Fock
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,cf12aa,nbfcn*dfnbasis,0d0,cf12aaf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nbec,1d0,fb(ncore+1,ncore+1),
     $nbfcn,cf12bb,nbfcn*dfnbasis,0d0,cf12bbf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,cf12ab,nbfcn*dfnbasis,0d0,cf12abf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nbec,1d0,fb(ncore+1,ncore+1),
     $nbfcn,cf12ba,nbfcn*dfnbasis,0d0,cf12baf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,vifaa,nbfcn*dfnbasis,0d0,vifaaf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nbec,1d0,fb(ncore+1,ncore+1),
     $nbfcn,vifbb,nbfcn*dfnbasis,0d0,vifbbf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,vifab,nbfcn*dfnbasis,0d0,vifabf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nbec,1d0,fb(ncore+1,ncore+1),
     $nbfcn,vifba,nbfcn*dfnbasis,0d0,vifbaf,nbfcn*dfnbasis)
C
       do i=1,nalc
         do j=1,nbec
C Intermediate B1, B2; V, F12/r12 contribution; X, F12^2 contribution
          e12=epab(i,j)
C Intermediate B, term 3
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,3d0/8d0, vifaa(1,1,i),
     $nbfcn,cf12bb(1,1,j),nbfcn,0d0,f12sy,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,3d0/8d0,cf12aa(1,1,i),
     $nbfcn, vifbb(1,1,j),nbfcn,1d0,f12sy,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0/8d0, vifba(1,1,j),
     $nbfcn,cf12ab(1,1,i),nbfcn,1d0,f12sy,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0/8d0,cf12ba(1,1,j),
     $nbfcn, vifab(1,1,i),nbfcn,1d0,f12sy,nbfcn)
          call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,exa,nbfcn,f12sy,
     $nbfcn,0d0,scr,nbfcn)
          call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,f12sy,nbfcn,exb,
     $nbfcn,1d0,scr,nbfcn)
C Intermediate B, term 4
          call dgemm('n','n',nbfcn,nbe,nbfcn,-1d0,fa,nbfcn,f12sy,
     $nbfcn,1d0,scr,nbfcn)
          call dgemm('n','n',nal,nbfcn,nbfcn,-1d0,f12sy,nbfcn,fb,
     $nbfcn,1d0,scr,nbfcn)
C Intermediate B, term 5
          call dgemm('n','n',nal,nbfan,nal,1d0,fa,nbfcn,
     $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
          call dgemm('n','n',nbfan,nbe,nbe,1d0,f12sy(nbasis+1,1),
     $nbfcn,fb,nbfcn,1d0,scr(nbasis+1,1),nbfcn)
C Intermediate B, term 6
          call dgemm('n','n',nbasis,nvirtbe,nbasis,-1d0,fa,nbfcn,
     $f12sy(1,nbe+1),nbfcn,1d0,scr(1,nbe+1),nbfcn)
          call dgemm('n','n',nvirtal,nbasis,nbasis,-1d0,
     $f12sy(nal+1,1),nbfcn,fb,nbfcn,1d0,scr(nal+1,1),nbfcn)
C Intermediate B, term 7
          call dgemm('n','n',nal,nbfan,nbfcn,-2d0,fa,nbfcn,
     $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
          call dgemm('n','n',nbfan,nbe,nbfcn,-2d0,f12sy(nbasis+1,1),
     $nbfcn,fb,nbfcn,1d0,scr(nbasis+1,1),nbfcn)
C Intermediate B, term 8
          call dgemm('n','n',nvirtal,nbasis,nbfan,1d0,f12sy(nal+1,
     $nbasis+1),nbfcn,fb(nbasis+1,1),nbfcn,0d0,rsa(nal+1,1),nbfcn)
          call dgemm('n','n',nbasis,nvirtbe,nbfan,1d0,fa(1,nbasis+1),
     $nbfcn,f12sy(nbasis+1,nbe+1),nbfcn,0d0,rsb(1,nbe+1),nbfcn)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(q)
           do q=1,nbasis
            scr(nal+1:nal+nbasis,q)=
     $      scr(nal+1:nal+nbasis,q)-2d0*rsa(nal+1:nal+nbasis,q)
           enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(b)
           do b=nbe+1,nbasis
            scr(1:nbasis,b)=scr(1:nbasis,b)-2d0*rsb(1:nbasis,b)
           enddo
C$OMP END PARALLEL DO
          e12=e12+ddot(nbfcn*nbfcn,f12sy,1,scr,1)
C Intermediate X, F12 contribution
          call dgemm('n','t',nbasis,nbasis,dfnbasis,3d0/8d0,
     $ vifaaf(1,1,i),nbfcn,cf12bb (1,1,j),nbfcn,0d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,3d0/8d0,
     $cf12aaf(1,1,i),nbfcn, vifbb (1,1,j),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0/8d0,
     $ vifba (1,1,j),nbfcn,cf12abf(1,1,i),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0/8d0,
     $cf12ba (1,1,j),nbfcn, vifabf(1,1,i),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,3d0/8d0,
     $ vifaa (1,1,i),nbfcn,cf12bbf(1,1,j),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,3d0/8d0,
     $cf12aa (1,1,i),nbfcn, vifbbf(1,1,j),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0/8d0,
     $ vifbaf(1,1,j),nbfcn,cf12ab (1,1,i),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0/8d0,
     $cf12baf(1,1,j),nbfcn, vifab (1,1,i),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbfan,nbe,dfnbasis,3d0/8d0,
     $ vifaaf(nbasis+1,1,i),nbfcn,cf12bb (1,1,j),nbfcn,0d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,3d0/8d0,
     $cf12aaf(nbasis+1,1,i),nbfcn, vifbb (1,1,j),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,1d0/8d0,
     $ vifba (nbasis+1,1,j),nbfcn,cf12abf(1,1,i),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,1d0/8d0,
     $cf12ba (nbasis+1,1,j),nbfcn, vifabf(1,1,i),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,3d0/8d0,
     $ vifaa (nbasis+1,1,i),nbfcn,cf12bbf(1,1,j),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,3d0/8d0,
     $cf12aa (nbasis+1,1,i),nbfcn, vifbbf(1,1,j),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,1d0/8d0,
     $ vifbaf(nbasis+1,1,j),nbfcn,cf12ab (1,1,i),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nbe,dfnbasis,1d0/8d0,
     $cf12baf(nbasis+1,1,j),nbfcn, vifab (1,1,i),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nal,nbfan,dfnbasis,3d0/8d0,
     $ vifaaf(1,1,i),nbfcn,cf12bb (nbasis+1,1,j),nbfcn,0d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,3d0/8d0,
     $cf12aaf(1,1,i),nbfcn, vifbb (nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0/8d0,
     $ vifba (1,1,j),nbfcn,cf12abf(nbasis+1,1,i),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0/8d0,
     $cf12ba (1,1,j),nbfcn, vifabf(nbasis+1,1,i),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,3d0/8d0,
     $ vifaa (1,1,i),nbfcn,cf12bbf(nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,3d0/8d0,
     $cf12aa (1,1,i),nbfcn, vifbbf(nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0/8d0,
     $ vifbaf(1,1,j),nbfcn,cf12ab (nbasis+1,1,i),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0/8d0,
     $cf12baf(1,1,j),nbfcn, vifab (nbasis+1,1,i),nbfcn,1d0,vintoa,nal)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q) REDUCTION(+:tmp)
           do p=1,nbasis
             do q=1,nbasis
              tmp=tmp+f12sy(p,q)*vintpq(p,q)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do o=1,nbe
             do a=nbasis+1,nbfcn
              tmp=tmp+f12sy(a,o)*vintao(a-nbasis,o)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do a=nbasis+1,nbfcn
             do o=1,nal
              tmp=tmp+f12sy(o,a)*vintoa(o,a-nbasis)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
C Intermediate V, F12 contribution
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,vifaa(1,1,i),
     $nbfcn,vifbb(1,1,j),nbfcn,0d0,vintpq,nbasis)
          call dgemm('n','t',nbfan,nbe,dfnbasis,1d0,vifaa(nbasis+1,1,i),
     $nbfcn,vifbb(1,1,j),nbfcn,0d0,vintao,nbfan)
          call dgemm('n','t',nal,nbfan,dfnbasis,1d0,vifaa(1,1,i),
     $nbfcn,vifbb(nbasis+1,1,j),nbfcn,0d0,vintoa,nal)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q) REDUCTION(+:tmp)
           do p=1,nbasis
             do q=1,nbasis
              tmp=tmp-f12sy(p,q)*vintpq(p,q)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+2.d0*tmp
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do a=nbasis+1,nbfcn
             do o=1,nal
              tmp=tmp-f12sy(o,a)*vintoa(o,a-nbasis)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+2.d0*tmp
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do a=nbasis+1,nbfcn
             do o=1,nbe
              tmp=tmp-f12sy(a,o)*vintao(a-nbasis,o)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+2.d0*tmp
C Intermediate C
          eij=eoa(ncore+i)+eob(ncore+j)
          eco=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,b,eijb,tmp) REDUCTION(+:eco)
           do b=nbe+1,nbasis
            eijb=eij-evb(b-nbe)
             do a=nal+1,nbasis
              tmp=rsb(a,b)+rsa(a,b)
              eco=eco+(vintpq(a,b)+0.5d0*tmp)*tmp/(eijb-eva(a-nal))
             enddo
           enddo
C$OMP END PARALLEL DO
          eco=2d0*eco
          ecoup=ecoup+eco
          e12=e12+eco
          ef12=ef12+e12
C Calculate MP2 energy
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,b,eijb) REDUCTION(+:tmp)
           do b=nbe+1,nbasis
            eijb=eij-evb(b-nbe)
             do a=nal+1,nbasis
              tmp=tmp+vintpq(a,b)**2/(eijb-eva(a-nal))
             enddo
           enddo
C$OMP END PARALLEL DO
          emp2=emp2+tmp
           if(verbosity.gt.2) write(iout,"(2(i4,a1),f20.12,2x,f20.12)")
     $ncore+i,'a',ncore+j,'b',tmp,e12
           if(e12.ge.cctol) write(iout,*)
     $'Warning! Non-negative pair energy!'
          tscalea(i,1)=tscalea(i,1)+0.5d0*e12
          tscalea(i,2)=tscalea(i,2)+0.5d0*tmp
          tscaleb(j,1)=tscaleb(j,1)+0.5d0*e12
          tscaleb(j,2)=tscaleb(j,2)+0.5d0*tmp
          epab(i,j)=0.d0
           if(dabs(tmp).gt.cctol) epab(i,j)=(tmp+e12)/tmp
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine epair_ss(nal,nalc,nbfcn,nbfan,nvirtal,ncore,nbasis,
     $dfnbasis,fa,exa,eoa,eva,iout,verbosity,ef12,emp2,ecoup,vifaa,
     $cf12aa,epaa,rsij,f12ij,f12sy,scr,vintpq,vintao,vintoa,cf12aaf,
     $vifaaf,spl,tscalea,cctol)
************************************************************************
* Calculate MP2-F12 same-spin pair energies
************************************************************************
      implicit none
      integer nbfcn,nbfan,nal,nalc,nvirtal,i,a,j,b,ncore,dfnbasis,ij
      integer o,p,q,nbasis,iout,verbosity
      real*8 emp2,tmp,ef12,ddot,temp,eijb,ecoup,tscalea(nalc,2),cctol
      real*8 eoa(*),eva(*),fa(nbfcn,nbfcn),exa(nbfcn,nbfcn),eij,e12,eco
      real*8 vifaa(nbfcn,dfnbasis,nalc),cf12aa(nbfcn,dfnbasis,nalc)
      real*8 rsij(nbfcn,nbfcn),f12ij(nbfcn,nbfcn),f12sy(nbfcn,nbfcn)
      real*8 scr(nbfcn,nbfcn),vintpq(nbasis,nbasis),vintao(nbfan,nal)
      real*8 vintoa(nal,nbfan),cf12aaf(nbfcn,dfnbasis,nalc),epaa(*)
      real*8 vifaaf(nbfcn,dfnbasis,nalc)
      character(len=1) spl
C Transform fitting coefficients with Fock
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,cf12aa,nbfcn*dfnbasis,0d0,cf12aaf,nbfcn*dfnbasis)
      call dsymm('r','l',nbfcn*dfnbasis,nalc,1d0,fa(ncore+1,ncore+1),
     $nbfcn,vifaa,nbfcn*dfnbasis,0d0,vifaaf,nbfcn*dfnbasis)
C Loop over pairs
      ij=0
       do i=1,nalc
         do j=1,i-1
C Intermediate B1, B2; V, F12/r12 contribution; X, F12^2 contribution
          ij=ij+1
          e12=epaa(ij)
C Intermediate B, term 3
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0, vifaa(1,1,i),
     $nbfcn,cf12aa(1,1,j),nbfcn,0d0,f12ij,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,cf12aa(1,1,i),
     $nbfcn, vifaa(1,1,j),nbfcn,1d0,f12ij,nbfcn)
          f12sy=-transpose(f12ij)
          call daxpy(nbfcn*nbfcn,1d0,f12ij,1,f12sy,1)
          call dgemm('n','n',nbfcn,nbfcn,nbfcn,1d0,exa,nbfcn,f12sy,
     $nbfcn,0d0,scr,nbfcn)
C Intermediate B, term 4
          call dgemm('n','n',nbfcn,nal,nbfcn,-1d0,fa,nbfcn,f12sy,
     $nbfcn,1d0,scr,nbfcn)
C Intermediate B, term 5
          call dgemm('n','n',nal,nbfan,nal,1d0,fa,nbfcn,
     $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
C Intermediate B, term 6
          call dgemm('n','n',nbasis,nvirtal,nbasis,-1d0,fa,nbfcn,
     $f12sy(1,nal+1),nbfcn,1d0,scr(1,nal+1),nbfcn)
C Intermediate B, term 7
          call dgemm('n','n',nal,nbfan,nbfcn,-2d0,fa,nbfcn,
     $f12sy(1,nbasis+1),nbfcn,1d0,scr(1,nbasis+1),nbfcn)
C Intermediate B, term 8
          rsij=0.d0
          call dgemm('n','n',nbasis,nvirtal,nbfan,1d0,fa(1,nbasis+1),
     $nbfcn,f12sy(nbasis+1,nal+1),nbfcn,0d0,rsij(1,nal+1),nbfcn)
          call daxpy(nbfcn*nvirtal,-2d0,rsij(1,nal+1),1,scr(1,nal+1),1)
          e12=e12+ddot(nbfcn*nbfcn,f12sy,1,scr,1)
C Intermediate X, F12 contribution
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0, vifaaf(1,1,i),
     $nbfcn,cf12aa (1,1,j),nbfcn,0d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,cf12aaf(1,1,i),
     $nbfcn, vifaa (1,1,j),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0, vifaa (1,1,i),
     $nbfcn,cf12aaf(1,1,j),nbfcn,1d0,vintpq,nbasis)
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,cf12aa (1,1,i),
     $nbfcn, vifaaf(1,1,j),nbfcn,1d0,vintpq,nbasis)
          tmp=0.d0
C$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(p,q)
C$OMP& REDUCTION(+:tmp)
           do p=1,nbasis
             do q=1,p-1
              tmp=tmp+f12sy(p,q)*(vintpq(p,q)-vintpq(q,p))
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $ vifaa (nbasis+1,1,i),nbfcn,cf12aaf(1,1,j),nbfcn,0d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $cf12aa (nbasis+1,1,i),nbfcn, vifaaf(1,1,j),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,-1d0,
     $cf12aaf(nbasis+1,1,j),nbfcn, vifaa (1,1,i),nbfcn,1d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,-1d0,
     $ vifaaf(nbasis+1,1,j),nbfcn,cf12aa (1,1,i),nbfcn,1d0,vintao,nbfan)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do o=1,nal
             do a=nbasis+1,nbfcn
              tmp=tmp+f12sy(a,o)*vintao(a-nbasis,o)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
          call dgemm('n','t',nal,nbfan,dfnbasis, 1d0, vifaaf(1,1,i),
     $nbfcn,cf12aa (nbasis+1,1,j),nbfcn,0d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis, 1d0,cf12aaf(1,1,i),
     $nbfcn, vifaa (nbasis+1,1,j),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,-1d0,
     $cf12aa (1,1,j),nbfcn, vifaaf(nbasis+1,1,i),nbfcn,1d0,vintoa,nal)
          call dgemm('n','t',nal,nbfan,dfnbasis,-1d0,
     $ vifaa (1,1,j),nbfcn,cf12aaf(nbasis+1,1,i),nbfcn,1d0,vintoa,nal)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do o=1,nal
             do a=nbasis+1,nbfcn
              tmp=tmp+f12sy(o,a)*vintoa(o,a-nbasis)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+tmp
C Intermediate V, F12 contribution
          call dgemm('n','t',nbasis,nbasis,dfnbasis,1d0,vifaa(1,1,i),
     $nbfcn,vifaa(1,1,j),nbfcn,0d0,vintpq,nbasis)
          tmp=0.d0
C$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(p,q,temp)
C$OMP& REDUCTION(+:tmp)
           do q=1,nbasis
             do p=1,q-1
              temp=vintpq(p,q)-vintpq(q,p)
              vintpq(p,q)=temp
              tmp=tmp-f12sy(p,q)*temp
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+8d0*tmp
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $vifaa(nbasis+1,1,i),nbfcn,vifaa(1,1,j),nbfcn,0d0,vintao,nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,-1d0,
     $vifaa(nbasis+1,1,j),nbfcn,vifaa(1,1,i),nbfcn,1d0,vintao,nbfan)
          tmp=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,o) REDUCTION(+:tmp)
           do o=1,nal
             do a=nbasis+1,nbfcn
              tmp=tmp-f12sy(a,o)*vintao(a-nbasis,o)
             enddo
           enddo
C$OMP END PARALLEL DO
          e12=e12+8d0*tmp
C Intermediate C
          eij=eoa(ncore+i)+eoa(ncore+j)
          eco=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(a,b,eijb,tmp) REDUCTION(+:eco)
           do b=nal+1,nal+nvirtal
            eijb=eij-eva(b-nal)
             do a=nal+1,b-1
              tmp=0.25d0*(rsij(a,b)-rsij(b,a))
              eco=eco+(2d0*vintpq(a,b)+tmp)*tmp/(eijb-eva(a-nal))
             enddo
           enddo
C$OMP END PARALLEL DO
          ecoup=ecoup+eco
          e12=e12/16d0+eco
          ef12=ef12+e12
C Calculate MP2 energy
          tmp=0.d0
C$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(a,b,eijb)
C$OMP& REDUCTION(+:tmp)
           do b=nal+1,nal+nvirtal
            eijb=eij-eva(b-nal)
             do a=nal+1,b-1
              tmp=tmp+vintpq(a,b)**2/(eijb-eva(a-nal))
             enddo
           enddo
C$OMP END PARALLEL DO
          emp2=emp2+tmp
           if(verbosity.gt.2) write(iout,"(2(i4,a1),f20.12,2x,f20.12)")
     $ncore+j,spl,ncore+i,spl,tmp,e12
           if(e12.ge.cctol) write(iout,*)
     $'Warning! Non-negative pair energy!'
          tscalea(i,1)=tscalea(i,1)+0.5d0*e12
          tscalea(i,2)=tscalea(i,2)+0.5d0*tmp
          tscalea(j,1)=tscalea(j,1)+0.5d0*e12
          tscalea(j,2)=tscalea(j,2)+0.5d0*tmp
          epaa(ij)=0.d0
           if(dabs(tmp).gt.cctol) epaa(ij)=(tmp+e12)/tmp
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine mp2sing(nbfcn,ncore,nocc,nvirt,f,epso,epsv,et1,tscale)
************************************************************************
* Calculate MP2 singles contribution
************************************************************************
      implicit none
      integer nbfcn,ncore,nocc,nvirt,a,i
      real*8 f(nbfcn,nbfcn),epso(*),epsv(*),et1,et2,tmp
      real*8 tscale(nocc-ncore)
C
      et2=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(PRIVATE) SHARED(nocc,nvirt,ncore,f,epsv,epso,tscale)
C$OMP& REDUCTION(+:et2)
       do i=ncore+1,nocc
        tmp=0.d0
         do a=nocc+1,nocc+nvirt
          tmp=tmp+f(a,i)**2/(epso(i)-epsv(a-nocc))
         enddo
        et2=et2+tmp
        tscale(i-ncore)=tscale(i-ncore)+tmp
       enddo
C$OMP END PARALLEL DO
      et1=et1+et2
C
      return
      end
C
************************************************************************
      subroutine cabssing(nbfcn,nbfan,nocc,nvirt,ecabs,iout,f,epso,epsv,
     $diisfile,errfile,ifltln,tol,ft,t0,t1,verbosity)
************************************************************************
* Calculate CABS singles correction
************************************************************************
      implicit none
      integer nbfcn,nbfan,nocc,nvirt,nit,a,b,c,i,iout,maxit,diisfile
      integer errfile,ifltln,verbosity
      parameter(maxit=30)
      real*8 epso(*),epsv(*),f(nbfcn,nbfcn),dnrm2,tmp,ecabs,norm,tol
      real*8 bmat(maxit,maxit),invbmat(maxit+1,maxit+1),ddot,eold,enew
      real*8 ft(nbfan,nocc),t1(nbfan,nocc),t0(nbfan,nocc)
C Construct f tilde
C$OMP PARALLEL DO
C$OMP& DEFAULT(PRIVATE)
C$OMP& SHARED(nocc,nvirt,nbfcn,f,epsv,epso,ft,t0)
       do b=nocc+nvirt+1,nbfcn
         do i=1,nocc
          tmp=f(b,i)
           do a=nocc+1,nocc+nvirt
            tmp=tmp-f(b,a)*f(a,i)/(epsv(a-nocc)-epso(i))
           enddo
          ft(b-nocc-nvirt,i)=tmp
          t0(b-nocc-nvirt,i)=tmp/(epso(i)-epsv(b-nocc))
         enddo 
       enddo 
C$OMP END PARALLEL DO
C Calculate singles amplitudes
      nit=0
      norm=1d0
      enew=ddot(nbfan*nocc,ft,1,t0,1)
      eold=enew+1.d0
       if(verbosity.ge.3) write(iout,*) 'Iter.  Convergence'
       do while(nit.lt.maxit.and.
     $(norm.gt.tol.or.dabs(enew-eold).gt.tol))
        nit=nit+1
C$OMP PARALLEL DO
C$OMP& DEFAULT(PRIVATE)
C$OMP& SHARED(nocc,nvirt,nbfan,f,epsv,epso,ft,t0,t1)
         do b=1,nbfan
           do i=1,nocc
            tmp=-ft(b,i)
             do c=1,nbfan
               do a=nocc+1,nocc+nvirt
                tmp=tmp+f(b+nocc+nvirt,a)*f(c+nocc+nvirt,a)*t0(c,i)
     $/(epsv(a-nocc)-epso(i))
               enddo
             enddo
            t1(b,i)=tmp/(epsv(b+nvirt)-epso(i))
           enddo
         enddo
C$OMP END PARALLEL DO
        t0=t0-t1
        norm=dnrm2(nbfan*nocc,t0,1)
        call diis(nit,nbfan*nocc,t1,t0,maxit,diisfile,errfile,ifltln,
     $bmat,invbmat)
        t0=t1
        eold=enew
        enew=ddot(nbfan*nocc,ft,1,t0,1)
         if(verbosity.ge.3) write(iout,"(i4,f15.8)") nit,norm
       enddo
       if(norm.gt.tol) then
        write(iout,*) 'Warning! Convergence has not been achieved!'
        write(iout,"(' Residual:  ',f15.8)") norm
        write(iout,"(' Threshold: ',f15.8)") tol
       endif
      ecabs=ecabs+enew
C
      return
      end
C
************************************************************************
      subroutine transf(nbfc,nbfcn,fock,cmo,scr,ev,lll)
************************************************************************
* Transform one-electron quantity to CABS MO basis
************************************************************************
      implicit none
      integer nbfc,nbfcn,i
      real*8 fock(*),cmo(*),scr(*),ev(*)
      logical lll
C
      call dsymm('l','l',nbfc,nbfcn,1d0,fock,nbfc,cmo,nbfc,0d0,scr,nbfc)
      call dgemm('t','n',nbfcn,nbfcn,nbfc,1.d0,cmo,nbfc,scr,nbfc,
     $0.d0,fock,nbfcn)
       if(lll) then
         do i=1,nbfcn
          ev(i)=fock((i-1)*nbfcn+i)
         enddo
       endif
C
      return
      end
C
************************************************************************
      subroutine cancabs(nbfc,nbfan,nbasis,fock,cmo,scr,ev,um,iout)
************************************************************************
* Canonicalize CABS virtuals
************************************************************************
      implicit none
      integer nbfc,nbfan,nbasis,iisyev,iout
      real*8 fock(*),cmo(nbfc,nbfc),scr(*),ev(nbfan),um(nbfan,nbfan)
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
C
      call dsymm('l','l',nbfc,nbfan,1d0,fock,nbfc,cmo(1,nbasis+1),nbfc,
     $0d0,scr,nbfc)
      call dgemm('t','n',nbfan,nbfan,nbfc,1.d0,cmo(1,nbasis+1),nbfc,scr,
     $nbfc,0.d0,um,nbfan)
      call dsyev('V','U',nbfan,um,nbfan,ev,scr,3*nbfan**2,isyev)
       if(isyev.ne.0) then
        write(iout,*)
     $'Fatal error at the canonicalization of CABS virtuals!'
        call mrccend(1)
       endif
      call dgemm('n','n',nbfc,nbfan,nbfan,1.d0,cmo(1,nbasis+1),nbfc,um,
     $nbfan,0.d0,scr,nbfc)
      call dcopy(nbfc*nbfan,scr,1,cmo(1,nbasis+1),1)
C
      return
      end
C
************************************************************************
      subroutine moread(nbfc,nbasis,scrfile1,conv,cmo,mo,scr)
************************************************************************
* Read MOs and convert them to CABS basis MOs
************************************************************************
      implicit none
      integer nbfc,nbasis,conv(nbfc),scrfile1,i
      real*8 cmo(nbfc,nbfc),mo(nbasis,nbasis),scr(*)
      call readmo(scr,scr,mo,scrfile1,nbasis,nbasis)
       do i=1,nbfc
         if(conv(i).gt.0) cmo(i,1:nbasis)=mo(conv(i),1:nbasis)
       enddo
C
      return
      end
C
************************************************************************
      subroutine aoproj(nbfc,nbfa,nbasis,conv,cmoa,s12)
************************************************************************
* Project out AO basis from CABS
************************************************************************
      implicit none
      integer nbfc,nbfa,nbasis,i,j,conv(nbfc)
      real*8 cmoa(nbfc,nbfc),s12(nbasis,nbfa)
C
      cmoa=0.d0
       do i=1,nbfc
         if(conv(i).lt.0) then
          cmoa(i,nbasis+iabs(conv(i)))=1.d0
           do j=1,nbfc
             if(conv(j).gt.0) cmoa(j,nbasis+iabs(conv(i)))=
     $cmoa(j,nbasis+iabs(conv(i)))-s12(conv(j),iabs(conv(i)))
           enddo
         endif
       enddo
C
      return
      end
C
************************************************************************
      subroutine ccimed_bl(nb,dfnbasis,nal,nbe,nval,nvbe,ncore,focka,
     $fockb,
     $nbfcn,nbfan,emp2f12,ecabs,scrfile1,scftype,vv,vaa,vbb,vab,vpia,
     $vpib,caia,caib,uaia,uaib,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,
     $cf12ab,cf12ba,iout,dcore,imem,imem1,maxcor,f12ij,f12pq,f12ao,
     $f12ai,verbosity,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,tscalea,
     $tscaleb,tfact,epaa,epbb,epab,ecoup,ifltln,vifaa2,cf12aa2,nel)
************************************************************************
* Calculate intermediates for CCSD calculation 
* BLOCK version
************************************************************************
      implicit none
      integer nb,i,j,k,l,p,q,r,s,a,b,nal,nbe,ncore,o,rs,nbfcn,iout,nnb
      integer nval,nvbe,c,kl,ij,jk,ab,dfnbasis,scrfile1,pq,imem,dblalloc
      integer icp,icm,ivs,iva,imn,rsdims,slo,sup,imem1,nel2
      integer maxcor,maxmem,imn2
      integer verbosity,nbfan,rsdima,aodim,rsdim,ijs,ija
      integer icabij,icabijaa,icabijbb,icabijab,ifltln
      integer ib,nel,kstart,ii,jj,nbb,nbbl
      integer freememspace,iel

      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn),emp2f12,ecabs,tmp

      real*8 vifaa2(nbfcn,dfnbasis,nel+1),cf12aa2(nbfcn,dfnbasis,nel+1)
c     real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
C Memory heavy matrices
      real*8 vifaa(nbfcn,dfnbasis,nel+1)
      real*8 cf12aa(nbfcn,dfnbasis,nel+1),tfact
      real*8 vifbb(nbfcn,dfnbasis,nb),cf12bb(nbfcn,dfnbasis,nbe),ecoup
      real*8 vifab(nbfcn,dfnbasis,nb),cf12ab(nbfcn,dfnbasis,nal)
      real*8 vifba(nbfcn,dfnbasis,nb),cf12ba(nbfcn,dfnbasis,nbe)

      real*8 cf12r12aa(nb,dfnbasis,nal),cf12r12ab(nb,dfnbasis,nal)
      real*8 cf12r12bb(nb,dfnbasis,nbe),cf12r12ba(nb,dfnbasis,nbe)

      real*8 vaa((nb-1)*nb/2,(nal-1)*nal/2),vab(nb,nb,nal,nbe),dcore(*)
      real*8 vbb((nb-1)*nb/2,(nbe-1)*nbe/2),vv(nb,nb,(nal+1)*nal/2)
      real*8 vpia(nb,nal),vpib(nb,nbe),uaia(nval,nal),uaib(nvbe,nbe)
      real*8 caia(nval,nal),caib(nvbe,nbe),f12pq(ncore+nb,ncore+nb)
      real*8 f12ij(nbfcn,nbfcn),f12ao(nbfan,nval),f12ai(nbfan,nval)
      real*8 tscalea(nal),tscaleb(nbe),epaa(nal*(nal+1)/2)
      real*8 epbb(nbe*(nbe+1)/2),epab(nal,nbe)
      character(len=5) scftype
      character(len=8) c8
C Calculate F12 intermediates


      write(iout,*)
      write(iout,*) 'Calculating CCSD F12 intermediates...'
      write(iout,*) '-------------------------------------'
      write(iout,*) '---------- BLOCK VERSION ------------'
      write(iout,*) '-------------------------------------'
      open(scrfile1,file='F12INTE',form='unformatted')
      write(scrfile1) ecabs,emp2f12,ecoup

      open(111,file='MPVIFAA2',access='direct',
     $recl=dfnbasis*nbfcn*ifltln,status='old')
      open(110,file='MPCF12AA',access='direct',
     $recl=dfnbasis*nbfcn*ifltln,status='old')

      nnb=ncore+nb
       if(scftype.eq.'rhf  ') then
C Restricted orbitals
C Intermediate V
C F12/r12 contribution
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12/r12 contribution...'
C block
C Ezt nagyon át kell gondolni!!!!
c     if(mod(nb,nel).ne.0)then
      if(mod(nal,nel).ne.0)then
       write(iout,*)
     $"# elem. doesn't divide last dim..." 
       nbbl = mod(nal,nel)
       nbb = 1+nal/nel
      else
       nbb = nal/nel
       nbbl = nel
      endif

      write(iout,"('G size (Mb):',14x,i6)")
     $(8*nbfcn*dfnbasis*nb)/(1024*1024)
      write(iout,"('Block size (Mb):',14x,i6)")
     $(8*nbfcn*dfnbasis*(nel+1))/(1024*1024)
      write(iout,*)'----'
      write(iout,*)'last dim. (nb): ', nb
      write(iout,*)'# blocks (nbb): ', nbb
      write(iout,*)'# elem. (nel): ', nel
      write(iout,*)'# last elem. (nbbl): ', nbbl
      write(iout,*)'----'

        call f12r12cs_bl(nbfcn,nb,dfnbasis,nal,ncore,vifaa2,cf12r12aa,
     $vv,f12ij,iout,nel)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution...'
C Calculate PPL-like contribution to intermediate V, term 1
C Build (anti)symmetrized F12 integrals
        icp=dblalloc((nnb+1)*nnb/2*(nal+1)*nal/2)
        icm=dblalloc((nnb-1)*nnb/2*(nal-1)*nal/2)
        ij=icm-1
        pq=icp-1

      do ib=1,nbb
C Read vifaa and cf12aa blocks
       kstart = (ib-1)*nel + 1

       nel2 = nel
       if(ib.eq.nbb)then
        nel2 = nbbl
       endif
       do k=kstart,kstart+nel2-1
        iel = k-(ib-1)*nel
        read(111,rec=k)vifaa2(1:nbfcn,1:dfnbasis,iel)
        read(110,rec=k)cf12aa2(1:nbfcn,1:dfnbasis,iel)
       enddo
       do jj=1,nel2 ! relative j index inside a block
        j = jj + (ib - 1)*nel ! absolute j index
           do i=1,j-1
C Check whether i is in memory
         ii = i - (ib - 1)*nel
         if(i.lt.kstart)then ! i is not in memory if it is lower than the start of the block 
          ii = nel+1
          read(111,rec=i)vifaa2(1:nbfcn,1:dfnbasis,ii)
          read(110,rec=i)cf12aa2(1:nbfcn,1:dfnbasis,ii)
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            call dgemm('n','t',nnb,nnb,dfnbasis,1d0, vifaa2(1,1,ii),
     $nbfcn,cf12aa2(1,1,jj),nbfcn,0d0,f12pq,nnb)
            call dgemm('n','t',nnb,nnb,dfnbasis,1d0,cf12aa2(1,1,ii),
     $nbfcn, vifaa2(1,1,jj),nbfcn,1d0,f12pq,nnb)
            call dscal(nnb*nnb,1d0/4d0,f12pq,1)
             do q=1,nnb
               do p=1,q
                pq=pq+1
                dcore(pq)= f12pq(p,q)+f12pq(q,p)
               enddo
              dcore(pq)=dcore(pq)/2d0
             enddo
             do q=1,nnb
               do p=1,q-1
                ij=ij+1
                dcore(ij)=(f12pq(p,q)-f12pq(q,p))/2d0
               enddo
             enddo
           enddo ! i
          call dsyr2k('U','N',nnb,dfnbasis,0.5d0,vifaa2(1,1,jj),nbfcn,
     $cf12aa2(1,1,jj),nbfcn,0.d0,f12pq,nnb)
            do q=1,nnb
             do p=1,q
              pq=pq+1
              dcore(pq)=f12pq(p,q)
             enddo
             dcore(pq)=dcore(pq)/2d0
            enddo
           enddo ! j
         enddo ! ib

C Loop over blocks
        maxmem=maxcor-(imem-imem1)
        slo=1
        pq=0
         do sup=1,nb
          rsdims=(sup-slo+1)*(slo+sup)/2
          rsdima=(sup-slo+1)*(slo+sup-2)/2
           if(sup.eq.nb.or.(rsdims+sup+1)*((nnb+1)*nnb/2+(nal+1)*nal/2)+
     $(rsdima+sup)*(nnb-1)*nnb/2.gt.maxmem) then
C ki kell bővíteni az ifet hogy a vifaa az slotol supig is beleferjen
             if(verbosity.ge.3) then
              pq=pq+1
              write(c8,'(i8)') pq
              write(iout,*) 'Term 1, block '//trim(adjustl(c8))//'...'
             endif
            ivs=dblalloc((nnb+1)*nnb/2*rsdims)
            iva=dblalloc((nnb-1)*nnb/2*rsdima)
            imn=dblalloc((nal+1)*nal/2*rsdims)
C slo,sup --> nb*nb blokkositva, vv a kimenete nb*nb*elektron par kombinaciok, fent eloall itt pedig rairodik valami kontribució
            write(iout,*) 'slo,sup: ',slo,sup
c           write(iout,*) 'entering vterm1_bl ', nel,nal
c           write(iout,*) 'ncore+nb*ncore+nb: ',(ncore+nb)*(ncore+nb)
c     write(iout,*) 'first ad outside: ',loc(f12pq(1,1))
c     write(iout,*) 'last ad outside: ',loc(f12pq(ncore+nb,ncore+nb))
            call vterm1_bl(nbfcn,nb,nnb,dfnbasis,nal,vifaa2,dcore(icp),
     $dcore(icm),dcore(ivs),dcore(iva),vv,f12pq,dcore(imn),dcore(imn),
     $rsdims,rsdima,slo,sup,nel,iout)
c           write(iout,*) 'leaving vterm1_bl '
            call dbldealloc(ivs)
            slo=sup+1
             if(verbosity.ge.3) call timer
           endif
         enddo

        call dbldealloc(icp)
C Calculate PPL-like contribution to intermediate V, term 2
C Build (anti)symmetrized F12 integrals
        icp=dblalloc(nbfan*(ncore+nal)*(nal+1)*nal/2)
        icm=dblalloc(nbfan*(ncore+nal)*(nal-1)*nal/2)
        aodim=nbfan*(ncore+nal)
        ijs=icp
        ija=icm

      do ib=1,nbb
C Read vifaa and cf12aa blocks
       kstart = (ib-1)*nel + 1

       nel2 = nel
       if(ib.eq.nbb)then
        nel2 = nbbl
       endif

       do k=kstart,kstart+nel2-1
        iel = k-(ib-1)*nel
        read(111,rec=k)vifaa2(1:nbfcn,1:dfnbasis,iel)
        read(110,rec=k)cf12aa2(1:nbfcn,1:dfnbasis,iel)
       enddo
       do jj=1,nel2 ! relative j index inside a block
        j = jj + (ib - 1)*nel ! absolute j index
           do i=1,j-1
C Check whether i is in memory
         ii = i - (ib - 1)*nel
         if(i.lt.kstart)then ! i is not in memory if it is lower than the start of the block 
          ii = nel+1
          read(111,rec=i)vifaa2(1:nbfcn,1:dfnbasis,ii)
          read(110,rec=i)cf12aa2(1:nbfcn,1:dfnbasis,ii)
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.25d0,
     $ vifaa2(nnb+1,1,ii),nbfcn,cf12aa2(1,1,jj),nbfcn,0d0,dcore(ijs),
     $ nbfan)
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.25d0,
     $cf12aa2(nnb+1,1,ii),nbfcn, vifaa2(1,1,jj),nbfcn,1d0,dcore(ijs),
     $ nbfan)
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $ vifaa2(nnb+1,1,jj),nbfcn,cf12aa2(1,1,ii),nbfcn,0d0,f12ao,nbfan)
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $cf12aa2(nnb+1,1,jj),nbfcn, vifaa2(1,1,ii),nbfcn,1d0,f12ao,nbfan)
            call dcopy(aodim,dcore(ijs),1,dcore(ija),1)
            call dscal(aodim,0.5d0,dcore(ija),1)
            call daxpy(aodim,-0.125d0,f12ao,1,dcore(ija),1)
            call daxpy(aodim,0.25d0,f12ao,1,dcore(ijs),1)
            ijs=ijs+aodim
            ija=ija+aodim
           enddo
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.5d0,
     $ vifaa2(nnb+1,1,jj),nbfcn,cf12aa2(1,1,jj),nbfcn,0d0,dcore(ijs),
     $ nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.5d0,
     $cf12aa2(nnb+1,1,jj),nbfcn, vifaa2(1,1,jj),nbfcn,1d0,dcore(ijs),
     $ nbfan)
          ijs=ijs+aodim
         enddo
         enddo

C Loop over blocks
        maxmem=maxcor-(imem-imem1)
        slo=1
        pq=0

         do sup=1,nb
          rsdims=(sup-slo+1)*(slo+sup)/2
          rsdima=(sup-slo+1)*(slo+sup-2)/2
           if(sup.eq.nb.or.(rsdims+sup+1)*(aodim+(nal+1)*nal/2)+
     $(rsdima+sup)*aodim.gt.maxmem) then
             if(verbosity.ge.3) then
              pq=pq+1
              write(c8,'(i8)') pq
              write(iout,*) 'Term 2, block '//trim(adjustl(c8))//'...'
             endif
            ivs=dblalloc(aodim*rsdims)
            iva=dblalloc(aodim*rsdima)
            imn=dblalloc((nal+1)*nal/2*rsdims)
            call vterm2_bl(nbfcn,nb,nnb,dfnbasis,nal,vifaa2,dcore(icp),
     $dcore(icm),dcore(ivs),dcore(iva),vv,f12ao,dcore(imn),dcore(imn),
     $rsdims,rsdima,slo,sup,ncore,nbfan,aodim,nel)
            call dbldealloc(ivs)
            slo=sup+1
             if(verbosity.ge.3) call timer
           endif
         enddo

        call dbldealloc(icp)
C Intermediate Vpi
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,k,tmp)

         do p=1,nb
           do i=1,nal
            tmp=vv(p,i,i*(i-1)/2+i)
             do k=1,i-1
              tmp=tmp+2d0*vv(k,p,i*(i-1)/2+k)-vv(p,k,i*(i-1)/2+k)
             enddo
             do k=i+1,nal
              tmp=tmp+2d0*vv(p,k,k*(k-1)/2+i)-vv(k,p,k*(k-1)/2+i)
             enddo
            vpia(p,i)=tmp
           enddo
         enddo

C$OMP END PARALLEL DO
C Intermediates Uai and Cai
         if(verbosity.ge.3) write(iout,*) 'Intermediate U...'
        iva=dblalloc(nbfan*nal)
        uaia(:,:)=vpia(nal+1:nal+nval,:)
        caia=0d0

      do ib=1,nbb
C Read vifaa and cf12aa blocks
       kstart = (ib-1)*nel + 1

       nel2 = nel
       if(ib.eq.nbb)then
        nel2 = nbbl
       endif

       do k=kstart,kstart+nel2-1
        iel = k-(ib-1)*nel
        read(111,rec=k)vifaa2(1:nbfcn,1:dfnbasis,iel)
        read(110,rec=k)cf12aa2(1:nbfcn,1:dfnbasis,iel)
       enddo
       do jj=1,nel2 ! relative j index inside a block
        j = jj + (ib - 1)*nel ! absolute j index
           do i=1,j
C Check whether i is in memory
         ii = i - (ib - 1)*nel
         if(i.lt.kstart)then ! i is not in memory if it is lower than the start of the block 
          ii = nel+1
          read(111,rec=i)vifaa2(1:nbfcn,1:dfnbasis,ii)
          read(110,rec=i)cf12aa2(1:nbfcn,1:dfnbasis,ii)
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $vifaa2(ncore+nb+1,1,ii),nbfcn,cf12aa2(ncore+nal+1,1,jj),nbfcn,0d0,
     $f12ai,nbfan)
            call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa2(ncore+nb+1,1,ii),nbfcn,vifaa2(ncore+nal+1,1,jj),nbfcn,1d0,
     $f12ai,nbfan)
            call dgemm('n','t',nbfan,nval,dfnbasis,-1d0,
     $vifaa2(ncore+nb+1,1,jj),nbfcn,cf12aa2(ncore+nal+1,1,ii),nbfcn,0d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nval,dfnbasis,-1d0,
     $cf12aa2(ncore+nb+1,1,jj),nbfcn,vifaa2(ncore+nal+1,1,ii),nbfcn,1d0,
     $f12ao,nbfan)
            call daxpy(nbfan*nval,5d0,f12ai,1,f12ao,1)
            call dgemv('t',nbfan,nval,1d0/8d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+i),1,1d0,caia(1,j),1)
            call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $vifaa2(ncore+nb+1,1,ii),nbfcn,vifaa2(ncore+1,1,jj),nbfcn,0d0,
     $dcore(iva),nbfan)
            call dgemm('t','n',nval,nal,nbfan,-1d0/4d0,f12ao,nbfan,
     $dcore(iva),nbfan,1d0,uaia,nval)
             if(i.ne.j) then
              call daxpy(nbfan*nval,-5d0,f12ai,1,f12ao,1)
              call dscal(nbfan*nval,-5d0,f12ao,1)
              call daxpy(nbfan*nval,-1d0,f12ai,1,f12ao,1)
              call dgemv('t',nbfan,nval,1d0/8d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+j),1,1d0,caia(1,i),1)
              call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $vifaa2(ncore+nb+1,1,jj),nbfcn,vifaa2(ncore+1,1,ii),nbfcn,0d0,
     $dcore(iva),nbfan)
              call dgemm('t','n',nval,nal,nbfan,-1d0/4d0,f12ao,nbfan,
     $dcore(iva),nbfan,1d0,uaia,nval)
             endif
           enddo
           enddo
         enddo

        call daxpy(nval*nal,1d0,uaia,1,caia,1)
        call dbldealloc(iva)
         if(verbosity.ge.3) call timer
C Intermediate Cabij
         if(verbosity.ge.3) write(iout,*) 'Intermediate C...'
        icabij=dblalloc(nval*nval*(nal+1)*nal/2)
        ivs=dblalloc(nbfan*nval*nal)
        iva=dblalloc(nbfan*nval*nal)
        imn=dblalloc(nbfan*nval*nal)

      freememspace = maxcor-imem+imem1
      if (nbfan*nval*nal.gt.freememspace) then
        call cabijcalc_bl(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,
     $cf12aa2,cf12r12aa,dcore(ivs),dcore(iva),dcore(imn),f12ao,focka,
     $dcore(icabij),vv,vifaa2,nel,0,dcore(imn),nbb,nbbl)
C                                  ---^ this is dummy
      else
        write(iout,*) 
     $'Mem. suffices, G^{P}_{ia}G^{P}_{jb} will be done with one I/O...'
        imn2=dblalloc(nbfan*nval*nal)
        call cabijcalc_bl(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,
     $cf12aa2,cf12r12aa,dcore(ivs),dcore(iva),dcore(imn),f12ao,focka,
     $dcore(icabij),vv,vifaa2,nel,1,dcore(imn2),nbb,nbbl)
      endif


        call dbldealloc(ivs)
         if(verbosity.ge.3) call timer
C Save intermediates
         if(verbosity.ge.3) write(iout,*) 'Saving intermediates...'
        call saveimcs(nb,nal,nval,scrfile1,dcore(icabij),vv,uaia,vpia,
     $caia,tscalea,tfact,epaa)
        call dbldealloc(icabij)
       else
C Unrestricted orbitals
C Intermediate V
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12/r12 contribution...'
        call f12r12os(nbfcn,nb,dfnbasis,nal,nbe,ncore,vifaa,vifbb,vifab,
     $vifba,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,vaa,vbb,vab,f12ij)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution, alpha-alpha block...'
        call vaacalc(nb,nnb,ncore,nal,nbfcn,nbfan,dfnbasis,dcore,
     $imem,imem1,maxcor,vifaa,cf12aa,f12pq,vaa)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution, beta-beta block...'
        call vaacalc(nb,nnb,ncore,nbe,nbfcn,nbfan,dfnbasis,dcore,
     $imem,imem1,maxcor,vifbb,cf12bb,f12pq,vbb)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution, alpha-beta block...'
        call vabcalc(nb,nnb,ncore,nal,nbe,nbfcn,nbfan,dfnbasis,
     $dcore,imem,imem1,maxcor,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,
     $cf12ab,cf12ba,vab)
         if(verbosity.ge.3) call timer
C Intermediate Vpi
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,k,tmp)
         do p=1,nb
           do i=1,nal
            tmp=0.d0
             do k=1,min(i-1,p-1)
              tmp=tmp+vaa((p-1)*(p-2)/2+k,(i-1)*(i-2)/2+k)
             enddo
             do k=i+1,min(nal,p-1)
              tmp=tmp-vaa((p-1)*(p-2)/2+k,(k-1)*(k-2)/2+i)
             enddo
             do k=max(1,p+1),i-1
              tmp=tmp-vaa((k-1)*(k-2)/2+p,(i-1)*(i-2)/2+k)
             enddo
             do k=max(i+1,p+1),nal
              tmp=tmp+vaa((k-1)*(k-2)/2+p,(k-1)*(k-2)/2+i)
             enddo
             do k=1,nbe
              tmp=tmp+vab(p,k,i,k)
             enddo
            vpia(p,i)=tmp
           enddo
         enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,k,tmp)
         do p=1,nb
           do i=1,nbe
            tmp=0.d0
             do k=1,min(i-1,p-1)
              tmp=tmp+vbb((p-1)*(p-2)/2+k,(i-1)*(i-2)/2+k)
             enddo
             do k=i+1,min(nbe,p-1)
              tmp=tmp-vbb((p-1)*(p-2)/2+k,(k-1)*(k-2)/2+i)
             enddo
             do k=max(1,p+1),i-1
              tmp=tmp-vbb((k-1)*(k-2)/2+p,(i-1)*(i-2)/2+k)
             enddo
             do k=max(i+1,p+1),nbe
              tmp=tmp+vbb((k-1)*(k-2)/2+p,(k-1)*(k-2)/2+i)
             enddo
             do k=1,nal
              tmp=tmp+vab(k,p,k,i)
             enddo
            vpib(p,i)=tmp
           enddo
         enddo
C$OMP END PARALLEL DO
C Intermediates Uai and Cai
         if(verbosity.ge.3) write(iout,*) 'Intermediate U...'
        call uaiaacalc(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,vifaa,
     $cf12aa,focka,f12ao,f12ai,vpia,uaia,caia)
        call uaiaacalc(nb,ncore,nbe,nvbe,nbfcn,nbfan,dfnbasis,vifbb,
     $cf12bb,fockb,f12ao,f12ai,vpib,uaib,caib)
C Alpha-beta
         do j=1,nbe
           do i=1,nal
            call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nal+1,1,i),nbfcn,cf12bb(ncore+nb+1,1,j),nbfcn,0d0,
     $f12ao,nval)
            call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nal+1,1,i),nbfcn, vifbb(ncore+nb+1,1,j),nbfcn,1d0,
     $f12ao,nval)
            call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0, 
     $ vifba(ncore+nal+1,1,j),nbfcn,cf12ab(ncore+nb+1,1,i),nbfcn,1d0,
     $f12ao,nval)
            call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nal+1,1,j),nbfcn, vifab(ncore+nb+1,1,i),nbfcn,1d0,
     $f12ao,nval)
            call dgemv('n',nval,nbfan,1d0,f12ao,nval,
     $fockb(ncore+nb+1,ncore+j),1,1d0,caia(1,i),1)
            call dgemm('n','t',nal,nbfan,dfnbasis,1d0,
     $vifaa(ncore+1,1,i),nbfcn,vifbb(ncore+nb+1,1,j),nbfcn,0d0,f12ai,
     $nal)
            call dgemm('n','t',nval,nal,nbfan,-2d0,f12ao,nval,f12ai,
     $nal,1d0,uaia,nval)
           enddo
         enddo
C
         do j=1,nbe
           do i=1,nal
            call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nb+1,1,i),nbfcn,cf12bb(ncore+nbe+1,1,j),nbfcn,0d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nb+1,1,i),nbfcn, vifbb(ncore+nbe+1,1,j),nbfcn,1d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $ vifba(ncore+nb+1,1,j),nbfcn,cf12ab(ncore+nbe+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nb+1,1,j),nbfcn, vifab(ncore+nbe+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
            call dgemv('t',nbfan,nvbe,1d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+i),1,1d0,caib(1,j),1)
            call dgemm('n','t',nbfan,nbe,dfnbasis,1d0,
     $vifaa(ncore+nb+1,1,i),nbfcn,vifbb(ncore+1,1,j),nbfcn,0d0,f12ai,
     $nbfan)
            call dgemm('t','n',nvbe,nbe,nbfan,-2d0,f12ao,nbfan,f12ai,
     $nbfan,1d0,uaib,nvbe)
           enddo
         enddo
        call daxpy(nval*nal,1d0,uaia,1,caia,1)
        call daxpy(nvbe*nbe,1d0,uaib,1,caib,1)
         if(verbosity.ge.3) call timer
C Intermediate Cabij
         if(verbosity.ge.3) write(iout,*) 'Intermediate C, alpha...'
        icabijaa=dblalloc((nval-1)*nval/2*(nal-1)*nal/2)
        icabijbb=dblalloc((nvbe-1)*nvbe/2*(nbe-1)*nbe/2)
        icabijab=dblalloc(nval*nvbe*nal*nbe)
        call cabijal(nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis,
     $vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,cf12r12aa,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,cf12r12ba,cf12r12ba,f12ao,
     $focka,fockb,dcore(icabijaa),dcore(icabijbb),dcore(icabijab),f12ij,
     $f12ij,vaa,vbb,vab)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*) 'Intermediate C, beta...'
        call cabijbe(nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis,
     $vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,cf12r12aa,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,cf12r12ba,cf12r12ba,f12ao,
     $focka,fockb,dcore(icabijaa),dcore(icabijbb),dcore(icabijab),f12ij,
     $f12ij)
         if(verbosity.ge.3) call timer
C Save intermediates
         if(verbosity.ge.3) write(iout,*) 'Saving intermediates...'
        call saveimos(nb,nal,nbe,nval,nvbe,scrfile1,dcore(icabijaa),
     $dcore(icabijbb),dcore(icabijab),vaa,vbb,vab,uaia,uaib,vpia,vpib,
     $caia,caib,tscalea,tscaleb,tfact,epaa,epbb,epab)
        call dbldealloc(icabijaa)
       endif
      close(scrfile1)
      call timer
C
C Close block files
      close(111)
      close(110)
      return
      end
C
************************************************************************
      subroutine ccimed(nb,dfnbasis,nal,nbe,nval,nvbe,ncore,focka,fockb,
     $nbfcn,nbfan,emp2f12,ecabs,scrfile1,scftype,vv,vaa,vbb,vab,vpia,
     $vpib,caia,caib,uaia,uaib,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,
     $cf12ab,cf12ba,iout,dcore,imem,imem1,maxcor,f12ij,f12pq,f12ao,
     $f12ai,verbosity,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,tscalea,
     $tscaleb,tfact,epaa,epbb,epab,ecoup)
************************************************************************
* Calculate intermediates for CCSD calculation
************************************************************************
      implicit none
      integer nb,i,j,k,l,p,q,r,s,a,b,nal,nbe,ncore,o,rs,nbfcn,iout,nnb
      integer nval,nvbe,c,kl,ij,jk,ab,dfnbasis,scrfile1,pq,imem,dblalloc
      integer icp,icm,ivs,iva,imn,rsdims,slo,sup,imem1,maxcor,maxmem
      integer verbosity,nbfan,rsdima,aodim,rsdim,ijs,ija
      integer icabij,icabijaa,icabijbb,icabijab
      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn),emp2f12,ecabs,tmp
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal),tfact
      real*8 vifbb(nbfcn,dfnbasis,nb),cf12bb(nbfcn,dfnbasis,nbe),ecoup
      real*8 vifab(nbfcn,dfnbasis,nb),cf12ab(nbfcn,dfnbasis,nal)
      real*8 vifba(nbfcn,dfnbasis,nb),cf12ba(nbfcn,dfnbasis,nbe)
      real*8 cf12r12aa(nb,dfnbasis,nal),cf12r12ab(nb,dfnbasis,nal)
      real*8 cf12r12bb(nb,dfnbasis,nbe),cf12r12ba(nb,dfnbasis,nbe)
      real*8 vaa((nb-1)*nb/2,(nal-1)*nal/2),vab(nb,nb,nal,nbe),dcore(*)
      real*8 vbb((nb-1)*nb/2,(nbe-1)*nbe/2),vv(nb,nb,(nal+1)*nal/2)
      real*8 vpia(nb,nal),vpib(nb,nbe),uaia(nval,nal),uaib(nvbe,nbe)
      real*8 caia(nval,nal),caib(nvbe,nbe),f12pq(ncore+nb,ncore+nb)
      real*8 f12ij(nbfcn,nbfcn),f12ao(nbfan,nval),f12ai(nbfan,nval)
      real*8 tscalea(nal),tscaleb(nbe),epaa(nal*(nal+1)/2)
      real*8 epbb(nbe*(nbe+1)/2),epab(nal,nbe)
      character(len=5) scftype
      character(len=8) c8
C Calculate F12 intermediates
      write(iout,*)
      write(iout,*) 'Calculating CCSD F12 intermediates...'
      open(scrfile1,file='F12INTE',form='unformatted')
      write(scrfile1) ecabs,emp2f12,ecoup
      nnb=ncore+nb
       if(scftype.eq.'rhf  ') then
C Restricted orbitals
C Intermediate V
C F12/r12 contribution
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12/r12 contribution...'
        call f12r12cs(nbfcn,nb,dfnbasis,nal,ncore,vifaa,cf12r12aa,vv,
     $f12ij)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution...'
C Calculate PPL-like contribution to intermediate V, term 1
C Build (anti)symmetrized F12 integrals
        icp=dblalloc((nnb+1)*nnb/2*(nal+1)*nal/2)
        icm=dblalloc((nnb-1)*nnb/2*(nal-1)*nal/2)
        ij=icm-1
        pq=icp-1
C ezt itt ugyan úgy lehet blokkosítani mint eddig
         do j=1,nal
           do i=1,j-1
            call dgemm('n','t',nnb,nnb,dfnbasis,1d0, vifaa(1,1,i),
     $nbfcn,cf12aa(1,1,j),nbfcn,0d0,f12pq,nnb)
            call dgemm('n','t',nnb,nnb,dfnbasis,1d0,cf12aa(1,1,i),
     $nbfcn, vifaa(1,1,j),nbfcn,1d0,f12pq,nnb)
            call dscal(nnb*nnb,1d0/4d0,f12pq,1)
             do q=1,nnb
               do p=1,q
                pq=pq+1
                dcore(pq)= f12pq(p,q)+f12pq(q,p)
               enddo
              dcore(pq)=dcore(pq)/2d0
             enddo
             do q=1,nnb
               do p=1,q-1
                ij=ij+1
                dcore(ij)=(f12pq(p,q)-f12pq(q,p))/2d0
               enddo
             enddo
           enddo
          call dsyr2k('U','N',nnb,dfnbasis,0.5d0,vifaa(1,1,j),nbfcn,
     $cf12aa(1,1,j),nbfcn,0.d0,f12pq,nnb)
           do q=1,nnb
             do p=1,q
              pq=pq+1
              dcore(pq)=f12pq(p,q)
             enddo
            dcore(pq)=dcore(pq)/2d0
           enddo
         enddo
C Loop over blocks
        maxmem=maxcor-(imem-imem1)
        slo=1
        pq=0
         do sup=1,nb
          rsdims=(sup-slo+1)*(slo+sup)/2
          rsdima=(sup-slo+1)*(slo+sup-2)/2
           if(sup.eq.nb.or.(rsdims+sup+1)*((nnb+1)*nnb/2+(nal+1)*nal/2)+
     $(rsdima+sup)*(nnb-1)*nnb/2.gt.maxmem) then
             if(verbosity.ge.3) then
              pq=pq+1
              write(c8,'(i8)') pq
              write(iout,*) 'Term 1, block '//trim(adjustl(c8))//'...'
             endif
            ivs=dblalloc((nnb+1)*nnb/2*rsdims)
            iva=dblalloc((nnb-1)*nnb/2*rsdima)
            imn=dblalloc((nal+1)*nal/2*rsdims)
            call vterm1(nbfcn,nb,nnb,dfnbasis,nal,vifaa,dcore(icp),
     $dcore(icm),dcore(ivs),dcore(iva),vv,f12pq,dcore(imn),dcore(imn),
     $rsdims,rsdima,slo,sup)
            call dbldealloc(ivs)
            slo=sup+1
             if(verbosity.ge.3) call timer
           endif
         enddo
        call dbldealloc(icp)
C Calculate PPL-like contribution to intermediate V, term 2
C Build (anti)symmetrized F12 integrals
        icp=dblalloc(nbfan*(ncore+nal)*(nal+1)*nal/2)
        icm=dblalloc(nbfan*(ncore+nal)*(nal-1)*nal/2)
        aodim=nbfan*(ncore+nal)
        ijs=icp
        ija=icm
         do j=1,nal
           do i=1,j-1
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.25d0,
     $ vifaa(nnb+1,1,i),nbfcn,cf12aa(1,1,j),nbfcn,0d0,dcore(ijs),nbfan)
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.25d0,
     $cf12aa(nnb+1,1,i),nbfcn, vifaa(1,1,j),nbfcn,1d0,dcore(ijs),nbfan)
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $ vifaa(nnb+1,1,j),nbfcn,cf12aa(1,1,i),nbfcn,0d0,f12ao,nbfan)
            call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $cf12aa(nnb+1,1,j),nbfcn, vifaa(1,1,i),nbfcn,1d0,f12ao,nbfan)
            call dcopy(aodim,dcore(ijs),1,dcore(ija),1)
            call dscal(aodim,0.5d0,dcore(ija),1)
            call daxpy(aodim,-0.125d0,f12ao,1,dcore(ija),1)
            call daxpy(aodim,0.25d0,f12ao,1,dcore(ijs),1)
            ijs=ijs+aodim
            ija=ija+aodim
           enddo
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.5d0,
     $ vifaa(nnb+1,1,j),nbfcn,cf12aa(1,1,j),nbfcn,0d0,dcore(ijs),nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,0.5d0,
     $cf12aa(nnb+1,1,j),nbfcn, vifaa(1,1,j),nbfcn,1d0,dcore(ijs),nbfan)
          ijs=ijs+aodim
         enddo
C Loop over blocks
        maxmem=maxcor-(imem-imem1)
        slo=1
        pq=0
         do sup=1,nb
          rsdims=(sup-slo+1)*(slo+sup)/2
          rsdima=(sup-slo+1)*(slo+sup-2)/2
           if(sup.eq.nb.or.(rsdims+sup+1)*(aodim+(nal+1)*nal/2)+
     $(rsdima+sup)*aodim.gt.maxmem) then
             if(verbosity.ge.3) then
              pq=pq+1
              write(c8,'(i8)') pq
              write(iout,*) 'Term 2, block '//trim(adjustl(c8))//'...'
             endif
            ivs=dblalloc(aodim*rsdims)
            iva=dblalloc(aodim*rsdima)
            imn=dblalloc((nal+1)*nal/2*rsdims)
            call vterm2(nbfcn,nb,nnb,dfnbasis,nal,vifaa,dcore(icp),
     $dcore(icm),dcore(ivs),dcore(iva),vv,f12ao,dcore(imn),dcore(imn),
     $rsdims,rsdima,slo,sup,ncore,nbfan,aodim)
            call dbldealloc(ivs)
            slo=sup+1
             if(verbosity.ge.3) call timer
           endif
         enddo
        call dbldealloc(icp)
C Intermediate Vpi
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,k,tmp)
         do p=1,nb
           do i=1,nal
            tmp=vv(p,i,i*(i-1)/2+i)
             do k=1,i-1
              tmp=tmp+2d0*vv(k,p,i*(i-1)/2+k)-vv(p,k,i*(i-1)/2+k)
             enddo
             do k=i+1,nal
              tmp=tmp+2d0*vv(p,k,k*(k-1)/2+i)-vv(k,p,k*(k-1)/2+i)
             enddo
            vpia(p,i)=tmp
           enddo
         enddo
C$OMP END PARALLEL DO
C Intermediates Uai and Cai
         if(verbosity.ge.3) write(iout,*) 'Intermediate U...'
        iva=dblalloc(nbfan*nal)
        uaia(:,:)=vpia(nal+1:nal+nval,:)
        caia=0d0
         do j=1,nal
           do i=1,j
            call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $ vifaa(ncore+nb+1,1,i),nbfcn,cf12aa(ncore+nal+1,1,j),nbfcn,0d0,
     $f12ai,nbfan)
            call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa(ncore+nb+1,1,i),nbfcn, vifaa(ncore+nal+1,1,j),nbfcn,1d0,
     $f12ai,nbfan)
            call dgemm('n','t',nbfan,nval,dfnbasis,-1d0,
     $ vifaa(ncore+nb+1,1,j),nbfcn,cf12aa(ncore+nal+1,1,i),nbfcn,0d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nval,dfnbasis,-1d0,
     $cf12aa(ncore+nb+1,1,j),nbfcn, vifaa(ncore+nal+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
            call daxpy(nbfan*nval,5d0,f12ai,1,f12ao,1)
            call dgemv('t',nbfan,nval,1d0/8d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+i),1,1d0,caia(1,j),1)
            call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $vifaa(ncore+nb+1,1,i),nbfcn,vifaa(ncore+1,1,j),nbfcn,0d0,
     $dcore(iva),nbfan)
            call dgemm('t','n',nval,nal,nbfan,-1d0/4d0,f12ao,nbfan,
     $dcore(iva),nbfan,1d0,uaia,nval)
             if(i.ne.j) then
              call daxpy(nbfan*nval,-5d0,f12ai,1,f12ao,1)
              call dscal(nbfan*nval,-5d0,f12ao,1)
              call daxpy(nbfan*nval,-1d0,f12ai,1,f12ao,1)
              call dgemv('t',nbfan,nval,1d0/8d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+j),1,1d0,caia(1,i),1)
              call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $vifaa(ncore+nb+1,1,j),nbfcn,vifaa(ncore+1,1,i),nbfcn,0d0,
     $dcore(iva),nbfan)
              call dgemm('t','n',nval,nal,nbfan,-1d0/4d0,f12ao,nbfan,
     $dcore(iva),nbfan,1d0,uaia,nval)
             endif
           enddo
         enddo
        call daxpy(nval*nal,1d0,uaia,1,caia,1)
        call dbldealloc(iva)
         if(verbosity.ge.3) call timer
C Intermediate Cabij
         if(verbosity.ge.3) write(iout,*) 'Intermediate C...'
        icabij=dblalloc(nval*nval*(nal+1)*nal/2)
        ivs=dblalloc(nbfan*nval*nal)
        iva=dblalloc(nbfan*nval*nal)
        imn=dblalloc(nbfan*nval*nal)
        call cabijcalc(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,vifaa,
     $cf12aa,cf12r12aa,dcore(ivs),dcore(iva),dcore(imn),f12ao,focka,
     $dcore(icabij),vv)
        call dbldealloc(ivs)
         if(verbosity.ge.3) call timer
C Save intermediates
         if(verbosity.ge.3) write(iout,*) 'Saving intermediates...'
        call saveimcs(nb,nal,nval,scrfile1,dcore(icabij),vv,uaia,vpia,
     $caia,tscalea,tfact,epaa)
        call dbldealloc(icabij)
       else
C Unrestricted orbitals
C Intermediate V
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12/r12 contribution...'
        call f12r12os(nbfcn,nb,dfnbasis,nal,nbe,ncore,vifaa,vifbb,vifab,
     $vifba,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,vaa,vbb,vab,f12ij)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution, alpha-alpha block...'
        call vaacalc(nb,nnb,ncore,nal,nbfcn,nbfan,dfnbasis,dcore,
     $imem,imem1,maxcor,vifaa,cf12aa,f12pq,vaa)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution, beta-beta block...'
        call vaacalc(nb,nnb,ncore,nbe,nbfcn,nbfan,dfnbasis,dcore,
     $imem,imem1,maxcor,vifbb,cf12bb,f12pq,vbb)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*)
     $'Intermediate V, F12 contribution, alpha-beta block...'
        call vabcalc(nb,nnb,ncore,nal,nbe,nbfcn,nbfan,dfnbasis,
     $dcore,imem,imem1,maxcor,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,
     $cf12ab,cf12ba,vab)
         if(verbosity.ge.3) call timer
C Intermediate Vpi
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,k,tmp)
         do p=1,nb
           do i=1,nal
            tmp=0.d0
             do k=1,min(i-1,p-1)
              tmp=tmp+vaa((p-1)*(p-2)/2+k,(i-1)*(i-2)/2+k)
             enddo
             do k=i+1,min(nal,p-1)
              tmp=tmp-vaa((p-1)*(p-2)/2+k,(k-1)*(k-2)/2+i)
             enddo
             do k=max(1,p+1),i-1
              tmp=tmp-vaa((k-1)*(k-2)/2+p,(i-1)*(i-2)/2+k)
             enddo
             do k=max(i+1,p+1),nal
              tmp=tmp+vaa((k-1)*(k-2)/2+p,(k-1)*(k-2)/2+i)
             enddo
             do k=1,nbe
              tmp=tmp+vab(p,k,i,k)
             enddo
            vpia(p,i)=tmp
           enddo
         enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,i,k,tmp)
         do p=1,nb
           do i=1,nbe
            tmp=0.d0
             do k=1,min(i-1,p-1)
              tmp=tmp+vbb((p-1)*(p-2)/2+k,(i-1)*(i-2)/2+k)
             enddo
             do k=i+1,min(nbe,p-1)
              tmp=tmp-vbb((p-1)*(p-2)/2+k,(k-1)*(k-2)/2+i)
             enddo
             do k=max(1,p+1),i-1
              tmp=tmp-vbb((k-1)*(k-2)/2+p,(i-1)*(i-2)/2+k)
             enddo
             do k=max(i+1,p+1),nbe
              tmp=tmp+vbb((k-1)*(k-2)/2+p,(k-1)*(k-2)/2+i)
             enddo
             do k=1,nal
              tmp=tmp+vab(k,p,k,i)
             enddo
            vpib(p,i)=tmp
           enddo
         enddo
C$OMP END PARALLEL DO
C Intermediates Uai and Cai
         if(verbosity.ge.3) write(iout,*) 'Intermediate U...'
        call uaiaacalc(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,vifaa,
     $cf12aa,focka,f12ao,f12ai,vpia,uaia,caia)
        call uaiaacalc(nb,ncore,nbe,nvbe,nbfcn,nbfan,dfnbasis,vifbb,
     $cf12bb,fockb,f12ao,f12ai,vpib,uaib,caib)
C Alpha-beta
         do j=1,nbe
           do i=1,nal
            call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nal+1,1,i),nbfcn,cf12bb(ncore+nb+1,1,j),nbfcn,0d0,
     $f12ao,nval)
            call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nal+1,1,i),nbfcn, vifbb(ncore+nb+1,1,j),nbfcn,1d0,
     $f12ao,nval)
            call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0, 
     $ vifba(ncore+nal+1,1,j),nbfcn,cf12ab(ncore+nb+1,1,i),nbfcn,1d0,
     $f12ao,nval)
            call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nal+1,1,j),nbfcn, vifab(ncore+nb+1,1,i),nbfcn,1d0,
     $f12ao,nval)
            call dgemv('n',nval,nbfan,1d0,f12ao,nval,
     $fockb(ncore+nb+1,ncore+j),1,1d0,caia(1,i),1)
            call dgemm('n','t',nal,nbfan,dfnbasis,1d0,
     $vifaa(ncore+1,1,i),nbfcn,vifbb(ncore+nb+1,1,j),nbfcn,0d0,f12ai,
     $nal)
            call dgemm('n','t',nval,nal,nbfan,-2d0,f12ao,nval,f12ai,
     $nal,1d0,uaia,nval)
           enddo
         enddo
C
         do j=1,nbe
           do i=1,nal
            call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nb+1,1,i),nbfcn,cf12bb(ncore+nbe+1,1,j),nbfcn,0d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nb+1,1,i),nbfcn, vifbb(ncore+nbe+1,1,j),nbfcn,1d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $ vifba(ncore+nb+1,1,j),nbfcn,cf12ab(ncore+nbe+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
            call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nb+1,1,j),nbfcn, vifab(ncore+nbe+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
            call dgemv('t',nbfan,nvbe,1d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+i),1,1d0,caib(1,j),1)
            call dgemm('n','t',nbfan,nbe,dfnbasis,1d0,
     $vifaa(ncore+nb+1,1,i),nbfcn,vifbb(ncore+1,1,j),nbfcn,0d0,f12ai,
     $nbfan)
            call dgemm('t','n',nvbe,nbe,nbfan,-2d0,f12ao,nbfan,f12ai,
     $nbfan,1d0,uaib,nvbe)
           enddo
         enddo
        call daxpy(nval*nal,1d0,uaia,1,caia,1)
        call daxpy(nvbe*nbe,1d0,uaib,1,caib,1)
         if(verbosity.ge.3) call timer
C Intermediate Cabij
         if(verbosity.ge.3) write(iout,*) 'Intermediate C, alpha...'
        icabijaa=dblalloc((nval-1)*nval/2*(nal-1)*nal/2)
        icabijbb=dblalloc((nvbe-1)*nvbe/2*(nbe-1)*nbe/2)
        icabijab=dblalloc(nval*nvbe*nal*nbe)
        call cabijal(nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis,
     $vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,cf12r12aa,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,cf12r12ba,cf12r12ba,f12ao,
     $focka,fockb,dcore(icabijaa),dcore(icabijbb),dcore(icabijab),f12ij,
     $f12ij,vaa,vbb,vab)
         if(verbosity.ge.3) call timer
         if(verbosity.ge.3) write(iout,*) 'Intermediate C, beta...'
        call cabijbe(nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis,
     $vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,cf12r12aa,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,cf12r12ba,cf12r12ba,f12ao,
     $focka,fockb,dcore(icabijaa),dcore(icabijbb),dcore(icabijab),f12ij,
     $f12ij)
         if(verbosity.ge.3) call timer
C Save intermediates
         if(verbosity.ge.3) write(iout,*) 'Saving intermediates...'
        call saveimos(nb,nal,nbe,nval,nvbe,scrfile1,dcore(icabijaa),
     $dcore(icabijbb),dcore(icabijab),vaa,vbb,vab,uaia,uaib,vpia,vpib,
     $caia,caib,tscalea,tscaleb,tfact,epaa,epbb,epab)
        call dbldealloc(icabijaa)
       endif
      close(scrfile1)
      call timer
C
      return
      end
C
************************************************************************
      subroutine mp2int_bl(nbf,nalc,nbec,ncore,nb,dfnbasis,nbfc,nbfcn,
     $cmoa,
     $cmob,epaa,epbb,epab,focka,fockb,hpja,hpjb,dcore,imem,natoms,iout,
     $imem1,maxcor,nangmax,ncontrmax,nprimmax,ncartmax,cartg,nsphermax,
     $nmboys,itol,nbfshmax,nang,ncontr,nprim,gexp,gcoef,coord,ctostr,cf,
     $boysval,indarr,gcn,pre,nshrange,thad,thcf2,scoord,rqqij,rqqkl,
     $hrec,nbset,spctostr,dfipre,icore,tcdfint,dfrqq,spre,natrange,
     $dfipra,logkc,kp,gck,cprea,cpreb,tcf12int,vifaa,vifbb,vifab,cf12aa,
     $cf12bb,cf12ab,cf122aa,cf122bb,cf122ab,vifba,cf12ba,cf122ba,
     $scrfile2,scftype,lcc,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,
     $verbosity,ifltln)
************************************************************************
c Calculate and transform two- and three-center integrals
c BLOCK version
************************************************************************
      implicit none
      integer nbf,nalc,nbec,ncore,dfnbasis,imem,natoms,iout,imem1,i,nb
      integer dblalloc,intalloc,icore(*),nbfc,nbfcn,scrfile2,nn,p,q,r,s
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp,maxcor
      integer nshrange(2,0:nangmax,natoms,nbset),ivhaia,ivhaib,ihai
      integer natrange(2,natoms,nbset),if12ij,iscr,nalx,nbex,verbosity
      integer ivifaa,j,k,ivifaa2,ifltln,freememspace,nlm,nlen,nbbl,nbb
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2),dfipra(natoms)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cprea(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),dfrqq,spre,dfipre(natoms,0:nangmax)
      real*8 vifaa(dfnbasis,nbfcn,nbfcn),vifbb(dfnbasis,nbfcn,nbfcn)
      real*8 vifab(dfnbasis,nbfcn,nbfcn),vifba(dfnbasis,nbfcn,nbfcn)
      real*8 cpreb(natoms,0:nangmax),tcf12int(dfnbasis,dfnbasis)
      real*8 cf12aa(nbfcn,dfnbasis,*),cf12bb(nbfcn,dfnbasis,*)
      real*8 cf12ab(nbfcn,dfnbasis,*),cf12ba(nbfcn,dfnbasis,*)
      real*8 cf122aa(nbfcn,dfnbasis,*),cf122bb(nbfcn,dfnbasis,*)
      real*8 cf122ab(nbfcn,dfnbasis,*),cf122ba(nbfcn,dfnbasis,*)
      real*8 epaa(*),epbb(*),epab(*)
      real*8 hpja(nbfcn,nbfcn),hpjb(nbfcn,nbfcn)
      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn)
      real*8 cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba
      character(len=5) scftype
      logical cartg,lcc,bver
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
      interface
      subroutine boys
      end
      subroutine boysf12
      end
      subroutine boysf122
      end
      subroutine boysdelf122
      end
      subroutine boysf12r12
      end
      end interface
C MO coefficient prescreening
       if(lcc) then
        nalx=nb
        nbex=nb
       else
        nalx=nalc
        nbex=nbec
       endif
      cprea=0.d0
      call dfcpre(cmoa(1,ncore+1),cprea,natoms,nangmax,ncontrmax,
     $nsphermax,cartg,indarr,nang,ncontr,nbfc,nalx,.true.)
       if(scftype.ne.'rhf  ') then
        cpreb=0.d0
        call dfcpre(cmob(1,ncore+1),cpreb,natoms,nangmax,ncontrmax,
     $nsphermax,cartg,indarr,nang,ncontr,nbfc,nbex,.true.)
       endif
C Two-center Coulomb integrals
      call df2int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,nmboys,
     $dcore,imem,0,0.00001d0*itol,nshrange(1,0,1,3),iout,imem1,maxcor,
     $thad,thcf2,scoord,rqqij,rqqkl,0,1,nang(1,3),ncontr(0,1,3),
     $nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfrqq,dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),tcdfint,dfipre,i,i,.true.,
     $scrfile2,tcdfint,spctostr,1,dfipra,tcmax,.false.,tcdfint,
     $0.d0,boys)
      iisyev=0
      call dpptrf('L',dfnbasis,tcdfint,isyev)
       if(isyev.ne.0) then
        write(iout,*) 'Fatal error at the Cholesky decomposition!'
        call mrccend(1)
       endif
C Three-center Coulomb integrals + first half transformation
      ihai=dblalloc(dfnbasis*max(4,nalc,nbec)*nbfc)

      write(iout,"('Bottleneck Haia size (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*nalx)/(1024*1024)

      freememspace = maxcor-imem+imem1
      bver=.false.
c     if(freememspace.lt.nbfc*dfnbasis*nalx*2)then
      if(.false.)then
       bver=.true.
       write(iout,*) 'Haia too big for RAM...'
C TODO calculate this
       nlm = 50
      write(iout,"('Memalloc attempt for Haia block (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*nlm)/(1024*1024)
       ivhaia=dblalloc(nbfc*dfnbasis*nlm)
      else
C ezt még számolni kell...
       nlen = nalc
       ivhaia=dblalloc(nbfc*dfnbasis*nlen)
      endif

       if(scftype.ne.'rhf  ') then
        ivhaib=dblalloc(nbfc*dfnbasis*nbex)
       else
        ivhaib=ivhaia
       endif
       if(lcc) then

c       if(freememspace.lt.nbfc*dfnbasis*nalx*2)then
        if(.false.)then
         write(iout,*) 'Haixx too big for RAM...'
      write(iout,"('Memalloc attempt for Haixx block (Mb):   ',14x,i6)")
     $(8*nbfc*dfnbasis*nlm)/(1024*1024)
         i=dblalloc(dfnbasis*nlm*nbfc) !haixx

         call dbldealloc(i)
        else

         write(iout,"('Memalloc attempt for Haixx (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*nlen)/(1024*1024)
         i=dblalloc(dfnbasis*nb*nlen) !haixx

         call dbldealloc(i)
        endif
       endif
      write(iout,*)
      write(iout,*) 'Calculating Coulomb integrals in blocks...'


C calc. # of blocks
      if(mod(nalx,nlen).ne.0)then
C if nal is not a multiple of nel
       nbbl = mod(nalx,nlen)
       nbb = 1+nalx/nlen
      else
C if nal is a multiple of nel
       nbb = nalx/nlen
       nbbl = nlen
      endif

      call gintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vifaa,vifab,cmoa,cmob,nbfcn,nalc,nalx,ncore,dfnbasis,spctostr,
     $cprea,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,icore,
     $scftype,dcore(ivhaia),dcore(ivhaia),lcc,verbosity,
     $ifltln,
     $nlen,
     $nbbl,
     $nbb)
c
      open(95,status='unknown',file='MPVIFAA',form='unformatted')
      write(95)nalx
      close(95)
c

       if(scftype.ne.'rhf  ') then
C Construct beta list
        call gintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vifbb,vifba,cmob,cmoa,nbfcn,nbec,nbex,ncore,dfnbasis,spctostr,
     $cpreb,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,icore,
     $scftype,dcore(ivhaib),dcore(ivhaib),dcore(imem),lcc,verbosity)
       endif
      call timer
C (DelF12)^2 integrals
      call f12init(nmboys)
      write(iout,*)
      write(iout,*) 'Calculating (DelF12)^2 integrals...'

      call fintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf122aa,cf122bb,cf122ab,cf122ba,boysdelf122,'delf122 ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
      call timer
C F12/r12 integrals
      write(iout,*)
      write(iout,*) 'Calculating F12/r12 integrals...'
      call fintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf12aa,cf12bb,cf12ab,cf12ba,boysf12r12,'f12r12  ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
C Calculate (DelF12)^2 and F12/r12 contributions to pair enegies
      call vb1_bl(cf122aa,cf122bb,cf122ab,cf122ba,cf12aa,cf12bb,cf12ab,
     $cf12ba,nalc,nbec,dfnbasis,epaa,epbb,epab,scftype,vifaa,vifbb,
     $vifab,vifba,ncore,nbfcn,dcore(ihai),dcore(ihai+dfnbasis),
     $dcore(ihai+2*dfnbasis),dcore(ihai+3*dfnbasis),ifltln)
      call timer
C F12^2 integrals
      write(iout,*)
      write(iout,*) 'Calculating F12^2 integrals...'
      call fintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf122aa,cf122bb,cf122ab,cf122ba,boysf122,'f122    ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
C Calculate F12^2 contribution epair_stf122epair_stf122to pair enegies
      iscr  =dblalloc(nbfcn*nbfcn)
      if12ij=dblalloc(nbfcn*nbfcn)
       if(scftype.eq.'rhf  ') then
C Closed-shell
        call epair_stf122_bl(nalc,nbfcn,ncore,dfnbasis,focka,hpja,vifaa,
     $cf122aa,epaa,dcore(if12ij),dcore(iscr),ifltln)
       else
C alpha-alpha
        call epair_ssf122(nalc,nbfcn,ncore,dfnbasis,focka,hpja,vifaa,
     $cf122aa,epaa,dcore(if12ij),dcore(iscr))
C beta-beta
        call epair_ssf122(nbec,nbfcn,ncore,dfnbasis,fockb,hpjb,vifbb,
     $cf122bb,epbb,dcore(if12ij),dcore(iscr))
C alpha-beta
        call epair_osf122(nalc,nbec,nbfcn,ncore,dfnbasis,focka,fockb,
     $hpja,hpjb,vifaa,vifbb,vifab,vifba,cf122aa,cf122bb,cf122ab,cf122ba,
     $epab,dcore(if12ij))
       endif
      call dbldealloc(iscr)
      call timer
C F12 integrals
      write(iout,*)
      write(iout,*) 'Calculating F12 integrals...'
      call fintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf12aa,cf12bb,cf12ab,cf12ba,boysf12,'f12     ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
      call timer
C
      call dbldealloc(ihai)
C
      return
      end
C
************************************************************************
      subroutine mp2int(nbf,nalc,nbec,ncore,nb,dfnbasis,nbfc,nbfcn,cmoa,
     $cmob,epaa,epbb,epab,focka,fockb,hpja,hpjb,dcore,imem,natoms,iout,
     $imem1,maxcor,nangmax,ncontrmax,nprimmax,ncartmax,cartg,nsphermax,
     $nmboys,itol,nbfshmax,nang,ncontr,nprim,gexp,gcoef,coord,ctostr,cf,
     $boysval,indarr,gcn,pre,nshrange,thad,thcf2,scoord,rqqij,rqqkl,
     $hrec,nbset,spctostr,dfipre,icore,tcdfint,dfrqq,spre,natrange,
     $dfipra,logkc,kp,gck,cprea,cpreb,tcf12int,vifaa,vifbb,vifab,cf12aa,
     $cf12bb,cf12ab,cf122aa,cf122bb,cf122ab,vifba,cf12ba,cf122ba,
     $scrfile2,scftype,lcc,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,
     $verbosity,ifltln)
************************************************************************
c Calculate and transform two- and three-center integrals
************************************************************************
      implicit none
      integer nbf,nalc,nbec,ncore,dfnbasis,imem,natoms,iout,imem1,i,nb
      integer dblalloc,intalloc,icore(*),nbfc,nbfcn,scrfile2,nn,p,q,r,s
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp,maxcor
      integer nshrange(2,0:nangmax,natoms,nbset),ivhaia,ivhaib,ihai
      integer natrange(2,natoms,nbset),if12ij,iscr,nalx,nbex,verbosity
      integer ivifaa,j,k,ivifaa2,ifltln,freememspace,nlm
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2),dfipra(natoms)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cprea(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),dfrqq,spre,dfipre(natoms,0:nangmax)
      real*8 vifaa(dfnbasis,nbfcn,nbfcn),vifbb(dfnbasis,nbfcn,nbfcn)
      real*8 vifab(dfnbasis,nbfcn,nbfcn),vifba(dfnbasis,nbfcn,nbfcn)
      real*8 cpreb(natoms,0:nangmax),tcf12int(dfnbasis,dfnbasis)
      real*8 cf12aa(nbfcn,dfnbasis,*),cf12bb(nbfcn,dfnbasis,*)
      real*8 cf12ab(nbfcn,dfnbasis,*),cf12ba(nbfcn,dfnbasis,*)
      real*8 cf122aa(nbfcn,dfnbasis,*),cf122bb(nbfcn,dfnbasis,*)
      real*8 cf122ab(nbfcn,dfnbasis,*),cf122ba(nbfcn,dfnbasis,*)
      real*8 epaa(*),epbb(*),epab(*)
      real*8 hpja(nbfcn,nbfcn),hpjb(nbfcn,nbfcn)
      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn)
      real*8 cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba
      character(len=5) scftype
      logical cartg,lcc,bver
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
      interface
      subroutine boys
      end
      subroutine boysf12
      end
      subroutine boysf122
      end
      subroutine boysdelf122
      end
      subroutine boysf12r12
      end
      end interface
C MO coefficient prescreening
       if(lcc) then
        nalx=nb
        nbex=nb
       else
        nalx=nalc
        nbex=nbec
       endif
      cprea=0.d0
      call dfcpre(cmoa(1,ncore+1),cprea,natoms,nangmax,ncontrmax,
     $nsphermax,cartg,indarr,nang,ncontr,nbfc,nalx,.true.)
       if(scftype.ne.'rhf  ') then
        cpreb=0.d0
        call dfcpre(cmob(1,ncore+1),cpreb,natoms,nangmax,ncontrmax,
     $nsphermax,cartg,indarr,nang,ncontr,nbfc,nbex,.true.)
       endif
C Two-center Coulomb integrals
      call df2int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,nmboys,
     $dcore,imem,0,0.00001d0*itol,nshrange(1,0,1,3),iout,imem1,maxcor,
     $thad,thcf2,scoord,rqqij,rqqkl,0,1,nang(1,3),ncontr(0,1,3),
     $nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfrqq,dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),tcdfint,dfipre,i,i,.true.,
     $scrfile2,tcdfint,spctostr,1,dfipra,tcmax,.false.,tcdfint,
     $0.d0,boys)
      iisyev=0
      call dpptrf('L',dfnbasis,tcdfint,isyev)
       if(isyev.ne.0) then
        write(iout,*) 'Fatal error at the Cholesky decomposition!'
        call mrccend(1)
       endif
C Three-center Coulomb integrals + first half transformation
      ihai=dblalloc(dfnbasis*max(4,nalc,nbec)*nbfc)

      write(iout,"('Bottleneck Haia size (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*nalx)/(1024*1024)

      freememspace = maxcor-imem+imem1
      bver=.false.
      if(freememspace.lt.nbfc*dfnbasis*nalx*2)then
       bver=.true.
       write(iout,*) 'Haia too big for RAM...'
C TODO calculate this
       nlm = 50
      write(iout,"('Memalloc attempt for Haia block (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*nlm)/(1024*1024)
       ivhaia=dblalloc(nbfc*dfnbasis*nlm)
      else
       ivhaia=dblalloc(nbfc*dfnbasis*nalx)
      endif

       if(scftype.ne.'rhf  ') then
        ivhaib=dblalloc(nbfc*dfnbasis*nbex)
       else
        ivhaib=ivhaia
       endif
       if(lcc) then

        if(freememspace.lt.nbfc*dfnbasis*nalx*2)then
         write(iout,*) 'Haixx too big for RAM...'
      write(iout,"('Memalloc attempt for Haixx block (Mb):   ',14x,i6)")
     $(8*nbfc*dfnbasis*nlm)/(1024*1024)
         i=dblalloc(dfnbasis*nlm*nbfc) !haixx

         call dbldealloc(i)
        else

         write(iout,"('Memalloc attempt for Haixx (Mb):    ',14x,i6)")
     $(8*nbfc*dfnbasis*nb)/(1024*1024)
         i=dblalloc(dfnbasis*nb*nbfc) !haixx

         call dbldealloc(i)
        endif
       endif
      write(iout,*)
      write(iout,*) 'Calculating Coulomb integrals...'

      call gintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vifaa,vifab,cmoa,cmob,nbfcn,nalc,nalx,ncore,dfnbasis,spctostr,
     $cprea,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,icore,
     $scftype,dcore(ivhaia),dcore(ivhaia),dcore(imem),lcc,verbosity)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Bence modification
c     open(ivifaa,status='unknown',file='MPVIFAA',form='unformatted')
c     write(ivifaa)nalx
c     close(ivifaa)

c     open(ivifaa2,status='unknown',file='MPVIFAA2',
c    $     access='direct',recl=dfnbasis*nbfcn*ifltln)
c     do k=1,nalx 
c      write(ivifaa2,rec=k)vifaa(1:dfnbasis,1:nbfcn,k)
c     enddo
c     close(ivifaa2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       if(scftype.ne.'rhf  ') then
C Construct beta list
        call gintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vifbb,vifba,cmob,cmoa,nbfcn,nbec,nbex,ncore,dfnbasis,spctostr,
     $cpreb,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,icore,
     $scftype,dcore(ivhaib),dcore(ivhaib),dcore(imem),lcc,verbosity)
       endif
      call timer
C (DelF12)^2 integrals
      call f12init(nmboys)
      write(iout,*)
      write(iout,*) 'Calculating (DelF12)^2 integrals...'

      call fintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf122aa,cf122bb,cf122ab,cf122ba,boysdelf122,'delf122 ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
      call timer
C F12/r12 integrals
      write(iout,*)
      write(iout,*) 'Calculating F12/r12 integrals...'
      call fintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf12aa,cf12bb,cf12ab,cf12ba,boysf12r12,'f12r12  ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
C Calculate (DelF12)^2 and F12/r12 contributions to pair enegies
      call vb1(cf122aa,cf122bb,cf122ab,cf122ba,cf12aa,cf12bb,cf12ab,
     $cf12ba,nalc,nbec,dfnbasis,epaa,epbb,epab,scftype,vifaa,vifbb,
     $vifab,vifba,ncore,nbfcn,dcore(ihai),dcore(ihai+dfnbasis),
     $dcore(ihai+2*dfnbasis),dcore(ihai+3*dfnbasis))
      call timer
C F12^2 integrals
      write(iout,*)
      write(iout,*) 'Calculating F12^2 integrals...'
      call fintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf122aa,cf122bb,cf122ab,cf122ba,boysf122,'f122    ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
C Calculate F12^2 contribution epair_stf122epair_stf122to pair enegies
      iscr  =dblalloc(nbfcn*nbfcn)
      if12ij=dblalloc(nbfcn*nbfcn)
       if(scftype.eq.'rhf  ') then
C Closed-shell
        call epair_stf122(nalc,nbfcn,ncore,dfnbasis,focka,hpja,vifaa,
     $cf122aa,epaa,dcore(if12ij),dcore(iscr))
       else
C alpha-alpha
        call epair_ssf122(nalc,nbfcn,ncore,dfnbasis,focka,hpja,vifaa,
     $cf122aa,epaa,dcore(if12ij),dcore(iscr))
C beta-beta
        call epair_ssf122(nbec,nbfcn,ncore,dfnbasis,fockb,hpjb,vifbb,
     $cf122bb,epbb,dcore(if12ij),dcore(iscr))
C alpha-beta
        call epair_osf122(nalc,nbec,nbfcn,ncore,dfnbasis,focka,fockb,
     $hpja,hpjb,vifaa,vifbb,vifab,vifba,cf122aa,cf122bb,cf122ab,cf122ba,
     $epab,dcore(if12ij))
       endif
      call dbldealloc(iscr)
      call timer
C F12 integrals
      write(iout,*)
      write(iout,*) 'Calculating F12 integrals...'
      call fintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $dcore(ivhaia),dcore(ivhaib),cmoa,cmob,nbfcn,nalc,nbec,ncore,
     $nb,dfnbasis,spctostr,cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,
     $logkc,kp,gck,nangmin,icore,scftype,dcore(ihai),tcf12int,dfrqq,
     $cf12aa,cf12bb,cf12ab,cf12ba,boysf12,'f12     ',lcc,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,ifltln)
      call timer
C
      call dbldealloc(ihai)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      return
      end
C
************************************************************************
      subroutine gintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vifaa,vifab,cmoa,cmob,nbfcn,nalc,nalx,ncore,dfnbasis,spctostr,
     $cpre,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,icore,
     $scftype,hai,haix,lcc,verbosity,ifltln,
     $nlen, ! number of elements in a block
     $nbbl, ! number of elements in the last block
     $nblock ! number of blocks
     $)
************************************************************************
c Calculate and transform three-center Coulomb integrals
c BLOCK version
************************************************************************
      implicit none
      integer dfnbasis,imem,natoms,iout,imem1,maxcor,i,nalc,ncore,nalx
      integer icore(*),nbfc,nbfcn,nn,p,q,r,s,natrange(2,natoms,nbset)
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp,verbosity
      integer nshrange(2,0:nangmax,natoms,nbset),j,k,ifltln,nlen,nlen2
      integer jblock,nblock,nbbl
      integer ii,jj,kk
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cpre(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),spre,dfipre(natoms,0:nangmax)
      real*8 dfipra(natoms),hai(dfnbasis,nlen,nbfc)
      real*8 vifaa(nbfcn,dfnbasis,nlen),vifab(nbfcn,dfnbasis,nalx)
      real*8 haix(dfnbasis,nalc,nbfc)
      character(len=5) scftype
      logical cartg,lcc,whai
      integer*4 isyev
      equivalence(isyev,iisyev) !for intel
      interface
      subroutine boys
      end
      end interface
c three-center coulomb integrals + first half transformation
      open(95,status='unknown',file='MPVIFAA2',
     $     access='direct',recl=dfnbasis*nbfcn*ifltln)

c     open(111,status='unknown',file='MPHAI',
c    $     access='direct',recl=dfnbasis*nbfc*ifltln)

      write(iout,*)'full length: ',nalx
      write(iout,*)'# of df3int calls: ',nblock
      write(iout,*)'# of elements/block: ',nlen
      write(iout,*)'# of elements in the last block: ',nbbl
      write(iout,"('bottleneck mem usage (Mb):    ',14x,i6)")
     $(8*(nbfcn*dfnbasis*nlen+2*dfnbasis*nlen*nbfc))/(1024*1024)

      whai = .true.
      j = 1
      do jblock=1,nblock
       write(iout,*)'jblock: ',jblock
c     if(jblock.eq.nblock) nlen = 
      nlen2 = nlen
      if(jblock.eq.nblock) nlen2 = nbbl
       write(iout,*)'calling df3int...'
      call df3int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,i,nbfc,gcn,dcore(imem),pre,5,itol,
     $nbfshmax,nshrange,i,iout,dcore,dcore,dcore,dcore,'     ',imem1,
     $maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,0,1,nang(1,3),
     $ncontr(0,1,3),nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),dcore(imem),tcdfint,dfipre,
     $dcore,hai,cmoa(1,ncore+j),nlen2,dcore,dfnbasis*nlen2*nbfc,i,i,i,i,
     $i,
     $i,i,itol,0,spctostr,1,nlen2,cpre,.true.,.true.,0,1.d0,.false.,0,0,
     $i,.false.,0.d0,dcore,0,nlen2,i,i,i,i,i,i,dcore,i,i,i,i,7,dcore,1,
     $.false.,dcore,.true.,nlen2,natrange(1,1,3),i,i,i,i,i,i,i,i,i,
     $dcore,i,dcore,dcore,dcore,dcore,dcore,spre,dfipra,tcmax,logkc,kp,
     $gck,dcore,0,0,i,i,i,i,i,i,natrange,i,i,i,0,icore,0.d0,boys,1,
     $dcore,i,i,i,.false.,dcore,i,dcore,i,dcore,dcore,i,nangmin,
     $verbosity.ge.3,i,i,dcore,dcore,dcore,dcore,dcore,dcore,dcore,
     $dcore,dcore,dcore,dcore,'off     ',i,.false.,dcore,.false.,
     $.false.,dcore,dcore,.false.)
c fitting
      nn=nlen2*nbfc
      call dtptrs('l','n','n',dfnbasis,nn,tcdfint,hai,dfnbasis,isyev)
       if(isyev.ne.0) then
        write(iout,*)
     $'fatal error at the fitting of the coulomb integrals!'
        call mrccend(1)
       endif
c second half transformation
c     call dgemm('t','t',nbfcn,dfnbasis*nalx,nbfc,1.d0,cmoa,nbfc,hai,
c    $dfnbasis*nalx,0.d0,vifaa,nbfcn)

      call dgemm('t','t',nbfcn,dfnbasis*nlen2,nbfc,1.d0,cmoa,nbfc,hai,
     $dfnbasis*nlen2,0.d0,vifaa,nbfcn)

      do ii=1,nlen2
       write(95,rec=(jblock-1)*nlen+ii)vifaa(1:nbfcn,1:dfnbasis,ii)
      enddo

C       if(whai)then
C        do ii=1,nlen2
CC Write hai only until nalc
C         if((jblock-1)*nlen+ii.le.nalc)then
C          write(111,rec=(jblock-1)*nlen+ii)hai(1:dfnbasis,ii,1:nbfc) ! write hai to disk until nalc
C         endif
C        enddo
CC Write hai only until nalc
C        if(jblock*nlen.gt.nalc)then
C         whai = .false.
C        endif 
C       endif ! if write hai

      j = j + nlen2
      enddo ! jblock

      close(95)
      
      if(scftype.ne.'rhf  ')
     $call dgemm('t','t',nbfcn,dfnbasis*nalx,nbfc,1.d0,cmob,nbfc,hai,
     $dfnbasis*nalx,0.d0,vifab,nbfcn)
c reorder hai in the case of cc
      if(lcc) then

      if(nlen.lt.nalc)write(iout,*)"nlen.lt.nalc(gintf_bl) BUG"
      write(iout,*)'last df3int until nalc...'
      call df3int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,i,nbfc,gcn,dcore(imem),pre,5,itol,
     $nbfshmax,nshrange,i,iout,dcore,dcore,dcore,dcore,'     ',imem1,
     $maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,0,1,nang(1,3),
     $ncontr(0,1,3),nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),dcore(imem),tcdfint,dfipre,
     $dcore,haix,cmoa(1,ncore+1),nalc,dcore,dfnbasis*nalc*nbfc,i,i,i,i,
     $i,
     $i,i,itol,0,spctostr,1,nalc,cpre,.true.,.true.,0,1.d0,.false.,0,0,
     $i,.false.,0.d0,dcore,0,nalc,i,i,i,i,i,i,dcore,i,i,i,i,7,dcore,1,
     $.false.,dcore,.true.,nalc,natrange(1,1,3),i,i,i,i,i,i,i,i,i,
     $dcore,i,dcore,dcore,dcore,dcore,dcore,spre,dfipra,tcmax,logkc,kp,
     $gck,dcore,0,0,i,i,i,i,i,i,natrange,i,i,i,0,icore,0.d0,boys,1,
     $dcore,i,i,i,.false.,dcore,i,dcore,i,dcore,dcore,i,nangmin,
     $verbosity.ge.3,i,i,dcore,dcore,dcore,dcore,dcore,dcore,dcore,
     $dcore,dcore,dcore,dcore,'off     ',i,.false.,dcore,.false.,
     $.false.,dcore,dcore,.false.)
c fitting
      call dtptrs('l','n','n',dfnbasis,nalc*nbfc,tcdfint,haix,dfnbasis,
     $     isyev)

C       do jblock=1,nblock
C        do ii=1,nlen2
C         if((jblock-1)*nlen+ii.le.nalc)then
C          read(111,rec=(jblock-1)*nlen+ii)haixx(1:dfnbasis,ii,1:nbfc)
C          do kk=1,dfnbasis
C           do jj=1,nbfc
C            haix(kk,(jblock-1)*nlen+ii,jj) = haixx(kk,ii,jj)
C           enddo ! loop on hai elements 1
C          enddo ! loop on hai elements 2
C         endif ! write if inside nalc block
C        enddo ! loop on elements in a block
C       enddo ! loop on blocks
      endif ! if lcc

c     close(111) ! close hai file on disk
C
      return
      end
C
************************************************************************
      subroutine gintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vifaa,vifab,cmoa,cmob,nbfcn,nalc,nalx,ncore,dfnbasis,spctostr,
     $cpre,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,icore,
     $scftype,hai,haix,haixx,lcc,verbosity)
************************************************************************
c Calculate and transform three-center Coulomb integrals
************************************************************************
      implicit none
      integer dfnbasis,imem,natoms,iout,imem1,maxcor,i,nalc,ncore,nalx
      integer icore(*),nbfc,nbfcn,nn,p,q,r,s,natrange(2,natoms,nbset)
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp,verbosity
      integer nshrange(2,0:nangmax,natoms,nbset),ivifaa,j,k
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cpre(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),spre,dfipre(natoms,0:nangmax)
      real*8 dfipra(natoms),hai(dfnbasis,nalx,nbfc)
      real*8 vifaa(nbfcn,dfnbasis,nalx),vifab(nbfcn,dfnbasis,nalx)
      real*8 haix(dfnbasis,nalc,nbfc),haixx(dfnbasis,nalx,nbfc)
      character(len=5) scftype
      logical cartg,lcc
      integer*4 isyev
      equivalence(isyev,iisyev) !for intel
      interface
      subroutine boys
      end
      end interface
c three-center coulomb integrals + first half transformation
      call df3int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,i,nbfc,gcn,dcore(imem),pre,5,itol,
     $nbfshmax,nshrange,i,iout,dcore,dcore,dcore,dcore,'     ',imem1,
     $maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,0,1,nang(1,3),
     $ncontr(0,1,3),nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),dcore(imem),tcdfint,dfipre,
     $dcore,hai,cmoa(1,ncore+1),nalx,dcore,dfnbasis*nalx*nbfc,i,i,i,i,i,
     $i,i,itol,0,spctostr,1,nalx,cpre,.true.,.true.,0,1.d0,.false.,0,0,
     $i,.false.,0.d0,dcore,0,nalx,i,i,i,i,i,i,dcore,i,i,i,i,7,dcore,1,
     $.false.,dcore,.true.,nalx,natrange(1,1,3),i,i,i,i,i,i,i,i,i,
     $dcore,i,dcore,dcore,dcore,dcore,dcore,spre,dfipra,tcmax,logkc,kp,
     $gck,dcore,0,0,i,i,i,i,i,i,natrange,i,i,i,0,icore,0.d0,boys,1,
     $dcore,i,i,i,.false.,dcore,i,dcore,i,dcore,dcore,i,nangmin,
     $verbosity.ge.3,i,i,dcore,dcore,dcore,dcore,dcore,dcore,dcore,
     $dcore,dcore,dcore,dcore,'off     ',i,.false.,dcore,.false.,
     $.false.,dcore,dcore,.false.)
c fitting
      nn=nalx*nbfc
      call dtptrs('l','n','n',dfnbasis,nn,tcdfint,hai,dfnbasis,isyev)
       if(isyev.ne.0) then
        write(iout,*)
     $'fatal error at the fitting of the coulomb integrals!'
        call mrccend(1)
       endif
c second half transformation
      call dgemm('t','t',nbfcn,dfnbasis*nalx,nbfc,1.d0,cmoa,nbfc,hai,
     $dfnbasis*nalx,0.d0,vifaa,nbfcn)
       if(scftype.ne.'rhf  ')
     $call dgemm('t','t',nbfcn,dfnbasis*nalx,nbfc,1.d0,cmob,nbfc,hai,
     $dfnbasis*nalx,0.d0,vifab,nbfcn)
c reorder hai in the case of cc
       if(lcc) then
        call dcopy(dfnbasis*nalx*nbfc,hai,1,haixx,1)
         do i=1,nalc
          haix(:,i,:)=haixx(:,i,:)
         enddo
       endif
c
      return
      end
C
************************************************************************
      subroutine fintf_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vhaia,vhaib,cmoa,cmob,nbfcn,nalc,nbec,ncore,nb,dfnbasis,spctostr,
     $cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,
     $icore,scftype,hai,tcf12int,dfrqq,cf12aa,cf12bb,cf12ab,cf12ba,boys,
     $intype,lcc,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,
     $ifltln)
************************************************************************
c Calculate and transform three-center Coulomb integrals
c BLOCK version
************************************************************************
      implicit none
      integer dfnbasis,imem,natoms,iout,imem1,maxcor,i,nb,verbosity
      integer icore(*),nbfc,nbfcn,nalc,nbec,ncore
      integer ii,jj,kk
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp
      integer nshrange(2,0:nangmax,natoms,nbset)
      integer natrange(2,natoms,nbset)
      integer nlm,ifltln,ibl,freememspace
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cprea(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),dfrqq,spre,dfipre(natoms,0:nangmax)
      real*8 dfipra(natoms)
C Mem. heavy tensors
C     real*8 hai(dfnbasis,nalc,nbfc)
      real*8 hai(dfnbasis,nalc,nbfc)
      real*8 cf12aa(nbfcn,dfnbasis,*),vhaia(dfnbasis,nalc,nbfc)
      real*8 cf12bb(nbfcn,dfnbasis,*),vhaib(dfnbasis,nalc,nbfc)
      real*8 cf12ab(nbfcn,dfnbasis,*),cf12ba(nbfcn,dfnbasis,*)
C
      real*8 tcf12int(dfnbasis,dfnbasis),cpreb(natoms,0:nangmax)
      real*8 cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba
      character(len=5) scftype
      character(len=8) intype
      logical cartg,lcc
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
      interface
      subroutine boys
      end
      end interface
C Two-center integrals
      write(iout,*)'2center int for block: ',ibl

      freememspace = maxcor-imem+imem1
      write(iout,"('M.left (Mb):',14x,i6)")(8*freememspace)/(1024*1024)


      i = 1

      call df2int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,nmboys,
     $dcore,imem,0,0.00001d0*itol,nshrange(1,0,1,3),iout,imem1,maxcor,
     $thad,thcf2,scoord,rqqij,rqqkl,0,1,nang(1,3),ncontr(0,1,3),
     $nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfrqq,dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),tcf12int,dfipre,i,i,.true.,
     $i,tcf12int,spctostr,1,dfipra,tcmax,.true.,tcf12int,1.d0,boys)
      call dtptrs('L','N','N',dfnbasis,dfnbasis,tcdfint,tcf12int,
     $dfnbasis,isyev)
       if(isyev.ne.0) then
        write(iout,*) 'Fatal error at the fitting of the F12 integrals!'
        call mrccend(1)
       endif
C Construct aa fitting coefficients

      call fintfa_bl(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $nbfcn,nalc,nbec,ncore,nb,dfnbasis,spctostr,natrange,spre,nbset,
     $dfipra,tcmax,logkc,kp,gck,nangmin,icore,hai,boys,tcf12int,vhaia,
     $cf12aa,cf12ab,cmoa,cmob,cprea,intype,scftype,cf12r12aa,cf12r12ab,
     $lcc,verbosity,ifltln)


       if(scftype.eq.'rhf  ') return
C Construct bb fitting coefficients
C TODO meg kell csinalni az OS blokkositast
      call fintfa(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $nbfcn,nbec,nalc,ncore,nb,dfnbasis,spctostr,natrange,spre,nbset,
     $dfipra,tcmax,logkc,kp,gck,nangmin,icore,hai,boys,tcf12int,vhaib,
     $cf12bb,cf12ba,cmob,cmoa,cpreb,intype,scftype,cf12r12bb,cf12r12ba,
     $lcc,verbosity,ifltln)
C
      return
      end
C
************************************************************************
      subroutine fintf(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $vhaia,vhaib,cmoa,cmob,nbfcn,nalc,nbec,ncore,nb,dfnbasis,spctostr,
     $cprea,cpreb,natrange,spre,nbset,dfipra,tcmax,logkc,kp,gck,nangmin,
     $icore,scftype,hai,tcf12int,dfrqq,cf12aa,cf12bb,cf12ab,cf12ba,boys,
     $intype,lcc,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,verbosity,
     $ifltln)
************************************************************************
c Calculate and transform three-center Coulomb integrals
************************************************************************
      implicit none
      integer dfnbasis,imem,natoms,iout,imem1,maxcor,i,nb,verbosity
      integer icore(*),nbfc,nbfcn,nalc,nbec,ncore
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp
      integer nshrange(2,0:nangmax,natoms,nbset)
      integer natrange(2,natoms,nbset)
      integer ifltln
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cprea(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),dfrqq,spre,dfipre(natoms,0:nangmax)
      real*8 dfipra(natoms),hai(dfnbasis,nalc,nbfc)
      real*8 cf12aa(nbfcn,dfnbasis,*),vhaia(dfnbasis,nalc,nbfc)
      real*8 cf12bb(nbfcn,dfnbasis,*),vhaib(dfnbasis,nalc,nbfc)
      real*8 cf12ab(nbfcn,dfnbasis,*),cf12ba(nbfcn,dfnbasis,*)
      real*8 tcf12int(dfnbasis,dfnbasis),cpreb(natoms,0:nangmax)
      real*8 cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba
      character(len=5) scftype
      character(len=8) intype
      logical cartg,lcc
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
      interface
      subroutine boys
      end
      end interface
C Two-center integrals
      call df2int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,nmboys,
     $dcore,imem,0,0.00001d0*itol,nshrange(1,0,1,3),iout,imem1,maxcor,
     $thad,thcf2,scoord,rqqij,rqqkl,0,1,nang(1,3),ncontr(0,1,3),
     $nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfrqq,dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),tcf12int,dfipre,i,i,.true.,
     $i,tcf12int,spctostr,1,dfipra,tcmax,.true.,tcf12int,1.d0,boys)
      call dtptrs('L','N','N',dfnbasis,dfnbasis,tcdfint,tcf12int,
     $dfnbasis,isyev)
       if(isyev.ne.0) then
        write(iout,*) 'Fatal error at the fitting of the F12 integrals!'
        call mrccend(1)
       endif
C Construct aa fitting coefficients
      call fintfa(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $nbfcn,nalc,nbec,ncore,nb,dfnbasis,spctostr,natrange,spre,nbset,
     $dfipra,tcmax,logkc,kp,gck,nangmin,icore,hai,boys,tcf12int,vhaia,
     $cf12aa,cf12ab,cmoa,cmob,cprea,intype,scftype,cf12r12aa,cf12r12ab,
     $lcc,verbosity,ifltln)
       if(scftype.eq.'rhf  ') return
C Construct bb fitting coefficients
      call fintfa(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $nbfcn,nbec,nalc,ncore,nb,dfnbasis,spctostr,natrange,spre,nbset,
     $dfipra,tcmax,logkc,kp,gck,nangmin,icore,hai,boys,tcf12int,vhaib,
     $cf12bb,cf12ba,cmob,cmoa,cpreb,intype,scftype,cf12r12bb,cf12r12ba,
     $lcc,verbosity,ifltln)
C
      return
      end
C
************************************************************************
      subroutine fintfa_bl(natoms,nangmax,ncontrmax,nprimmax,nang,
     $ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $nbfcn,nalc,nbec,ncore,nb,dfnbasis,spctostr,natrange,spre,nbset,
     $dfipra,tcmax,logkc,kp,gck,nangmin,icore,hai,boys,tcf12int,vhaia,
     $cf12aa,cf12ab,cmoa,cmob,cpre,intype,scftype,cf12r12aa,cf12r12ab,
     $lcc,verbosity,ifltln)
************************************************************************
c Calculate and transform three-center Coulomb integrals
************************************************************************
      implicit none
      integer dfnbasis,imem,natoms,iout,imem1,maxcor,i,nb,verbosity
      integer icore(*),nbfc,nbfcn,nn,nalc,nbec,ncore
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp
      integer nshrange(2,0:nangmax,natoms,nbset)
      integer natrange(2,natoms,nbset)
      integer ifltln,k
      logical fits_in_memory
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cpre(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),spre,dfipre(natoms,0:nangmax)
      real*8 dfipra(natoms),hai(dfnbasis,nalc,nbfc)
      real*8 vhaia(dfnbasis,nalc,nbfc),cf12aa(nbfcn,dfnbasis,nalc)
      real*8 tcf12int(dfnbasis,dfnbasis),cf12ab(nbfcn,dfnbasis,nalc)
      real*8 cf12r12aa,cf12r12ab
      logical cartg,lcc
      character(len=5) scftype
      character(len=8) intype
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
      interface
      subroutine boys
      end
      end interface
C Three-center integrals + first half transformation
      call df3int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,i,nbfc,gcn,dcore(imem),pre,5,itol,
     $nbfshmax,nshrange,i,iout,dcore,dcore,dcore,dcore,'     ',imem1,
     $maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,0,1,nang(1,3),
     $ncontr(0,1,3),nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),dcore(imem),tcdfint,dfipre,
     $dcore,hai,cmoa(1,ncore+1),nalc,dcore,dfnbasis*nalc*nbfc,i,i,i,i,i,
     $i,i,itol,0,spctostr,1,nalc,cpre,.true.,.true.,0,1.d0,.false.,0,0,
     $i,.false.,0.d0,dcore,0,nalc,i,i,i,i,i,i,dcore,i,i,i,i,7,dcore,1,
     $.false.,dcore,.true.,nalc,natrange(1,1,3),i,i,i,i,i,i,i,i,i,
     $dcore,i,dcore,dcore,dcore,dcore,dcore,spre,dfipra,tcmax,logkc,kp,
     $gck,dcore,0,0,i,i,i,i,i,i,natrange,i,i,i,0,icore,1.d0,boys,1,
     $dcore,i,i,i,.false.,dcore,i,dcore,i,dcore,dcore,i,nangmin,
     $verbosity.ge.3,i,i,dcore,dcore,dcore,dcore,dcore,dcore,dcore,
     $dcore,dcore,dcore,dcore,'off     ',i,.false.,dcore,.false.,
     $.false.,dcore,dcore,.false.)
C Fitting
      nn=nalc*nbfc
      call dgemm('t','n',dfnbasis,nn,dfnbasis,-0.5d0,tcf12int,dfnbasis,
     $vhaia,dfnbasis,1.d0,hai,dfnbasis)
      call dtptrs('L','N','N',dfnbasis,nn,tcdfint,hai,dfnbasis,isyev)
       if(isyev.ne.0) then
        write(iout,*) 'Fatal error at the fitting of the F12 integrals!'
        call mrccend(1)
       endif
C Second half transformation
       if(intype.eq.'f12     '.or.intype.eq.'f122    ') then
        call dgemm('t','t',nbfcn,dfnbasis*nalc,nbfc,1.d0,cmoa,nbfc,hai,
     $dfnbasis*nalc,0.d0,cf12aa,nbfcn)
C Write cf12aa if out-of-memory
      fits_in_memory = .false.
      if(.not.fits_in_memory)then

       open(109,status='unknown',file='MPCF12AA',
     $      access='direct',recl=dfnbasis*nbfcn*ifltln)
       do k=1,nalc 
        write(109,rec=k)cf12aa(1:nbfcn,1:dfnbasis,k)
       enddo
       close(109)
       
      endif
C
       else
        call dgemm('n','n',dfnbasis*nalc,nalc,nbfc,1.d0,hai,
     $dfnbasis*nalc,cmoa(1,ncore+1),nbfc,0.d0,cf12aa,dfnbasis*nalc)
         if(intype.eq.'f12r12  '.and.lcc)
     $  call dgemm('t','t',nb,dfnbasis*nalc,nbfc,1.d0,cmoa(1,ncore+1),
     $nbfc,hai,dfnbasis*nalc,0.d0,cf12r12aa,nb)
       endif
       if(scftype.eq.'rhf  ') return
       if(intype.eq.'f12     '.or.intype.eq.'f122    ') then
        call dgemm('t','t',nbfcn,dfnbasis*nalc,nbfc,1.d0,cmob,nbfc,hai,
     $dfnbasis*nalc,0.d0,cf12ab,nbfcn)
       else
        call dgemm('n','n',dfnbasis*nalc,nbec,nbfc,1.d0,hai,
     $dfnbasis*nalc,cmob(1,ncore+1),nbfc,0.d0,cf12ab,dfnbasis*nalc)
         if(intype.eq.'f12r12  '.and.lcc)
     $  call dgemm('t','t',nb,dfnbasis*nalc,nbfc,1.d0,cmob(1,ncore+1),
     $nbfc,hai,dfnbasis*nalc,0.d0,cf12r12ab,nb)
       endif
C
      return
      end
C
************************************************************************
      subroutine fintfa(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,
     $nprim,gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,nbfc,gcn,pre,itol,nbfshmax,nshrange,iout,
     $imem1,maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,tcdfint,dfipre,
     $nbfcn,nalc,nbec,ncore,nb,dfnbasis,spctostr,natrange,spre,nbset,
     $dfipra,tcmax,logkc,kp,gck,nangmin,icore,hai,boys,tcf12int,vhaia,
     $cf12aa,cf12ab,cmoa,cmob,cpre,intype,scftype,cf12r12aa,cf12r12ab,
     $lcc,verbosity,ifltln)
************************************************************************
c Calculate and transform three-center Coulomb integrals
************************************************************************
      implicit none
      integer dfnbasis,imem,natoms,iout,imem1,maxcor,i,nb,verbosity
      integer icore(*),nbfc,nbfcn,nn,nalc,nbec,ncore
      integer nangmax,ncontrmax,nprimmax,ncartmax,nsphermax,nmboys,nbset
      integer nbfshmax,nang(natoms,nbset),ncontr(0:nangmax,natoms,nbset)
      integer nprim(0:nangmax,natoms,nbset),nangmin(natoms,nbset),iisyev
      integer indarr(natoms,0:nangmax,ncontrmax,nsphermax,nbset),logkc
      integer gcn(2,ncontrmax,0:nangmax,natoms,nbset),thad,kp
      integer nshrange(2,0:nangmax,natoms,nbset)
      integer natrange(2,natoms,nbset)
      integer ifltln,k
      logical fits_in_memory
      real*8 dcore(*),tcdfint((dfnbasis+1)*dfnbasis/2)
      real*8 cmoa(nbfc,nbfcn),cmob(nbfc,nbfcn),cpre(natoms,0:nangmax)
      real*8 itol,gexp(nprimmax,0:nangmax,natoms,nbset),ctostr,cf,gck
      real*8 gcoef(nprimmax,ncontrmax,0:nangmax,natoms,nbset)
      real*8 rqqij,rqqkl,hrec,spctostr,boysval,pre,thcf2,scoord,tcmax
      real*8 coord(3,natoms),spre,dfipre(natoms,0:nangmax)
      real*8 dfipra(natoms),hai(dfnbasis,nalc,nbfc)
      real*8 vhaia(dfnbasis,nalc,nbfc),cf12aa(nbfcn,dfnbasis,nalc)
      real*8 tcf12int(dfnbasis,dfnbasis),cf12ab(nbfcn,dfnbasis,nalc)
      real*8 cf12r12aa,cf12r12ab
      logical cartg,lcc
      character(len=5) scftype
      character(len=8) intype
      integer*4 isyev
      equivalence(isyev,iisyev) !For Intel
      interface
      subroutine boys
      end
      end interface
C Three-center integrals + first half transformation
      call df3int(natoms,nangmax,ncontrmax,nprimmax,nang,ncontr,nprim,
     $gexp,gcoef,coord,ncartmax,ctostr,cartg,nsphermax,cf,boysval,
     $nmboys,dcore,imem,indarr,i,nbfc,gcn,dcore(imem),pre,5,itol,
     $nbfshmax,nshrange,i,iout,dcore,dcore,dcore,dcore,'     ',imem1,
     $maxcor,thad,thcf2,scoord,rqqij,rqqkl,hrec,0,1,nang(1,3),
     $ncontr(0,1,3),nprim(0,1,3),gexp(1,0,1,3),gcn(1,1,0,1,3),dfnbasis,
     $gcoef(1,1,0,1,3),indarr(1,0,1,1,3),dcore(imem),tcdfint,dfipre,
     $dcore,hai,cmoa(1,ncore+1),nalc,dcore,dfnbasis*nalc*nbfc,i,i,i,i,i,
     $i,i,itol,0,spctostr,1,nalc,cpre,.true.,.true.,0,1.d0,.false.,0,0,
     $i,.false.,0.d0,dcore,0,nalc,i,i,i,i,i,i,dcore,i,i,i,i,7,dcore,1,
     $.false.,dcore,.true.,nalc,natrange(1,1,3),i,i,i,i,i,i,i,i,i,
     $dcore,i,dcore,dcore,dcore,dcore,dcore,spre,dfipra,tcmax,logkc,kp,
     $gck,dcore,0,0,i,i,i,i,i,i,natrange,i,i,i,0,icore,1.d0,boys,1,
     $dcore,i,i,i,.false.,dcore,i,dcore,i,dcore,dcore,i,nangmin,
     $verbosity.ge.3,i,i,dcore,dcore,dcore,dcore,dcore,dcore,dcore,
     $dcore,dcore,dcore,dcore,'off     ',i,.false.,dcore,.false.,
     $.false.,dcore,dcore,.false.)
C Fitting
      nn=nalc*nbfc
      call dgemm('t','n',dfnbasis,nn,dfnbasis,-0.5d0,tcf12int,dfnbasis,
     $vhaia,dfnbasis,1.d0,hai,dfnbasis)
      call dtptrs('L','N','N',dfnbasis,nn,tcdfint,hai,dfnbasis,isyev)
       if(isyev.ne.0) then
        write(iout,*) 'Fatal error at the fitting of the F12 integrals!'
        call mrccend(1)
       endif
C Second half transformation
       if(intype.eq.'f12     '.or.intype.eq.'f122    ') then
        call dgemm('t','t',nbfcn,dfnbasis*nalc,nbfc,1.d0,cmoa,nbfc,hai,
     $dfnbasis*nalc,0.d0,cf12aa,nbfcn)
C Write cf12aa if out-of-memory
      fits_in_memory = .false.
      if(.not.fits_in_memory)then

       open(109,status='unknown',file='MPCF12AA',
     $      access='direct',recl=dfnbasis*nbfcn*ifltln)
       do k=1,nalc 
        write(109,rec=k)cf12aa(1:nbfcn,1:dfnbasis,k)
       enddo
       close(109)
       
      endif
C
       else
        call dgemm('n','n',dfnbasis*nalc,nalc,nbfc,1.d0,hai,
     $dfnbasis*nalc,cmoa(1,ncore+1),nbfc,0.d0,cf12aa,dfnbasis*nalc)
         if(intype.eq.'f12r12  '.and.lcc)
     $  call dgemm('t','t',nb,dfnbasis*nalc,nbfc,1.d0,cmoa(1,ncore+1),
     $nbfc,hai,dfnbasis*nalc,0.d0,cf12r12aa,nb)
       endif
       if(scftype.eq.'rhf  ') return
       if(intype.eq.'f12     '.or.intype.eq.'f122    ') then
        call dgemm('t','t',nbfcn,dfnbasis*nalc,nbfc,1.d0,cmob,nbfc,hai,
     $dfnbasis*nalc,0.d0,cf12ab,nbfcn)
       else
        call dgemm('n','n',dfnbasis*nalc,nbec,nbfc,1.d0,hai,
     $dfnbasis*nalc,cmob(1,ncore+1),nbfc,0.d0,cf12ab,dfnbasis*nalc)
         if(intype.eq.'f12r12  '.and.lcc)
     $  call dgemm('t','t',nb,dfnbasis*nalc,nbfc,1.d0,cmob(1,ncore+1),
     $nbfc,hai,dfnbasis*nalc,0.d0,cf12r12ab,nb)
       endif
C
      return
      end
C
************************************************************************
      subroutine vb1_bl(cdelf122aa,cdelf122bb,cdelf122ab,cdelf122ba,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,nalc,nbec,dfnbasis,epaa,
     $epbb,epab,scftype,vifaa,vifbb,vifab,vifba,ncore,nbfcn,vifii,vifjj,
     $vifij,vifji,ifltln)
************************************************************************
c Calculate intermediate B1 + intermediate V, F12/r12 contribution
c BLOCK version
************************************************************************
      implicit none
      integer nalc,nbec,dfnbasis,i,j,ij,nbfc,nbfcn,ncore,ifltln
      real*8 cdelf122aa(dfnbasis,nalc,nalc),epaa(*),epbb(*),ddot
      real*8 cdelf122bb(dfnbasis,nbec,nbec),vifaa(nbfcn,dfnbasis,2)
      real*8 cdelf122ab(dfnbasis,nalc,nbec),vifbb(nbfcn,dfnbasis,nbec)
      real*8 cf12r12aa(dfnbasis,nalc,nalc),cf12r12bb(dfnbasis,nbec,nbec)
      real*8 cf12r12ab(dfnbasis,nalc,nbec),vifii(dfnbasis)
      real*8 vifab(nbfcn,dfnbasis,nalc),vifjj(dfnbasis),vifij(dfnbasis)
      real*8 vifba(nbfcn,dfnbasis,nbec),cdelf122ba(dfnbasis,nbec,nalc)
      real*8 cf12r12ba(dfnbasis,nbec,nalc),epab(nalc,nbec)
      real*8 vifji(dfnbasis)
      character(len=5) scftype
C
      open(95,file='MPVIFAA2',access='direct',
     $recl=dfnbasis*nbfcn*ifltln,status='old')

       if(scftype.eq.'rhf  ') then
        ij=0
         do i=1,nalc

          read(95,rec=i)vifaa(1:nbfcn,1:dfnbasis,1)

          vifii=vifaa(ncore+i,1:dfnbasis,1)
           do j=1,i
            read(95,rec=j)vifaa(1:nbfcn,1:dfnbasis,2)
            vifjj=vifaa(ncore+j,1:dfnbasis,2)
            vifij=vifaa(ncore+i,1:dfnbasis,2)
            ij=ij+1
            epaa(ij)=
C Intermediate B, term 1 (= tau)
     $      7d0*(ddot(dfnbasis,vifii            ,1,cdelf122aa(1,j,j),1)+
     $           ddot(dfnbasis,cdelf122aa(1,i,i),1,vifjj            ,1))
     $          +ddot(dfnbasis,vifij            ,1,cdelf122aa(1,i,j),1)
     $          +ddot(dfnbasis,cdelf122aa(1,i,j),1,vifij            ,1)
C Intermediate V, F12/r12 contribution
     $+8d0*(5d0*(ddot(dfnbasis,vifii           ,1,cf12r12aa(1,j,j),1)+
     $           ddot(dfnbasis,cf12r12aa(1,i,i),1,vifjj           ,1))
     $          -ddot(dfnbasis,vifij           ,1,cf12r12aa(1,i,j),1) 
     $          -ddot(dfnbasis,cf12r12aa(1,i,j),1,vifij           ,1))
           enddo
         enddo
       else
C alpha-alpha
        ij=0
         do i=1,nalc
          vifii=vifaa(ncore+i,1:dfnbasis,i)
           do j=1,i-1
            vifjj=vifaa(ncore+j,1:dfnbasis,j)
            vifij=vifaa(ncore+i,1:dfnbasis,j)
            ij=ij+1
            epaa(ij)=
C Intermediate B, term 1 (= tau)
     $           ddot(dfnbasis,vifii            ,1,cdelf122aa(1,j,j),1)
     $          +ddot(dfnbasis,cdelf122aa(1,i,i),1,vifjj            ,1)
     $          -ddot(dfnbasis,vifij            ,1,cdelf122aa(1,i,j),1)
     $          -ddot(dfnbasis,cdelf122aa(1,i,j),1,vifij            ,1)
C Intermediate V, F12/r12 contribution
     $     +8d0*(ddot(dfnbasis,vifii           ,1,cf12r12aa(1,j,j),1) 
     $          +ddot(dfnbasis,cf12r12aa(1,i,i),1,vifjj           ,1) 
     $          -ddot(dfnbasis,vifij           ,1,cf12r12aa(1,i,j),1)
     $          -ddot(dfnbasis,cf12r12aa(1,i,j),1,vifij           ,1))
           enddo
         enddo
C beta-beta
        ij=0
         do i=1,nbec
          vifii=vifbb(ncore+i,1:dfnbasis,i)
           do j=1,i-1
            vifjj=vifbb(ncore+j,1:dfnbasis,j)
            vifij=vifbb(ncore+i,1:dfnbasis,j)
            ij=ij+1
            epbb(ij)=
C Intermediate B, term 1 (= tau)
     $           ddot(dfnbasis,vifii            ,1,cdelf122bb(1,j,j),1)
     $          +ddot(dfnbasis,cdelf122bb(1,i,i),1,vifjj            ,1)
     $          -ddot(dfnbasis,vifij            ,1,cdelf122bb(1,i,j),1)
     $          -ddot(dfnbasis,cdelf122bb(1,i,j),1,vifij            ,1)
C Intermediate V, F12/r12 contribution
     $     +8d0*(ddot(dfnbasis,vifii           ,1,cf12r12bb(1,j,j),1) 
     $          +ddot(dfnbasis,cf12r12bb(1,i,i),1,vifjj           ,1) 
     $          -ddot(dfnbasis,vifij           ,1,cf12r12bb(1,i,j),1)
     $          -ddot(dfnbasis,cf12r12bb(1,i,j),1,vifij           ,1))
           enddo
         enddo
C alpha-beta
         do i=1,nalc
          vifii=vifaa(ncore+i,1:dfnbasis,i)
           do j=1,nbec
            vifjj=vifbb(ncore+j,1:dfnbasis,j)
            vifij=vifba(ncore+i,1:dfnbasis,j)
            vifji=vifab(ncore+j,1:dfnbasis,i)
            epab(i,j)=
C Intermediate B, term 1 (= tau)
     $ 5d0/32d0*(ddot(dfnbasis,vifii            ,1,cdelf122bb(1,j,j),1)
     $          +ddot(dfnbasis,cdelf122aa(1,i,i),1,vifjj            ,1))
     $+3d0/32d0*(ddot(dfnbasis,vifij            ,1,cdelf122ab(1,i,j),1)
     $          +ddot(dfnbasis,cdelf122ba(1,j,i),1,vifji            ,1))
C Intermediate V, F12/r12 contribution
     $ +3d0/4d0*(ddot(dfnbasis,vifii           ,1,cf12r12bb(1,j,j),1)
     $          +ddot(dfnbasis,cf12r12aa(1,i,i),1,vifjj           ,1))
     $ +1d0/4d0*(ddot(dfnbasis,vifij           ,1,cf12r12ab(1,i,j),1)
     $          +ddot(dfnbasis,cf12r12ba(1,j,i),1,vifji           ,1))
           enddo
         enddo
       endif
C
      close(95)
      return
      end
C
************************************************************************
      subroutine vb1(cdelf122aa,cdelf122bb,cdelf122ab,cdelf122ba,
     $cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,nalc,nbec,dfnbasis,epaa,
     $epbb,epab,scftype,vifaa,vifbb,vifab,vifba,ncore,nbfcn,vifii,vifjj,
     $vifij,vifji)
************************************************************************
c Calculate intermediate B1 + intermediate V, F12/r12 contribution
************************************************************************
      implicit none
      integer nalc,nbec,dfnbasis,i,j,ij,nbfc,nbfcn,ncore
      real*8 cdelf122aa(dfnbasis,nalc,nalc),epaa(*),epbb(*),ddot
      real*8 cdelf122bb(dfnbasis,nbec,nbec),vifaa(nbfcn,dfnbasis,nalc)
      real*8 cdelf122ab(dfnbasis,nalc,nbec),vifbb(nbfcn,dfnbasis,nbec)
      real*8 cf12r12aa(dfnbasis,nalc,nalc),cf12r12bb(dfnbasis,nbec,nbec)
      real*8 cf12r12ab(dfnbasis,nalc,nbec),vifii(dfnbasis)
      real*8 vifab(nbfcn,dfnbasis,nalc),vifjj(dfnbasis),vifij(dfnbasis)
      real*8 vifba(nbfcn,dfnbasis,nbec),cdelf122ba(dfnbasis,nbec,nalc)
      real*8 cf12r12ba(dfnbasis,nbec,nalc),epab(nalc,nbec)
      real*8 vifji(dfnbasis)
      character(len=5) scftype
C
       if(scftype.eq.'rhf  ') then
        ij=0
         do i=1,nalc
          vifii=vifaa(ncore+i,1:dfnbasis,i)
           do j=1,i
            vifjj=vifaa(ncore+j,1:dfnbasis,j)
            vifij=vifaa(ncore+i,1:dfnbasis,j)
            ij=ij+1
            epaa(ij)=
C Intermediate B, term 1 (= tau)
     $      7d0*(ddot(dfnbasis,vifii            ,1,cdelf122aa(1,j,j),1)+
     $           ddot(dfnbasis,cdelf122aa(1,i,i),1,vifjj            ,1))
     $          +ddot(dfnbasis,vifij            ,1,cdelf122aa(1,i,j),1)
     $          +ddot(dfnbasis,cdelf122aa(1,i,j),1,vifij            ,1)
C Intermediate V, F12/r12 contribution
     $+8d0*(5d0*(ddot(dfnbasis,vifii           ,1,cf12r12aa(1,j,j),1)+
     $           ddot(dfnbasis,cf12r12aa(1,i,i),1,vifjj           ,1))
     $          -ddot(dfnbasis,vifij           ,1,cf12r12aa(1,i,j),1) 
     $          -ddot(dfnbasis,cf12r12aa(1,i,j),1,vifij           ,1))
           enddo
         enddo
       else
C alpha-alpha
        ij=0
         do i=1,nalc
          vifii=vifaa(ncore+i,1:dfnbasis,i)
           do j=1,i-1
            vifjj=vifaa(ncore+j,1:dfnbasis,j)
            vifij=vifaa(ncore+i,1:dfnbasis,j)
            ij=ij+1
            epaa(ij)=
C Intermediate B, term 1 (= tau)
     $           ddot(dfnbasis,vifii            ,1,cdelf122aa(1,j,j),1)
     $          +ddot(dfnbasis,cdelf122aa(1,i,i),1,vifjj            ,1)
     $          -ddot(dfnbasis,vifij            ,1,cdelf122aa(1,i,j),1)
     $          -ddot(dfnbasis,cdelf122aa(1,i,j),1,vifij            ,1)
C Intermediate V, F12/r12 contribution
     $     +8d0*(ddot(dfnbasis,vifii           ,1,cf12r12aa(1,j,j),1) 
     $          +ddot(dfnbasis,cf12r12aa(1,i,i),1,vifjj           ,1) 
     $          -ddot(dfnbasis,vifij           ,1,cf12r12aa(1,i,j),1)
     $          -ddot(dfnbasis,cf12r12aa(1,i,j),1,vifij           ,1))
           enddo
         enddo
C beta-beta
        ij=0
         do i=1,nbec
          vifii=vifbb(ncore+i,1:dfnbasis,i)
           do j=1,i-1
            vifjj=vifbb(ncore+j,1:dfnbasis,j)
            vifij=vifbb(ncore+i,1:dfnbasis,j)
            ij=ij+1
            epbb(ij)=
C Intermediate B, term 1 (= tau)
     $           ddot(dfnbasis,vifii            ,1,cdelf122bb(1,j,j),1)
     $          +ddot(dfnbasis,cdelf122bb(1,i,i),1,vifjj            ,1)
     $          -ddot(dfnbasis,vifij            ,1,cdelf122bb(1,i,j),1)
     $          -ddot(dfnbasis,cdelf122bb(1,i,j),1,vifij            ,1)
C Intermediate V, F12/r12 contribution
     $     +8d0*(ddot(dfnbasis,vifii           ,1,cf12r12bb(1,j,j),1) 
     $          +ddot(dfnbasis,cf12r12bb(1,i,i),1,vifjj           ,1) 
     $          -ddot(dfnbasis,vifij           ,1,cf12r12bb(1,i,j),1)
     $          -ddot(dfnbasis,cf12r12bb(1,i,j),1,vifij           ,1))
           enddo
         enddo
C alpha-beta
         do i=1,nalc
          vifii=vifaa(ncore+i,1:dfnbasis,i)
           do j=1,nbec
            vifjj=vifbb(ncore+j,1:dfnbasis,j)
            vifij=vifba(ncore+i,1:dfnbasis,j)
            vifji=vifab(ncore+j,1:dfnbasis,i)
            epab(i,j)=
C Intermediate B, term 1 (= tau)
     $ 5d0/32d0*(ddot(dfnbasis,vifii            ,1,cdelf122bb(1,j,j),1)
     $          +ddot(dfnbasis,cdelf122aa(1,i,i),1,vifjj            ,1))
     $+3d0/32d0*(ddot(dfnbasis,vifij            ,1,cdelf122ab(1,i,j),1)
     $          +ddot(dfnbasis,cdelf122ba(1,j,i),1,vifji            ,1))
C Intermediate V, F12/r12 contribution
     $ +3d0/4d0*(ddot(dfnbasis,vifii           ,1,cf12r12bb(1,j,j),1)
     $          +ddot(dfnbasis,cf12r12aa(1,i,i),1,vifjj           ,1))
     $ +1d0/4d0*(ddot(dfnbasis,vifij           ,1,cf12r12ab(1,i,j),1)
     $          +ddot(dfnbasis,cf12r12ba(1,j,i),1,vifji           ,1))
           enddo
         enddo
       endif
C
      return
      end
C
************************************************************************
      subroutine epair_stf122_bl(nalc,nbfcn,ncore,dfnbasis,fa,hpja,
     $vifaa,cf122aa,epaa,f12ij,scr,ifltln)
************************************************************************
* Calculate F12^2 contribution to MP2-F12 pair energies for RHF
* BLOCK version
************************************************************************
      implicit none
      integer nbfcn,nalc,i,j,ncore,ij,dfnbasis,ifltln
      real*8 fa(nbfcn,nbfcn),hpja(nbfcn,nbfcn),e12,epaa(*),ddot
      real*8 vifaa(nbfcn,dfnbasis,2),cf122aa(nbfcn,dfnbasis,nalc)
      real*8 f12ij(nbfcn,nbfcn),scr(nbfcn,nbfcn)
C Loop over pairs
      open(95,file='MPVIFAA2',access='direct',
     $recl=dfnbasis*nbfcn*ifltln,status='old')

      ij=0
       do i=1,nalc
        read(95,rec=i)vifaa(1:nbfcn,1:dfnbasis,1)
         do j=1,i
C Intermediate B, term 2 (= pi)
          read(95,rec=j)vifaa(1:nbfcn,1:dfnbasis,2)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,  vifaa(1,1,1),
     $nbfcn,cf122aa(1,1,j),nbfcn,0d0,scr,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,cf122aa(1,1,i),
     $nbfcn,  vifaa(1,1,2),nbfcn,1d0,scr,nbfcn)
          f12ij=transpose(scr)
          call daxpy(nbfcn*nbfcn,7d0,scr,1,f12ij,1)
          e12=    ddot(nbfcn,f12ij(1,ncore+j),1    ,hpja(1,ncore+i),1)
          e12=e12+ddot(nbfcn,f12ij(ncore+i,1),nbfcn,hpja(1,ncore+j),1)
C Intermediate X, F12^2 contribution
          e12=e12-ddot(nalc,f12ij(ncore+1,ncore+j),1    ,
     $fa(ncore+1,ncore+i),1)
          e12=e12-ddot(nalc,f12ij(ncore+i,ncore+1),nbfcn,
     $fa(ncore+1,ncore+j),1)
          ij=ij+1
          epaa(ij)=epaa(ij)+e12
         enddo
       enddo
      close(95)
C
      return
      end
C
************************************************************************
      subroutine epair_stf122(nalc,nbfcn,ncore,dfnbasis,fa,hpja,vifaa,
     $cf122aa,epaa,f12ij,scr)
************************************************************************
* Calculate F12^2 contribution to MP2-F12 pair energies for RHF
************************************************************************
      implicit none
      integer nbfcn,nalc,i,j,ncore,ij,dfnbasis
      real*8 fa(nbfcn,nbfcn),hpja(nbfcn,nbfcn),e12,epaa(*),ddot
      real*8 vifaa(nbfcn,dfnbasis,nalc),cf122aa(nbfcn,dfnbasis,nalc)
      real*8 f12ij(nbfcn,nbfcn),scr(nbfcn,nbfcn)
C Loop over pairs
      ij=0
       do i=1,nalc
         do j=1,i
C Intermediate B, term 2 (= pi)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,  vifaa(1,1,i),
     $nbfcn,cf122aa(1,1,j),nbfcn,0d0,scr,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,cf122aa(1,1,i),
     $nbfcn,  vifaa(1,1,j),nbfcn,1d0,scr,nbfcn)
          f12ij=transpose(scr)
          call daxpy(nbfcn*nbfcn,7d0,scr,1,f12ij,1)
          e12=    ddot(nbfcn,f12ij(1,ncore+j),1    ,hpja(1,ncore+i),1)
          e12=e12+ddot(nbfcn,f12ij(ncore+i,1),nbfcn,hpja(1,ncore+j),1)
C Intermediate X, F12^2 contribution
          e12=e12-ddot(nalc,f12ij(ncore+1,ncore+j),1    ,
     $fa(ncore+1,ncore+i),1)
          e12=e12-ddot(nalc,f12ij(ncore+i,ncore+1),nbfcn,
     $fa(ncore+1,ncore+j),1)
          ij=ij+1
          epaa(ij)=epaa(ij)+e12
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine epair_ssf122(nalc,nbfcn,ncore,dfnbasis,fa,hpja,vifaa,
     $cf122aa,epaa,f12ij,scr)
************************************************************************
* Calculate same-spin F12^2 contribution to MP2-F12 pair energies
************************************************************************
      implicit none
      integer nbfcn,nalc,i,j,ncore,ij,dfnbasis
      real*8 fa(nbfcn,nbfcn),hpja(nbfcn,nbfcn),e12,epaa(*),ddot
      real*8 vifaa(nbfcn,dfnbasis,nalc),cf122aa(nbfcn,dfnbasis,nalc)
      real*8 f12ij(nbfcn,nbfcn),scr(nbfcn,nbfcn)
C Loop over pairs
      ij=0
       do i=1,nalc
         do j=1,i-1
C Intermediate B, term 2 (= pi)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,  vifaa(1,1,i),
     $nbfcn,cf122aa(1,1,j),nbfcn,0d0,scr,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,1d0,cf122aa(1,1,i),
     $nbfcn,  vifaa(1,1,j),nbfcn,1d0,scr,nbfcn)
          f12ij=-transpose(scr)
          call daxpy(nbfcn*nbfcn,1d0,scr,1,f12ij,1)
          e12=    ddot(nbfcn,f12ij(1,ncore+j),1    ,hpja(1,ncore+i),1)
          e12=e12+ddot(nbfcn,f12ij(ncore+i,1),nbfcn,hpja(1,ncore+j),1)
C Intermediate X, F12^2 contribution
          e12=e12-ddot(nalc,f12ij(ncore+1,ncore+j),1    ,
     $fa(ncore+1,ncore+i),1)
          e12=e12-ddot(nalc,f12ij(ncore+i,ncore+1),nbfcn,
     $fa(ncore+1,ncore+j),1)
          ij=ij+1
          epaa(ij)=epaa(ij)+e12
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine epair_osf122(nalc,nbec,nbfcn,ncore,dfnbasis,fa,fb,hpja,
     $hpjb,vifaa,vifbb,vifab,vifba,cf122aa,cf122bb,cf122ab,cf122ba,epab,
     $f12ij)
************************************************************************
* Calculate opposite-spin F12^2 contribution to MP2-F12 pair energies
************************************************************************
      implicit none
      integer nbfcn,nalc,nbec,ncore,i,j,dfnbasis
      real*8 fa(nbfcn,nbfcn),fb(nbfcn,nbfcn),epab(nalc,nbec)
      real*8 hpja(nbfcn,nbfcn),hpjb(nbfcn,nbfcn),e12,ddot
      real*8 vifaa(nbfcn,dfnbasis,nalc),cf122aa(nbfcn,dfnbasis,nalc)
      real*8 vifbb(nbfcn,dfnbasis,nbec),cf122bb(nbfcn,dfnbasis,nbec)
      real*8 vifab(nbfcn,dfnbasis,nalc),cf122ab(nbfcn,dfnbasis,nalc)
      real*8 vifba(nbfcn,dfnbasis,nbec),cf122ba(nbfcn,dfnbasis,nbec)
      real*8 f12ij(nbfcn,nbfcn)
C
       do i=1,nalc
         do j=1,nbec
C Intermediate B, term 2 (= pi)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,5d0/32d0,
     $  vifaa(1,1,i),nbfcn,cf122bb(1,1,j),nbfcn,0d0,f12ij,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,5d0/32d0,
     $cf122aa(1,1,i),nbfcn,  vifbb(1,1,j),nbfcn,1d0,f12ij,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,3d0/32d0,
     $  vifba(1,1,j),nbfcn,cf122ab(1,1,i),nbfcn,1d0,f12ij,nbfcn)
          call dgemm('n','t',nbfcn,nbfcn,dfnbasis,3d0/32d0,
     $cf122ba(1,1,j),nbfcn,  vifab(1,1,i),nbfcn,1d0,f12ij,nbfcn)
          e12=    ddot(nbfcn,f12ij(1,ncore+j),1    ,hpja(1,ncore+i),1)
          e12=e12+ddot(nbfcn,f12ij(ncore+i,1),nbfcn,hpjb(1,ncore+j),1)
C Intermediate X, F12^2 contribution
          e12=e12-ddot(nalc,f12ij(ncore+1,ncore+j),1    ,
     $fa(ncore+1,ncore+i),1)
          e12=e12-ddot(nbec,f12ij(ncore+i,ncore+1),nbfcn,
     $fb(ncore+1,ncore+j),1)
          epab(i,j)=epab(i,j)+e12
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine f12r12cs_bl(nbfcn,nb,dfnbasis,nal,ncore,vifaa2,
     $cf12r12aa,vv,f12ij,iout,nel)
************************************************************************
* Assemble F12/r12 integrals for intermediate V, closed-shell
* BLOCK version
************************************************************************
      implicit none
      integer nbfcn,nb,dfnbasis,nal,ncore,i,j,ij,iout
      integer ib,nbb,kstart,nel,jj,k,ii,nbbl,nel2,iel
      real*8 vifaa2(nbfcn,dfnbasis,nel+1),cf12r12aa(nb,dfnbasis,nal)
      real*8 vv(nb,nb,(nal+1)*nal/2),f12ij(nb,nb)
C
      ij=0

C calc. # of blocks
      if(mod(nal,nel).ne.0)then
C if nal is not a multiple of nel
       nbbl = mod(nal,nel)
       nbb = 1+nal/nel
      else
C if nal is a multiple of nel
       nbb = nal/nel
       nbbl = nel
      endif

      do ib=1,nbb
C Read vifaa and cf12aa blocks
       kstart = (ib-1)*nel + 1

       nel2 = nel
       if(ib.eq.nbb)then
        nel2 = nbbl
       endif

       do k=kstart,kstart+nel2-1
        iel = k-(ib-1)*nel
        read(111,rec=k)vifaa2(1:nbfcn,1:dfnbasis,iel)
       enddo
c      write(iout,*)'ib: ',ib
       do jj=1,nel2 ! relative j index inside a block
        j = jj + (ib - 1)*nel ! absolute j index
           do i=1,j
C Check whether i is in memory
         ii = i - (ib - 1)*nel
         if(i.lt.kstart)then ! i is not in memory if it is lower than the start of the block 
          ii = nel+1
c         write(iout,*)'i read: ',i
          read(111,rec=i)vifaa2(1:nbfcn,1:dfnbasis,ii)
c         write(iout,*)'i read done: ',i
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ij=ij+1
c         write(iout,*)'dgemm 1'
         call dgemm('n','t',nb,nb,dfnbasis,1d0,vifaa2(ncore+1,1,ii),
     $nbfcn,cf12r12aa(1,1,j)  ,nb   ,0d0,f12ij,nb)
c         write(iout,*)'dgemm 2'
          call dgemm('n','t',nb,nb,dfnbasis,1d0,  cf12r12aa(1,1,i),
     $nb   ,vifaa2(ncore+1,1,jj),nbfcn,1d0,f12ij,nb)
c         write(iout,*)'dscal'
          call dscal(nb*nb,1d0/8d0,f12ij,1)
c         write(iout,*)'transpose: ', ij
          vv(:,:,ij)=transpose(f12ij)
c         write(iout,*)'daxpy'
          call daxpy(nb*nb,3d0,f12ij,1,vv(1,1,ij),1)
c         write(iout,*)'done'

c         call dgemm('n','t',nb,nb,dfnbasis,1d0,vifaa(ncore+1,1,i),
c    $nbfcn,cf12r12aa(1,1,j)  ,nb   ,0d0,f12ij,nb)
c         call dgemm('n','t',nb,nb,dfnbasis,1d0,  cf12r12aa(1,1,i),
c    $nb   ,vifaa(ncore+1,1,j),nbfcn,1d0,f12ij,nb)
c         call dscal(nb*nb,1d0/8d0,f12ij,1)
c         vv(:,:,ij)=transpose(f12ij)
c         call daxpy(nb*nb,3d0,f12ij,1,vv(1,1,ij),1)

         enddo
       enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine f12r12cs(nbfcn,nb,dfnbasis,nal,ncore,vifaa,cf12r12aa,
     $vv,f12ij)
************************************************************************
* Assemble F12/r12 integrals for intermediate V, closed-shell
************************************************************************
      implicit none
      integer nbfcn,nb,dfnbasis,nal,ncore,i,j,ij
      real*8 vifaa(nbfcn,dfnbasis,nal),cf12r12aa(nb,dfnbasis,nal)
      real*8 vv(nb,nb,(nal+1)*nal/2),f12ij(nb,nb)
C
      ij=0
       do j=1,nal
         do i=1,j
          ij=ij+1
          call dgemm('n','t',nb,nb,dfnbasis,1d0,vifaa(ncore+1,1,i),
     $nbfcn,cf12r12aa(1,1,j)  ,nb   ,0d0,f12ij,nb)
          call dgemm('n','t',nb,nb,dfnbasis,1d0,  cf12r12aa(1,1,i),
     $nb   ,vifaa(ncore+1,1,j),nbfcn,1d0,f12ij,nb)
          call dscal(nb*nb,1d0/8d0,f12ij,1)
          vv(:,:,ij)=transpose(f12ij)
          call daxpy(nb*nb,3d0,f12ij,1,vv(1,1,ij),1)
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine f12r12os(nbfcn,nb,dfnbasis,nal,nbe,ncore,vifaa,vifbb,
     $vifab,vifba,cf12r12aa,cf12r12bb,cf12r12ab,cf12r12ba,vaa,vbb,vab,
     $f12ij)
************************************************************************
* Assemble F12/r12 integrals for intermediate V, open-shell
************************************************************************
      implicit none
      integer nbfcn,nb,dfnbasis,nal,nbe,ncore,i,j,ij,p,q,pq
      real*8 vifaa(nbfcn,dfnbasis,nal),cf12r12aa(nb,dfnbasis,nal)
      real*8 vifbb(nbfcn,dfnbasis,nbe),cf12r12bb(nb,dfnbasis,nbe)
      real*8 vifab(nbfcn,dfnbasis,nal),cf12r12ab(nb,dfnbasis,nal)
      real*8 vifba(nbfcn,dfnbasis,nbe),cf12r12ba(nb,dfnbasis,nbe)
      real*8 vaa((nb-1)*nb/2,(nal-1)*nal/2),vab(nb,nb,nal,nbe)
      real*8 vbb((nb-1)*nb/2,(nbe-1)*nbe/2),f12ij(nb,nb)
C
      ij=0
       do j=1,nal
         do i=1,j-1
          ij=ij+1
          call dgemm('n','t',nb,nb,dfnbasis,1d0,vifaa(ncore+1,1,i),
     $nbfcn,cf12r12aa(1,1,j)  ,nb   ,0d0,f12ij,nb)
          call dgemm('n','t',nb,nb,dfnbasis,1d0,  cf12r12aa(1,1,i),
     $nb   ,vifaa(ncore+1,1,j),nbfcn,1d0,f12ij,nb)
          pq=0
           do q=1,nb
             do p=1,q-1
              pq=pq+1
              vaa(pq,ij)=(f12ij(p,q)-f12ij(q,p))/4d0
             enddo
           enddo
         enddo
       enddo
      ij=0
       do j=1,nbe
         do i=1,j-1
          ij=ij+1
          call dgemm('n','t',nb,nb,dfnbasis,1d0,vifbb(ncore+1,1,i),
     $nbfcn,cf12r12bb(1,1,j)  ,nb   ,0d0,f12ij,nb)
          call dgemm('n','t',nb,nb,dfnbasis,1d0,  cf12r12bb(1,1,i),
     $nb   ,vifbb(ncore+1,1,j),nbfcn,1d0,f12ij,nb)
          pq=0
           do q=1,nb
             do p=1,q-1
              pq=pq+1
              vbb(pq,ij)=(f12ij(p,q)-f12ij(q,p))/4d0
             enddo
           enddo
         enddo
       enddo
       do i=1,nal
         do j=1,nbe
         call dgemm('n','t',nb,nb,dfnbasis,3d0/8d0,vifaa(ncore+1,1,i),
     $nbfcn,cf12r12bb(1,1,j)  ,nb   ,0d0,vab(1,1,i,j),nb)
         call dgemm('n','t',nb,nb,dfnbasis,3d0/8d0,  cf12r12aa(1,1,i),
     $nb   ,vifbb(ncore+1,1,j),nbfcn,1d0,vab(1,1,i,j),nb)
         call dgemm('n','t',nb,nb,dfnbasis,1d0/8d0,vifba(ncore+1,1,j),
     $nbfcn,cf12r12ab(1,1,i),  nb   ,1d0,vab(1,1,i,j),nb)
         call dgemm('n','t',nb,nb,dfnbasis,1d0/8d0,  cf12r12ba(1,1,j),
     $nb   ,vifab(ncore+1,1,i),nbfcn,1d0,vab(1,1,i,j),nb)
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine vterm1_bl(nbfcn,nb,nnb,dfnbasis,nal,vifaa2,cp,cm,vs,
     $va,vv,f12pq,mns,mna,rsdims,rsdima,slo,sup,nel,iout)
************************************************************************
* Calculate PPL-like contribution to intermediate V, term 1
* BLOCK version
************************************************************************
      implicit none
      integer nbfcn,nb,nnb,dfnbasis,nal,i,j,ij,p,q,pq,r,s,rs,rsdims,slo
      integer sup,rsdima,rss,rsa,ijs,ija,nel,iout,ii,jj,kk
      real*8 vifaa2(nbfcn,dfnbasis,nel+1),tmp
      real*8 cp((nnb+1)*nnb/2*(nal+1)*nal/2),vv(nb,nb,(nal+1)*nal/2)
      real*8 cm((nnb-1)*nnb/2*(nal-1)*nal/2),f12pq(nnb,nnb)
      real*8 vs((nnb+1)*nnb/2,rsdims)
      real*8 va((nnb-1)*nnb/2,rsdima)
      real*8 mns(rsdims,(nal+1)*nal/2),mna(rsdima,(nal+1)*nal/2)

C Build (anti)symmetrized Coulomb integrals
c     write(iout,*) 'inside vterm1_bl '
      rss=0
      rsa=0

c     write(iout,*) 'slo, sup ',slo,sup
c     write(iout,*) 'f12pq dec size: ', nnb*nnb
c     write(iout,*) 'first address: ',loc(f12pq(1,1))
c     write(iout,*) 'last address: ',loc(f12pq(nnb,nnb))
       do s=slo,sup
        read(111,rec=s)vifaa2(1:nbfcn,1:dfnbasis,2)
         do r=1,s-1
          read(111,rec=r)vifaa2(1:nbfcn,1:dfnbasis,1)
c         write(iout,*) 's,r ',s,r
C ha r kisebb mint slo akkor be kell huzni, de ha a blokk slo-tol sup-ig bent van, akkor a blokkbol kell kulon olvasni
          rss=rss+1
          rsa=rsa+1

c          do ii=1,nnb
c          do jj=1,nnb
c         f12pq(ii,jj) = 0.d0
c          do kk=1,dfnbasis
c         f12pq(ii,jj) = f12pq(ii,jj) + 
c    $                   vifaa2(ii,kk,1)*vifaa2(kk,jj,2)
c          enddo
c          enddo
c          enddo
          call dgemm('n','t',nnb,nnb,dfnbasis,1d0,vifaa2(1,1,1),
     $nbfcn,vifaa2(1,1,2),nbfcn,0d0,f12pq,nnb)

           pq=0
           do q=1,nnb
             do p=1,q
              pq=pq+1
              vs(pq,rss)=f12pq(p,q)+f12pq(q,p)

             enddo
           enddo

           pq=0
           do q=1,nnb
             do p=1,q-1
              pq=pq+1
              va(pq,rsa)=f12pq(p,q)-f12pq(q,p)
             enddo
           enddo
         enddo ! r
        rss=rss+1
c         write(iout,*) 'calling dsyrk '
        call dsyrk('U','N',nnb,dfnbasis,2d0,vifaa2(1,1,2),nbfcn,0d0,
     $f12pq,nnb)
         pq=0
         do q=1,nnb
           do p=1,q
            pq=pq+1
            vs(pq,rss)=f12pq(p,q)
           enddo
         enddo
       enddo
C Term 2
c         write(iout,*) 'calling vppl '
      call vppl(nb,nal,cp,cm,vs,va,vv,mns,mna,rsdims,rsdima,slo,
     $sup,(nnb+1)*nnb/2,(nnb-1)*nnb/2)
C
      return
      end
C
C
************************************************************************
      subroutine vterm1(nbfcn,nb,nnb,dfnbasis,nal,vifaa,cp,cm,vs,
     $va,vv,f12pq,mns,mna,rsdims,rsdima,slo,sup)
************************************************************************
* Calculate PPL-like contribution to intermediate V, term 1
************************************************************************
      implicit none
      integer nbfcn,nb,nnb,dfnbasis,nal,i,j,ij,p,q,pq,r,s,rs,rsdims,slo
      integer sup,rsdima,rss,rsa,ijs,ija
      real*8 vifaa(nbfcn,dfnbasis,nb),tmp
      real*8 cp((nnb+1)*nnb/2*(nal+1)*nal/2),vv(nb,nb,(nal+1)*nal/2)
      real*8 cm((nnb-1)*nnb/2*(nal-1)*nal/2),f12pq(nnb,nnb)
      real*8 vs((nnb+1)*nnb/2,rsdims)
      real*8 va((nnb-1)*nnb/2,rsdima)
      real*8 mns(rsdims,(nal+1)*nal/2),mna(rsdima,(nal+1)*nal/2)
C Build (anti)symmetrized Coulomb integrals
      rss=0
      rsa=0
       do s=slo,sup
         do r=1,s-1
          rss=rss+1
          rsa=rsa+1
          call dgemm('n','t',nnb,nnb,dfnbasis,1d0,vifaa(1,1,r),
     $nbfcn,vifaa(1,1,s),nbfcn,0d0,f12pq,nnb)
           pq=0
           do q=1,nnb
             do p=1,q
              pq=pq+1
              vs(pq,rss)=f12pq(p,q)+f12pq(q,p)
             enddo
           enddo
           pq=0
           do q=1,nnb
             do p=1,q-1
              pq=pq+1
              va(pq,rsa)=f12pq(p,q)-f12pq(q,p)
             enddo
           enddo
         enddo
        rss=rss+1
        call dsyrk('U','N',nnb,dfnbasis,2d0,vifaa(1,1,s),nbfcn,0d0,
     $f12pq,nnb)
         pq=0
         do q=1,nnb
           do p=1,q
            pq=pq+1
            vs(pq,rss)=f12pq(p,q)
           enddo
         enddo
       enddo
C Term 2
      call vppl(nb,nal,cp,cm,vs,va,vv,mns,mna,rsdims,rsdima,slo,
     $sup,(nnb+1)*nnb/2,(nnb-1)*nnb/2)
C
      return
      end
C
************************************************************************
      subroutine vterm2_bl(nbfcn,nb,nnb,dfnbasis,nal,vifaa2,cp,cm,vs,
     $va,vv,f12ao,mns,mna,rsdims,rsdima,slo,sup,ncore,nbfan,aodim,nel)
************************************************************************
* Calculate PPL-like contribution to intermediate V, term 2
* BLOCK version
************************************************************************
      implicit none
      integer nbfcn,nb,nnb,dfnbasis,nal,i,j,ij,p,q,pq,r,s,rs,rsdims,slo
      integer sup,rsdima,rss,rsa,ijs,ija,a,o,ncore,nbfan,aodim,nel

      real*8 vifaa2(nbfcn,dfnbasis,nel+1),tmp

      real*8 cp(nbfan,ncore+nal,(nal+1)*nal/2),vv(nb,nb,(nal+1)*nal/2)
      real*8 cm(nbfan,ncore+nal,(nal-1)*nal/2),f12ao(nbfan,ncore+nal)
      real*8 vs(nbfan,ncore+nal,rsdims)
      real*8 va(nbfan,ncore+nal,rsdima)
      real*8 mns(rsdims,(nal+1)*nal/2),mna(rsdima,(nal+1)*nal/2)

C Build (anti)symmetrized Coulomb integrals
      rss=0
      rsa=0
       do s=slo,sup
        read(111,rec=s)vifaa2(1:nbfcn,1:dfnbasis,2)
         do r=1,s-1
          read(111,rec=r)vifaa2(1:nbfcn,1:dfnbasis,1)
          rss=rss+1
          rsa=rsa+1
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $vifaa2(nnb+1,1,1),nbfcn,vifaa2(1,1,2),nbfcn,0d0,vs(1,1,rss),nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $vifaa2(nnb+1,1,2),nbfcn,vifaa2(1,1,1),nbfcn,0d0,f12ao      ,nbfan)
          call dcopy(aodim, vs(1,1,rss),1,va(1,1,rsa),1)
          call daxpy(aodim,   1d0,f12ao,1,vs(1,1,rss),1)
          call daxpy(aodim,  -1d0,f12ao,1,va(1,1,rsa),1)
         enddo
        rss=rss+1
        call dgemm('n','t',nbfan,ncore+nal,dfnbasis,2d0,
     $vifaa2(nnb+1,1,2),nbfcn,vifaa2(1,1,2),nbfcn,0d0,vs(1,1,rss),nbfan)
       enddo
C Term 2
      call vppl(nb,nal,cp,cm,vs,va,vv,mns,mna,rsdims,rsdima,slo,
     $sup,aodim,aodim)
C
      return
      end
C
************************************************************************
      subroutine vterm2(nbfcn,nb,nnb,dfnbasis,nal,vifaa,cp,cm,vs,
     $va,vv,f12ao,mns,mna,rsdims,rsdima,slo,sup,ncore,nbfan,aodim)
************************************************************************
* Calculate PPL-like contribution to intermediate V, term 2
************************************************************************
      implicit none
      integer nbfcn,nb,nnb,dfnbasis,nal,i,j,ij,p,q,pq,r,s,rs,rsdims,slo
      integer sup,rsdima,rss,rsa,ijs,ija,a,o,ncore,nbfan,aodim
      real*8 vifaa(nbfcn,dfnbasis,nb),tmp
      real*8 cp(nbfan,ncore+nal,(nal+1)*nal/2),vv(nb,nb,(nal+1)*nal/2)
      real*8 cm(nbfan,ncore+nal,(nal-1)*nal/2),f12ao(nbfan,ncore+nal)
      real*8 vs(nbfan,ncore+nal,rsdims)
      real*8 va(nbfan,ncore+nal,rsdima)
      real*8 mns(rsdims,(nal+1)*nal/2),mna(rsdima,(nal+1)*nal/2)
C Build (anti)symmetrized Coulomb integrals
      rss=0
      rsa=0
       do s=slo,sup
         do r=1,s-1
          rss=rss+1
          rsa=rsa+1
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $vifaa(nnb+1,1,r),nbfcn,vifaa(1,1,s),nbfcn,0d0,vs(1,1,rss),nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,1d0,
     $vifaa(nnb+1,1,s),nbfcn,vifaa(1,1,r),nbfcn,0d0,f12ao      ,nbfan)
          call dcopy(aodim, vs(1,1,rss),1,va(1,1,rsa),1)
          call daxpy(aodim,   1d0,f12ao,1,vs(1,1,rss),1)
          call daxpy(aodim,  -1d0,f12ao,1,va(1,1,rsa),1)
         enddo
        rss=rss+1
        call dgemm('n','t',nbfan,ncore+nal,dfnbasis,2d0,
     $vifaa(nnb+1,1,s),nbfcn,vifaa(1,1,s),nbfcn,0d0,vs(1,1,rss),nbfan)
       enddo
C Term 2
      call vppl(nb,nal,cp,cm,vs,va,vv,mns,mna,rsdims,rsdima,slo,
     $sup,aodim,aodim)
C
      return
      end
C
************************************************************************
      subroutine vppl(nb,nal,cp,cm,vs,va,vv,mns,mna,rsdims,rsdima,slo,
     $sup,aodims,aodima)
************************************************************************
* Calculate PPL-like contribution to intermediate V
************************************************************************
      implicit none
      integer nb,nal,i,j,ij,r,s,rs,rsdims,slo,sup,rsdima,ijs,ija,aodims
      integer aodima
      real*8 tmp,vv(nb,nb,(nal+1)*nal/2)
      real*8 cp(aodims,(nal+1)*nal/2),vs(aodims,rsdims)
      real*8 cm(aodima,(nal-1)*nal/2),va(aodima,rsdima)
      real*8 mns(rsdims,(nal+1)*nal/2),mna(rsdima,(nal+1)*nal/2)
C
      call dgemm('t','n',rsdims,(nal+1)*nal/2,aodims,1d0,vs,aodims,cp,
     $aodims,0d0,mns,rsdims)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,rs,r,s,tmp)
       do ij=1,(nal+1)*nal/2
        rs=0
         do s=slo,sup
           do r=1,s-1
            rs=rs+1
            tmp=mns(rs,ij)
            vv(r,s,ij)=vv(r,s,ij)-tmp
            vv(s,r,ij)=vv(s,r,ij)-tmp
           enddo
          rs=rs+1
          tmp=mns(rs,ij)
          vv(s,s,ij)=vv(s,s,ij)-tmp
         enddo
       enddo
C$OMP END PARALLEL DO
      call dgemm('t','n',rsdima,(nal-1)*nal/2,aodima,1d0,va,aodima,cm,
     $aodima,0d0,mna,rsdima)
C$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,ijs,ija,rs,r,s,tmp)
       do j=1,nal
        ijs=(j-1)*j/2
        ija=(j-1)*(j-2)/2
         do i=1,j-1
          ijs=ijs+1
          ija=ija+1
          rs=0
           do s=slo,sup
             do r=1,s-1
              rs=rs+1
              tmp=mna(rs,ija)
              vv(r,s,ijs)=vv(r,s,ijs)-tmp
              vv(s,r,ijs)=vv(s,r,ijs)+tmp
             enddo
           enddo
         enddo
       enddo
C$OMP END PARALLEL DO
C
      return
      end
C
************************************************************************
      subroutine cabijcalc_bl(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,
     $cf12aa2,f12caj,f12acj,f12abk,vinabk,f12ao,focka,cabij,vv,
     $vifaa2,nel,vinlog,vinabk2,nbb,nbbl)
************************************************************************
* Calculate Cabij intermediate, closed-shell 
* BLOCK version
************************************************************************
      implicit none
      integer nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,i,j,k,ij,jk,kl,a,b
      integer nel,jblock,jj,kk,kblock,vinlog,jind,nbb,nel2,nel3,nbbl
      real*8 focka(nbfcn,nbfcn)

c     real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)

      real*8 cf12aa2(nbfcn,dfnbasis,nel+1)

      real*8 vifaa2(nbfcn,dfnbasis,nel+1)

C ezeket is blokkositani
      real*8 f12caj(nbfan,nval,nal)
      real*8 f12acj(nbfan,nval,nal)
      real*8 f12abk(nbfan,nval,nal)
      real*8 vinabk(nval,nbfan,nal)
      real*8 vinabk2(nval,nbfan,nal)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 f12ao(nbfan,nval)
C
      real*8 vv(nb,nb,(nal+1)*nal/2)
      real*8 cabij(nval,nval,(nal+1)*nal/2)
C
      cabij(:,:,:)=vv(nal+1:nal+nval,nal+1:nal+nval,:)

      write(6,*)'nbfan: ', nbfan
      write(6,*)'nval: ', nval
      write(6,*)'nal: ', nal
      write(6,*)'nbfcn,dfnbasis,nal: ', nbfcn*dfnbasis*nal

      jk=0

      do kblock=1,nbb

       nel2 = nel
       if(kblock.eq.nbb)then
        nel2 = nbbl
       endif

       do kk=1,nel2
        k = (kblock - 1) * nel + kk
        read(111,rec=k)vifaa2(1:nbfcn,1:dfnbasis,kk)
        read(110,rec=k)cf12aa2(1:nbfcn,1:dfnbasis,kk)
       enddo

       do kk=1,nel2
        k = (kblock - 1) * nel + kk

         do jblock=1,nbb

       nel3 = nel
       if(jblock.eq.nbb)then
        nel3 = nbbl
       endif

          do jj=1,nel3
           j = (jblock - 1) * nel + jj

           if(kblock.eq.jblock)then
            jind = jj
           else
            jind = nel + 1
            read(111,rec=j)vifaa2(1:nbfcn,1:dfnbasis,jind)
            read(110,rec=j)cf12aa2(1:nbfcn,1:dfnbasis,jind)
           endif


          if(vinlog.eq.1)then ! if there is enough mem for vinabk2
           call dgemm('n','t',nval,nbfan,dfnbasis,1d0,
     $vifaa2(ncore+nal+1,1,jind),nbfcn,vifaa2(ncore+nb+1,1,kk),
     $nbfcn,0d0,vinabk2(1,1,j),nval)
          endif

          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $vifaa2(ncore+nb+1,1,jind),nbfcn,cf12aa2(ncore+nal+1,1,kk),nbfcn,
     $0d0,f12caj(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa2(ncore+nb+1,1,jind),nbfcn, 
     $vifaa2(ncore+nal+1,1,kk),nbfcn,1d0,f12caj(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0, 
     $vifaa2(ncore+nb+1,1,kk),nbfcn,cf12aa2(ncore+nal+1,1,jind),nbfcn,
     $0d0,f12acj(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa2(ncore+nb+1,1,kk),nbfcn, vifaa2(ncore+nal+1,1,jind),nbfcn,
     $1d0,f12acj(1,1,j),nbfan)
          enddo ! jj
         enddo ! jblock

        call dscal(nbfan*nval*nal,1d0/8d0,f12caj,1)
        call dscal(nbfan*nval*nal,1d0/8d0,f12acj,1)
        call dcopy(nbfan*nval*nal,f12acj,1,f12abk,1)
        call daxpy(nbfan*nval*nal,3d0,f12caj,1,f12abk,1)

         do a=1,nval
          read(111,rec=nal+a)vifaa2(1:nbfcn,1:dfnbasis,nel+1)
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $vifaa2(ncore+nb+1,1,nel+1),nbfcn,vifaa2(ncore+1,1,kk),nbfcn,0d0,
     $f12ao,nbfan) 
          call dcopy(nbfan*nal,f12ao,1,vinabk(a,1,1),nval)
         enddo

        ij=0
        kl=jk

         do j=1,nal
           if(j.le.k) then
            jk=jk+1
            call dgemm('n','n',nval,nval,nbfan,1d0,
     $focka(ncore+nal+1,ncore+nb+1),nbfcn,f12abk(1,1,j),nbfan,1d0,
     $cabij(1,1,jk),nval)
           endif
           do i=1,j
            ij=ij+1
            call dgemm('n','n',nval,nval,nbfan,-1d0,
     $vinabk(1,1,j),nval,f12abk(1,1,i),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,-1d0,
     $f12abk(1,1,j),nbfan,vinabk(1,1,i),nval,1d0,cabij(1,1,ij),nval)
           enddo
         enddo

        call dcopy(nbfan*nval*nal,f12caj,1,f12abk,1)
        call daxpy(nbfan*nval*nal,3d0,f12acj,1,f12abk,1)
        ij=0
        jk=kl
         do j=1,nal
           if(j.le.k) then
            jk=jk+1
            call dgemm('t','n',nval,nval,nbfan,1d0,
     $f12abk(1,1,j),nbfan,focka(ncore+nb+1,ncore+nal+1),nbfcn,1d0,
     $cabij(1,1,jk),nval)
           endif
           do i=1,j
            ij=ij+1
            call dgemm('n','n',nval,nval,nbfan,-1d0,
     $vinabk(1,1,i),nval,f12abk(1,1,j),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,-1d0,
     $f12abk(1,1,i),nbfan,vinabk(1,1,j),nval,1d0,cabij(1,1,ij),nval)
           enddo
         enddo
        call dcopy(nbfan*nval*nal,f12caj,1,f12abk,1)
        f12abk=-f12abk
        call daxpy(nbfan*nval*nal,5d0,f12acj,1,f12abk,1)
        ij=0


         do jblock=1,nbb

       nel3 = nel
       if(jblock.eq.nbb)then
        nel3 = nbbl
       endif

          do jj=1,nel3
           j = (jblock - 1) * nel + jj

          if(vinlog.eq.0)then ! if there is not enough mem for vinabk2

           if(kblock.eq.jblock)then
            jind = jj
           else
            jind = nel + 1
            read(111,rec=j)vifaa2(1:nbfcn,1:dfnbasis,jind)
           endif

          call dgemm('n','t',nval,nbfan,dfnbasis,1d0,
     $vifaa2(ncore+nal+1,1,jind),nbfcn,vifaa2(ncore+nb+1,1,kk),
     $nbfcn,0d0,vinabk(1,1,j),nval)
          endif

           do i=1,j
            ij=ij+1

          if(vinlog.eq.0)then ! if there is not enough mem for vinabk2
            call dgemm('n','n',nval,nval,nbfan,1d0,
     $vinabk(1,1,i),nval,f12abk(1,1,j),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,1d0,
     $f12abk(1,1,i),nbfan,vinabk(1,1,j),nval,1d0,cabij(1,1,ij),nval)
            endif

            if(vinlog.eq.1)then ! if there is enough mem for vinabk2
            call dgemm('n','n',nval,nval,nbfan,1d0,
     $vinabk2(1,1,i),nval,f12abk(1,1,j),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,1d0,
     $f12abk(1,1,i),nbfan,vinabk2(1,1,j),nval,1d0,cabij(1,1,ij),nval)
            endif

           enddo
         enddo ! jj
         enddo ! jblock

        enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine cabijcalc(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,vifaa,
     $cf12aa,f12caj,f12acj,f12abk,vinabk,f12ao,focka,cabij,vv)
************************************************************************
* Calculate Cabij intermediate, closed-shell
************************************************************************
      implicit none
      integer nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,i,j,k,ij,jk,kl,a,b
      real*8 focka(nbfcn,nbfcn)
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
C ezeket is blokkositani
      real*8 f12caj(nbfan,nval,nal),f12ao(nbfan,nval)
      real*8 f12acj(nbfan,nval,nal)
      real*8 f12abk(nbfan,nval,nal),vinabk(nval,nbfan,nal)
C
      real*8 vv(nb,nb,(nal+1)*nal/2)
      real*8 cabij(nval,nval,(nal+1)*nal/2)
C
      cabij(:,:,:)=vv(nal+1:nal+nval,nal+1:nal+nval,:)
      jk=0
       do k=1,nal
         do j=1,nal
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $ vifaa(ncore+nb+1,1,j),nbfcn,cf12aa(ncore+nal+1,1,k),nbfcn,0d0,
     $f12caj(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa(ncore+nb+1,1,j),nbfcn, vifaa(ncore+nal+1,1,k),nbfcn,1d0,
     $f12caj(1,1,j),nbfan)
         enddo
        call dscal(nbfan*nval*nal,1d0/8d0,f12caj,1)
         do j=1,nal
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0, 
     $ vifaa(ncore+nb+1,1,k),nbfcn,cf12aa(ncore+nal+1,1,j),nbfcn,0d0,
     $f12acj(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa(ncore+nb+1,1,k),nbfcn, vifaa(ncore+nal+1,1,j),nbfcn,1d0,
     $f12acj(1,1,j),nbfan)
         enddo
        call dscal(nbfan*nval*nal,1d0/8d0,f12acj,1)
        call dcopy(nbfan*nval*nal,f12acj,1,f12abk,1)
        call daxpy(nbfan*nval*nal,3d0,f12caj,1,f12abk,1)
         do a=1,nval
          call dgemm('n','t',nbfan,nal,dfnbasis,1d0,
     $vifaa(ncore+nb+1,1,nal+a),nbfcn,vifaa(ncore+1,1,k),nbfcn,0d0,
     $f12ao,nbfan)
          call dcopy(nbfan*nal,f12ao,1,vinabk(a,1,1),nval)
         enddo
        ij=0
        kl=jk
         do j=1,nal
           if(j.le.k) then
            jk=jk+1
            call dgemm('n','n',nval,nval,nbfan,1d0,
     $focka(ncore+nal+1,ncore+nb+1),nbfcn,f12abk(1,1,j),nbfan,1d0,
     $cabij(1,1,jk),nval)
           endif
           do i=1,j
            ij=ij+1
            call dgemm('n','n',nval,nval,nbfan,-1d0,
     $vinabk(1,1,j),nval,f12abk(1,1,i),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,-1d0,
     $f12abk(1,1,j),nbfan,vinabk(1,1,i),nval,1d0,cabij(1,1,ij),nval)
           enddo
         enddo
        call dcopy(nbfan*nval*nal,f12caj,1,f12abk,1)
        call daxpy(nbfan*nval*nal,3d0,f12acj,1,f12abk,1)
        ij=0
        jk=kl
         do j=1,nal
           if(j.le.k) then
            jk=jk+1
            call dgemm('t','n',nval,nval,nbfan,1d0,
     $f12abk(1,1,j),nbfan,focka(ncore+nb+1,ncore+nal+1),nbfcn,1d0,
     $cabij(1,1,jk),nval)
           endif
           do i=1,j
            ij=ij+1
            call dgemm('n','n',nval,nval,nbfan,-1d0,
     $vinabk(1,1,i),nval,f12abk(1,1,j),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,-1d0,
     $f12abk(1,1,i),nbfan,vinabk(1,1,j),nval,1d0,cabij(1,1,ij),nval)
           enddo
         enddo
        call dcopy(nbfan*nval*nal,f12caj,1,f12abk,1)
        f12abk=-f12abk
        call daxpy(nbfan*nval*nal,5d0,f12acj,1,f12abk,1)
        ij=0
         do j=1,nal
          call dgemm('n','t',nval,nbfan,dfnbasis,1d0,
     $vifaa(ncore+nal+1,1,j),nbfcn,vifaa(ncore+nb+1,1,k),nbfcn,0d0,
     $vinabk(1,1,j),nval)
           do i=1,j
            ij=ij+1
            call dgemm('n','n',nval,nval,nbfan,1d0,
     $vinabk(1,1,i),nval,f12abk(1,1,j),nbfan,1d0,cabij(1,1,ij),nval)
            call dgemm('t','t',nval,nval,nbfan,1d0,
     $f12abk(1,1,i),nbfan,vinabk(1,1,j),nval,1d0,cabij(1,1,ij),nval)
           enddo
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine vaacalc(nb,nnb,ncore,nal,nbfcn,nbfan,dfnbasis,dcore,
     $imem,imem1,maxcor,vifaa,cf12aa,f12pq,vaa)
************************************************************************
* Calculate V intermediate, same-spin block
************************************************************************
      implicit none
      integer nb,ncore,nal,nbfcn,nbfan,dfnbasis,imem,imem1,maxcor,s,r,rs
      integer i,j,ij,p,q,pq,aodim,ijdim,rsdim,icm,iva,dblalloc,slo,sup
      integer maxmem,nnb
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
      real*8 dcore(*),f12pq(nnb,nnb),vaa((nb-1)*nb/2,(nal-1)*nal/2)
C Calculate PPL-like contribution to intermediate V, term 1
C Build antisymmetrized F12 integrals
      aodim=(nnb-1)*nnb/2
      ijdim=(nal-1)*nal/2
      icm=dblalloc(aodim*ijdim)
      pq=icm-1
       do j=1,nal
         do i=1,j-1
          call dgemm('n','t',nnb,nnb,dfnbasis,1d0, vifaa(1,1,i),
     $nbfcn,cf12aa(1,1,j),nbfcn,0d0,f12pq,nnb)
          call dgemm('n','t',nnb,nnb,dfnbasis,1d0,cf12aa(1,1,i),
     $nbfcn, vifaa(1,1,j),nbfcn,1d0,f12pq,nnb)
           do q=1,nnb
             do p=1,q-1
              pq=pq+1
              dcore(pq)=(f12pq(p,q)-f12pq(q,p))/4d0
             enddo
           enddo
         enddo
       enddo
C Loop over blocks
      maxmem=maxcor-(imem-imem1)
      slo=1
      rs=1
       do sup=2,nb
        rsdim=(sup-slo+1)*(slo+sup-2)/2
         if(sup.eq.nb.or.(rsdim+sup)*aodim.gt.maxmem) then
          iva=dblalloc(aodim*rsdim)
          pq=iva-1
           do s=slo,sup
             do r=1,s-1
              call dgemm('n','t',nnb,nnb,dfnbasis,1d0,vifaa(1,1,r),
     $nbfcn,vifaa(1,1,s),nbfcn,0d0,f12pq,nnb)
               do q=1,nnb
                 do p=1,q-1
                  pq=pq+1
                  dcore(pq)=f12pq(p,q)-f12pq(q,p)
                 enddo
               enddo
             enddo
           enddo
          call dgemm('t','n',rsdim,ijdim,aodim,-1d0,dcore(iva),
     $aodim,dcore(icm),aodim,1d0,vaa(rs,1),(nb-1)*nb/2)
          call dbldealloc(iva)
          slo=sup+1
          rs=rs+rsdim
         endif
       enddo
      call dbldealloc(icm)
C Calculate PPL-like contribution to intermediate V, term 2
C Build antisymmetrized F12 integrals
      aodim=nbfan*(ncore+nal)
      icm=dblalloc(aodim*ijdim)
      ij=icm
       do j=1,nal
         do i=1,j-1
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis, 1d0,
     $ vifaa(nnb+1,1,i),nbfcn,cf12aa(1,1,j),nbfcn,0d0,dcore(ij),nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis, 1d0,
     $cf12aa(nnb+1,1,i),nbfcn, vifaa(1,1,j),nbfcn,1d0,dcore(ij),nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,-1d0,
     $ vifaa(nnb+1,1,j),nbfcn,cf12aa(1,1,i),nbfcn,1d0,dcore(ij),nbfan)
          call dgemm('n','t',nbfan,ncore+nal,dfnbasis,-1d0,
     $cf12aa(nnb+1,1,j),nbfcn, vifaa(1,1,i),nbfcn,1d0,dcore(ij),nbfan)
          call dscal(aodim,0.25d0,dcore(ij),1)
          ij=ij+aodim
         enddo
       enddo
C Loop over blocks
      maxmem=maxcor-(imem-imem1)
      slo=1
      rs=1
       do sup=2,nb
        rsdim=(sup-slo+1)*(slo+sup-2)/2
         if(sup.eq.nb.or.(rsdim+sup)*aodim.gt.maxmem) then
          iva=dblalloc(aodim*rsdim)
          pq=iva
           do s=slo,sup
             do r=1,s-1
              call dgemm('n','t',nbfan,ncore+nal,dfnbasis, 1d0,
     $vifaa(nnb+1,1,r),nbfcn,vifaa(1,1,s),nbfcn,0d0,dcore(pq),nbfan)
              call dgemm('n','t',nbfan,ncore+nal,dfnbasis,-1d0,
     $vifaa(nnb+1,1,s),nbfcn,vifaa(1,1,r),nbfcn,1d0,dcore(pq),nbfan)
              pq=pq+aodim
             enddo
           enddo
          call dgemm('t','n',rsdim,ijdim,aodim,-1d0,dcore(iva),
     $aodim,dcore(icm),aodim,1d0,vaa(rs,1),(nb-1)*nb/2)
          call dbldealloc(iva)
          slo=sup+1
          rs=rs+rsdim
         endif
       enddo
      call dbldealloc(icm)
C
      return
      end
C
************************************************************************
      subroutine vabcalc(nb,nnb,ncore,nal,nbe,nbfcn,nbfan,dfnbasis,
     $dcore,imem,imem1,maxcor,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,
     $cf12ab,cf12ba,vab)
************************************************************************
* Calculate V intermediate, opposite-spin block
************************************************************************
      implicit none
      integer nb,ncore,nal,nbfcn,nbfan,dfnbasis,imem,imem1,maxcor,s,r,rs
      integer i,j,ij,p,q,pq,aodim,ijdim,rsdim,icm,iva,dblalloc,slo,sup
      integer maxmem,nnb,nbe
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
      real*8 vifbb(nbfcn,dfnbasis,nb),cf12bb(nbfcn,dfnbasis,nbe)
      real*8 vifab(nbfcn,dfnbasis,nb),cf12ab(nbfcn,dfnbasis,nal)
      real*8 vifba(nbfcn,dfnbasis,nb),cf12ba(nbfcn,dfnbasis,nbe)
      real*8 dcore(*),vab(nb,nb,nal,nbe)
C Calculate PPL-like contribution to intermediate V, term 1
C Build F12 integrals
      aodim=nnb*nnb
      ijdim=nal*nbe
      icm=dblalloc(aodim*ijdim)
      pq=icm
       do j=1,nbe
         do i=1,nal
          call dgemm('n','t',nnb,nnb,dfnbasis,3d0/8d0, vifaa(1,1,i),
     $nbfcn,cf12bb(1,1,j),nbfcn,0d0,dcore(pq),nnb)
          call dgemm('n','t',nnb,nnb,dfnbasis,3d0/8d0,cf12aa(1,1,i),
     $nbfcn, vifbb(1,1,j),nbfcn,1d0,dcore(pq),nnb)
          call dgemm('n','t',nnb,nnb,dfnbasis,1d0/8d0, vifba(1,1,j),
     $nbfcn,cf12ab(1,1,i),nbfcn,1d0,dcore(pq),nnb)
          call dgemm('n','t',nnb,nnb,dfnbasis,1d0/8d0,cf12ba(1,1,j),
     $nbfcn, vifab(1,1,i),nbfcn,1d0,dcore(pq),nnb)
          pq=pq+aodim
         enddo
       enddo
C Loop over blocks
      maxmem=maxcor-(imem-imem1)
      slo=1
      rs=1
       do sup=1,nb
        rsdim=(sup-slo+1)*nb
         if(sup.eq.nb.or.(rsdim+nb)*aodim.gt.maxmem) then
          iva=dblalloc(aodim*rsdim)
          pq=iva
           do s=slo,sup
             do r=1,nb
              call dgemm('n','t',nnb,nnb,dfnbasis,1d0,vifaa(1,1,r),
     $nbfcn,vifbb(1,1,s),nbfcn,0d0,dcore(pq),nnb)
              pq=pq+aodim
             enddo
           enddo
          call dgemm('t','n',rsdim,ijdim,aodim,-1d0,dcore(iva),
     $aodim,dcore(icm),aodim,1d0,vab(1,slo,1,1),nb*nb)
          call dbldealloc(iva)
          slo=sup+1
          rs=rs+rsdim
         endif
       enddo
      call dbldealloc(icm)
C Calculate PPL-like contribution to intermediate V, term 2
C Build F12 integrals
      aodim=nbfan*(ncore+nbe)
      icm=dblalloc(aodim*ijdim)
      pq=icm
       do j=1,nbe
         do i=1,nal
          call dgemm('n','t',nbfan,ncore+nbe,dfnbasis,3d0/8d0,
     $ vifaa(nnb+1,1,i),nbfcn,cf12bb(1,1,j),nbfcn,0d0,dcore(pq),nbfan)
          call dgemm('n','t',nbfan,ncore+nbe,dfnbasis,3d0/8d0,
     $cf12aa(nnb+1,1,i),nbfcn, vifbb(1,1,j),nbfcn,1d0,dcore(pq),nbfan)
          call dgemm('n','t',nbfan,ncore+nbe,dfnbasis,1d0/8d0,
     $ vifba(nnb+1,1,j),nbfcn,cf12ab(1,1,i),nbfcn,1d0,dcore(pq),nbfan)
          call dgemm('n','t',nbfan,ncore+nbe,dfnbasis,1d0/8d0,
     $cf12ba(nnb+1,1,j),nbfcn, vifab(1,1,i),nbfcn,1d0,dcore(pq),nbfan)
          pq=pq+aodim
         enddo
       enddo
C Loop over blocks
      maxmem=maxcor-(imem-imem1)
      slo=1
      rs=1
       do sup=1,nb
        rsdim=(sup-slo+1)*nb
         if(sup.eq.nb.or.(rsdim+nb)*aodim.gt.maxmem) then
          iva=dblalloc(aodim*rsdim)
          pq=iva
           do s=slo,sup
             do r=1,nb
              call dgemm('n','t',nbfan,ncore+nbe,dfnbasis,1d0,
     $vifaa(nnb+1,1,r),nbfcn,vifbb(1,1,s),nbfcn,0d0,dcore(pq),nbfan)
              pq=pq+aodim
             enddo
           enddo
          call dgemm('t','n',rsdim,ijdim,aodim,-1d0,dcore(iva),
     $aodim,dcore(icm),aodim,1d0,vab(1,slo,1,1),nb*nb)
          call dbldealloc(iva)
          slo=sup+1
          rs=rs+rsdim
         endif
       enddo
      call dbldealloc(icm)
C Calculate PPL-like contribution to intermediate V, term 3
C Build F12 integrals
      aodim=nbfan*(ncore+nal)
      icm=dblalloc(aodim*ijdim)
      pq=icm
       do j=1,nbe
         do i=1,nal
          call dgemm('n','t',ncore+nal,nbfan,dfnbasis,3d0/8d0,
     $ vifaa(1,1,i),nbfcn,cf12bb(nnb+1,1,j),nbfcn,0d0,dcore(pq),
     $ncore+nal)
          call dgemm('n','t',ncore+nal,nbfan,dfnbasis,3d0/8d0,
     $cf12aa(1,1,i),nbfcn, vifbb(nnb+1,1,j),nbfcn,1d0,dcore(pq),
     $ncore+nal)
          call dgemm('n','t',ncore+nal,nbfan,dfnbasis,1d0/8d0,
     $ vifba(1,1,j),nbfcn,cf12ab(nnb+1,1,i),nbfcn,1d0,dcore(pq),
     $ncore+nal)
          call dgemm('n','t',ncore+nal,nbfan,dfnbasis,1d0/8d0,
     $cf12ba(1,1,j),nbfcn, vifab(nnb+1,1,i),nbfcn,1d0,dcore(pq),
     $ncore+nal)
          pq=pq+aodim
         enddo
       enddo
C Loop over blocks
      maxmem=maxcor-(imem-imem1)
      slo=1
      rs=1
       do sup=1,nb
        rsdim=(sup-slo+1)*nb
         if(sup.eq.nb.or.(rsdim+nb)*aodim.gt.maxmem) then
          iva=dblalloc(aodim*rsdim)
          pq=iva
           do s=slo,sup
             do r=1,nb
              call dgemm('n','t',ncore+nal,nbfan,dfnbasis,1d0,
     $vifaa(1,1,r),nbfcn,vifbb(nnb+1,1,s),nbfcn,0d0,dcore(pq),ncore+nal)
              pq=pq+aodim
             enddo
           enddo
          call dgemm('t','n',rsdim,ijdim,aodim,-1d0,dcore(iva),
     $aodim,dcore(icm),aodim,1d0,vab(1,slo,1,1),nb*nb)
          call dbldealloc(iva)
          slo=sup+1
          rs=rs+rsdim
         endif
       enddo
      call dbldealloc(icm)
C
      return
      end
C
************************************************************************
      subroutine uaiaacalc(nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,vifaa,
     $cf12aa,focka,f12ao,f12ai,vpia,uaia,caia)
************************************************************************
* Calculate Uai and Cai intermediates, same-spin block
************************************************************************
      implicit none
      integer nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,i,j
      real*8 vpia(nb,nal),uaia(nval,nal),caia(nval,nal)
      real*8 f12ao(nbfan,nval),f12ai(nbfan,nval),focka(nbfcn,nbfcn)
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
C
      uaia(:,:)=vpia(nal+1:nal+nval,:)
      caia=0d0
       do j=1,nal
         do i=1,j-1
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $ vifaa(ncore+nb+1,1,i),nbfcn,cf12aa(ncore+nal+1,1,j),nbfcn,0d0,
     $f12ao,nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa(ncore+nb+1,1,i),nbfcn, vifaa(ncore+nal+1,1,j),nbfcn,1d0,
     $f12ao,nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,-1d0,
     $ vifaa(ncore+nb+1,1,j),nbfcn,cf12aa(ncore+nal+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
          call dgemm('n','t',nbfan,nval,dfnbasis,-1d0,
     $cf12aa(ncore+nb+1,1,j),nbfcn, vifaa(ncore+nal+1,1,i),nbfcn,1d0,
     $f12ao,nbfan)
          call dscal(nbfan*nval,0.25d0,f12ao,1)
          call dgemv('t',nbfan,nval, 1d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+i),1,1d0,caia(1,j),1)
          call dgemv('t',nbfan,nval,-1d0,f12ao,nbfan,
     $focka(ncore+nb+1,ncore+j),1,1d0,caia(1,i),1)
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $vifaa(ncore+nb+1,1,i),nbfcn,vifaa(ncore+1,1,j),nbfcn,0d0,f12ai,
     $nbfan)
          call dgemm('n','t',nbfan,nal,dfnbasis,-1d0,
     $vifaa(ncore+nb+1,1,j),nbfcn,vifaa(ncore+1,1,i),nbfcn,1d0,f12ai,
     $nbfan)
          call dgemm('t','n',nval,nal,nbfan,-2d0,f12ao,nbfan,f12ai,
     $nbfan,1d0,uaia,nval)
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine cabijal(nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,
     $dfnbasis,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,
     $f12caj,f12abb,f12acj,f12abk,vinabk,vinabb,vinaab,f12ao,focka,
     $fockb,cabijaa,cabijbb,cabijab,f12pq,f12rs,vaa,vbb,vab)
************************************************************************
* Calculate Cabij intermediate, open-shell, alpha loop out
************************************************************************
      implicit none
      integer nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis
      integer i,j,k,ij,jk,kl,a,b,c,ab
      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn),f12ao(nbfan,nval),tmp
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
      real*8 vifbb(nbfcn,dfnbasis,nb),cf12bb(nbfcn,dfnbasis,nbe)
      real*8 vifab(nbfcn,dfnbasis,nb),cf12ab(nbfcn,dfnbasis,nal)
      real*8 vifba(nbfcn,dfnbasis,nb),cf12ba(nbfcn,dfnbasis,nbe)
      real*8 f12caj(nbfan,nval,nal),f12acj(nbfan,nval,nal)
      real*8 f12abk(nbfan,nval,nal),vinabk(nval,nbfan,nal)
      real*8 vinabb(nvbe,nbfan,nbe),vinaab(nvbe,nbfan,nal)
      real*8 cabijaa((nval-1)*nval/2,(nal-1)*nal/2),f12pq(nval,nval)
      real*8 cabijbb((nvbe-1)*nvbe/2,(nbe-1)*nbe/2),f12rs(nvbe,nvbe)
      real*8 cabijab(nval,nvbe,nal,nbe),f12abb(nbfan,nvbe,nbe)
      real*8 vaa((nb-1)*nb/2,(nal-1)*nal/2),vab(nb,nb,nal,nbe)
      real*8 vbb((nb-1)*nb/2,(nbe-1)*nbe/2)
C Initialize C
      ab=0
       do b=1,nval
         do a=1,b-1
          ab=ab+1
          cabijaa(ab,:)=vaa((nal+b-1)*(nal+b-2)/2+nal+a,:)
         enddo
       enddo
      ab=0
       do b=1,nvbe
         do a=1,b-1
          ab=ab+1
          cabijbb(ab,:)=vbb((nbe+b-1)*(nbe+b-2)/2+nbe+a,:)
         enddo
       enddo
      cabijab(:,:,:,:)=vab(nal+1:nal+nval,nbe+1:nbe+nvbe,:,:)
C Alpha loop
      jk=0
       do k=1,nal
        call cabijss(k,jk,nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,vifaa,
     $cf12aa,f12caj,f12acj,f12abk,vinabk,f12ao,focka,cabijaa,f12pq)
         do j=1,nbe
          call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nb+1,1,k),nbfcn,cf12bb(ncore+nbe+1,1,j),nbfcn,0d0,
     $f12abb(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nb+1,1,k),nbfcn, vifbb(ncore+nbe+1,1,j),nbfcn,1d0,
     $f12abb(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $ vifba(ncore+nb+1,1,j),nbfcn,cf12ab(ncore+nbe+1,1,k),nbfcn,1d0,
     $f12abb(1,1,j),nbfan)
          call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nb+1,1,j),nbfcn, vifab(ncore+nbe+1,1,k),nbfcn,1d0,
     $f12abb(1,1,j),nbfan)
           do i=1,nal
            call dgemm('n','n',nval,nvbe,nbfan,1d0,
     $vinabk(1,1,i),nval,f12abb(1,1,j),nbfan,1d0,cabijab(1,1,i,j),nval)
           enddo
         enddo
         do a=1,nvbe
          call dgemm('n','t',nbfan,nbe,dfnbasis, 1d0,
     $vifaa(ncore+nb+1,1,k),nbfcn,vifbb(ncore+1,1,nbe+a),nbfcn,0d0,
     $f12ao,nbfan)
          call dcopy(nbfan*nbe,f12ao,1,vinabb(a,1,1),nvbe)
         enddo
        ij=0
         do j=1,nbe
           do i=1,j-1
            ij=ij+1
            call dgemm('n','n',nvbe,nvbe,nbfan, 1d0,
     $vinabb(1,1,i),nvbe,f12abb(1,1,j),nbfan,0d0,f12rs,nvbe)
            call dgemm('n','n',nvbe,nvbe,nbfan,-1d0,
     $vinabb(1,1,j),nvbe,f12abb(1,1,i),nbfan,1d0,f12rs,nvbe)
            ab=0
             do b=1,nvbe
               do a=1,b-1
                ab=ab+1
                cabijbb(ab,ij)=cabijbb(ab,ij)+f12rs(a,b)-f12rs(b,a)
               enddo
             enddo
           enddo
         enddo
         do j=1,nbe
           do i=1,nal
            call dgemm('t','t',nval,nvbe,nbfan,1d0,
     $f12abk(1,1,i),nbfan,vinabb(1,1,j),nvbe,1d0,cabijab(1,1,i,j),nval)
           enddo
         enddo
         do a=1,nvbe
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $vifbb(ncore+nb+1,1,nbe+a),nbfcn,vifaa(ncore+1,1,k),nbfcn,0d0,
     $f12ao,nbfan)
          call dcopy(nbfan*nal,f12ao,1,vinaab(a,1,1),nvbe)
         enddo
         do j=1,nbe
          call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nal+1,1,k),nbfcn,cf12bb(ncore+nb+1,1,j),nbfcn,0d0,
     $f12ao,nval)
          call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nal+1,1,k),nbfcn, vifbb(ncore+nb+1,1,j),nbfcn,1d0,
     $f12ao,nval)
          call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0,
     $ vifba(ncore+nal+1,1,j),nbfcn,cf12ab(ncore+nb+1,1,k),nbfcn,1d0,
     $f12ao,nval)
          call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nal+1,1,j),nbfcn, vifab(ncore+nb+1,1,k),nbfcn,1d0,
     $f12ao,nval)
           do i=1,nal
            call dgemm('n','t',nval,nvbe,nbfan,-1d0,
     $f12ao,nval,vinaab(1,1,i),nvbe,1d0,cabijab(1,1,i,j),nval)
           enddo
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine cabijss(k,jk,nb,ncore,nal,nval,nbfcn,nbfan,dfnbasis,
     $vifaa,cf12aa,f12caj,f12acj,f12abk,vinabk,f12ao,focka,cabijaa,
     $f12pq)
************************************************************************
* Calculate Cabij intermediate, same-spin contribution
************************************************************************
      implicit none
      integer nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis
      integer i,j,k,ij,jk,kl,a,b,c,ab
      real*8 focka(nbfcn,nbfcn),f12ao(nbfan,nval),f12pq(nval,nval),tmp
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
      real*8 f12caj(nbfan,nval,nal),f12acj(nbfan,nval,nal)
      real*8 f12abk(nbfan,nval,nal),vinabk(nval,nbfan,nal)
      real*8 cabijaa((nval-1)*nval/2,(nal-1)*nal/2)
C
       do j=1,nal
        call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $ vifaa(ncore+nb+1,1,j),nbfcn,cf12aa(ncore+nal+1,1,k),nbfcn,0d0,
     $f12caj(1,1,j),nbfan)
        call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa(ncore+nb+1,1,j),nbfcn, vifaa(ncore+nal+1,1,k),nbfcn,1d0,
     $f12caj(1,1,j),nbfan)
       enddo
      call dscal(nbfan*nval*nal,1d0/4d0,f12caj,1)
       do j=1,nal
        call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $ vifaa(ncore+nb+1,1,k),nbfcn,cf12aa(ncore+nal+1,1,j),nbfcn,0d0,
     $f12acj(1,1,j),nbfan)
        call dgemm('n','t',nbfan,nval,dfnbasis,1d0,
     $cf12aa(ncore+nb+1,1,k),nbfcn, vifaa(ncore+nal+1,1,j),nbfcn,1d0,
     $f12acj(1,1,j),nbfan)
       enddo
      call dscal(nbfan*nval*nal,1d0/4d0,f12acj,1)
      call dcopy(nbfan*nval*nal,f12acj,1,f12abk,1)
      f12abk=-f12abk
      call daxpy(nbfan*nval*nal,1d0,f12caj,1,f12abk,1)
       do j=1,k-1
        jk=jk+1
        call dgemm('n','n',nval,nval,nbfan, 1d0,
     $focka(ncore+nal+1,ncore+nb+1),nbfcn,f12abk(1,1,j),nbfan,0d0,
     $f12pq,nval)
        ab=0
         do b=1,nval
           do a=1,b-1
            ab=ab+1
            cabijaa(ab,jk)=cabijaa(ab,jk)+f12pq(a,b)-f12pq(b,a)
           enddo
         enddo
       enddo
      call dcopy(nbfan*nval*nal,f12caj,1,f12abk,1)
      f12abk=-f12abk
      call daxpy(nbfan*nval*nal,1d0,f12acj,1,f12abk,1)
       do a=1,nval
        call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $vifaa(ncore+nb+1,1,k),nbfcn,vifaa(ncore+1,1,nal+a),nbfcn,0d0,
     $f12ao,nbfan)
        call dgemm('n','t',nbfan,nal,dfnbasis,-1d0,
     $vifaa(ncore+nb+1,1,nal+a),nbfcn,vifaa(ncore+1,1,k),nbfcn,1d0,
     $f12ao,nbfan)
        call dcopy(nbfan*nal,f12ao,1,vinabk(a,1,1),nval)
       enddo
      ij=0
       do j=1,nal
         do i=1,j-1
          ij=ij+1
          call dgemm('n','n',nval,nval,nbfan, 1d0,
     $vinabk(1,1,i),nval,f12abk(1,1,j),nbfan,0d0,f12pq,nval)
          call dgemm('n','n',nval,nval,nbfan,-1d0,
     $vinabk(1,1,j),nval,f12abk(1,1,i),nbfan,1d0,f12pq,nval)
          ab=0
           do b=1,nval
             do a=1,b-1
              ab=ab+1
              cabijaa(ab,ij)=cabijaa(ab,ij)+f12pq(a,b)-f12pq(b,a)
             enddo
           enddo
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine cabijbe(nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,
     $dfnbasis,vifaa,vifbb,vifab,vifba,cf12aa,cf12bb,cf12ab,cf12ba,
     $f12caj,f12abb,f12acj,f12abk,vinabk,vinabb,vinaab,f12ao,focka,
     $fockb,cabijaa,cabijbb,cabijab,f12pq,f12rs)
************************************************************************
* Calculate Cabij intermediate, open-shell, beta loop out
************************************************************************
      implicit none
      integer nb,ncore,nal,nbe,nval,nvbe,nbfcn,nbfan,dfnbasis
      integer i,j,k,ij,jk,kl,a,b,c,ab
      real*8 focka(nbfcn,nbfcn),fockb(nbfcn,nbfcn),f12ao(nbfan,nval),tmp
      real*8 vifaa(nbfcn,dfnbasis,nb),cf12aa(nbfcn,dfnbasis,nal)
      real*8 vifbb(nbfcn,dfnbasis,nb),cf12bb(nbfcn,dfnbasis,nbe)
      real*8 vifab(nbfcn,dfnbasis,nb),cf12ab(nbfcn,dfnbasis,nal)
      real*8 vifba(nbfcn,dfnbasis,nb),cf12ba(nbfcn,dfnbasis,nbe)
      real*8 f12caj(nbfan,nvbe,nbe),f12acj(nbfan,nvbe,nbe)
      real*8 f12abk(nbfan,nvbe,nbe),vinabk(nvbe,nbfan,nbe)
      real*8 vinabb(nval,nbfan,nal),vinaab(nval,nbfan,nal)
      real*8 cabijaa((nval-1)*nval/2,(nal-1)*nal/2),f12pq(nvbe,nvbe)
      real*8 cabijbb((nvbe-1)*nvbe/2,(nbe-1)*nbe/2),f12rs(nval,nval)
      real*8 cabijab(nval,nvbe,nal,nbe),f12abb(nval,nbfan,nal)
C Beta loop
      jk=0
       do k=1,nbe
        call cabijss(k,jk,nb,ncore,nbe,nvbe,nbfcn,nbfan,dfnbasis,vifbb,
     $cf12bb,f12caj,f12acj,f12abk,vinabk,f12ao,fockb,cabijbb,f12pq)
         do j=1,nal
          call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nal+1,1,j),nbfcn,cf12bb(ncore+nb+1,1,k),nbfcn,0d0,
     $f12abb(1,1,j),nval)
          call dgemm('n','t',nval,nbfan,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nal+1,1,j),nbfcn, vifbb(ncore+nb+1,1,k),nbfcn,1d0,
     $f12abb(1,1,j),nval)
          call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0,
     $ vifba(ncore+nal+1,1,k),nbfcn,cf12ab(ncore+nb+1,1,j),nbfcn,1d0,
     $f12abb(1,1,j),nval)
          call dgemm('n','t',nval,nbfan,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nal+1,1,k),nbfcn, vifab(ncore+nb+1,1,j),nbfcn,1d0,
     $f12abb(1,1,j),nval)
         enddo
         do j=1,nal
          call dgemm('n','n',nval,nvbe,nbfan,1d0,f12abb(1,1,j),nval,
     $fockb(ncore+nb+1,ncore+nbe+1),nbfcn,1d0,cabijab(1,1,j,k),nval)
           do i=1,nal
            call dgemm('n','t',nval,nvbe,nbfan,1d0,
     $f12abb(1,1,i),nval,vinabk(1,1,j),nvbe,1d0,cabijab(1,1,i,j),nval)
           enddo
         enddo
         do a=1,nval
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $vifbb(ncore+nb+1,1,k),nbfcn,vifaa(ncore+1,1,nal+a),nbfcn,0d0,
     $f12ao,nbfan)
          call dcopy(nbfan*nal,f12ao,1,vinabb(a,1,1),nval)
         enddo
        ij=0
         do j=1,nal
           do i=1,j-1
            ij=ij+1
            call dgemm('n','t',nval,nval,nbfan, 1d0,
     $vinabb(1,1,i),nval,f12abb(1,1,j),nval,0d0,f12rs,nval)
            call dgemm('n','t',nval,nval,nbfan,-1d0,
     $vinabb(1,1,j),nval,f12abb(1,1,i),nval,1d0,f12rs,nval)
            ab=0
             do b=1,nval
               do a=1,b-1
                ab=ab+1
                cabijaa(ab,ij)=cabijaa(ab,ij)+f12rs(a,b)-f12rs(b,a)
               enddo
             enddo
           enddo
         enddo
         do j=1,nbe
           do i=1,nal
            call dgemm('n','n',nval,nvbe,nbfan,1d0,
     $vinabb(1,1,i),nval,f12abk(1,1,j),nbfan,1d0,cabijab(1,1,i,j),nval)
           enddo
         enddo
         do a=1,nval
          call dgemm('n','t',nbfan,nal,dfnbasis, 1d0,
     $vifaa(ncore+nb+1,1,nal+a),nbfcn,vifbb(ncore+1,1,k),nbfcn,0d0,
     $f12ao,nbfan)
          call dcopy(nbfan*nal,f12ao,1,vinaab(a,1,1),nval)
         enddo
         do i=1,nal
          call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $ vifaa(ncore+nb+1,1,i),nbfcn,cf12bb(ncore+nbe+1,1,k),nbfcn,0d0,
     $f12caj(1,1,i),nbfan)
          call dgemm('n','t',nbfan,nvbe,dfnbasis,3d0/8d0,
     $cf12aa(ncore+nb+1,1,i),nbfcn, vifbb(ncore+nbe+1,1,k),nbfcn,1d0,
     $f12caj(1,1,i),nbfan)
          call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $ vifba(ncore+nb+1,1,k),nbfcn,cf12ab(ncore+nbe+1,1,i),nbfcn,1d0,
     $f12caj(1,1,i),nbfan)
          call dgemm('n','t',nbfan,nvbe,dfnbasis,1d0/8d0,
     $cf12ba(ncore+nb+1,1,k),nbfcn, vifab(ncore+nbe+1,1,i),nbfcn,1d0,
     $f12caj(1,1,i),nbfan)
          call dgemm('n','n',nval,nvbe,nbfan,1d0,
     $focka(ncore+nal+1,ncore+nb+1),nbfcn,f12caj(1,1,i),nbfan,1d0,
     $cabijab(1,1,i,k),nval)
         enddo
         do j=1,nbe
           do i=1,nal
            call dgemm('n','n',nval,nvbe,nbfan,-1d0,
     $vinaab(1,1,j),nval,f12caj(1,1,i),nbfan,1d0,cabijab(1,1,i,j),nval)
           enddo
         enddo
       enddo
C
      return
      end
C
************************************************************************
      subroutine saveimos(nb,nal,nbe,nval,nvbe,scrfile1,cabijaa,
     $cabijbb,cabijab,vaa,vbb,vab,uaia,uaib,vpia,vpib,caia,caib,tscalea,
     $tscaleb,tfact,epaa,epbb,epab)
************************************************************************
* Save intermediates, open-shell
************************************************************************
      implicit none
      integer nb,nal,nbe,nval,nvbe,scrfile1,i,j,k,ij,kl,a,b
      real*8 vaa((nb-1)*nb/2,(nal-1)*nal/2),vab(nb,nb,nal,nbe)
      real*8 vbb((nb-1)*nb/2,(nbe-1)*nbe/2)
      real*8 vpia(nb,nal),vpib(nb,nbe),uaia(nval,nal),uaib(nvbe,nbe)
      real*8 caia(nval,nal),cabijaa((nval-1)*nval/2,(nal-1)*nal/2)
      real*8 caib(nvbe,nbe),cabijbb((nvbe-1)*nvbe/2,(nbe-1)*nbe/2)
      real*8 cabijab(nval,nvbe,nal,nbe),tscalea(nal),tscaleb(nbe),tfact
      real*8 epaa(nal*(nal-1)/2),epbb(nbe*(nbe-1)/2),epab(nal,nbe)
C
      write(scrfile1) tfact,tscalea,tscaleb,epaa,epbb,epab
      if(nal.gt.1.and.nval.gt.1) write(scrfile1) cabijaa
      if(nbe.gt.1.and.nvbe.gt.1) write(scrfile1) cabijbb
      if(nal.gt.0.and.nval.gt.0.and.nbe.gt.0.and.nvbe.gt.0)
     $write(scrfile1) cabijab
      if(nal.gt.1.and.nval.gt.1) write(scrfile1)
     $(((vaa((b-1)*(b-2)/2+a,ij),a=nal+1,b-1),b=nal+1,nb),
     $ij=1,(nal-1)*nal/2)
      if(nbe.gt.1.and.nvbe.gt.1) write(scrfile1)
     $(((vbb((b-1)*(b-2)/2+a,ij),a=nbe+1,b-1),b=nbe+1,nb),
     $ij=1,(nbe-1)*nbe/2)
      if(nal.gt.0.and.nval.gt.0.and.nbe.gt.0.and.nvbe.gt.0)
     $write(scrfile1) (((vab(nal+1:nb,nbe+b,i,j),b=1,nvbe),i=1,nal),
     $j=1,nbe)  ! before: vab(nal+1:nb,nbe+1:nb,:,:)
      if(nal.gt.0.and.nval.gt.0) write(scrfile1) uaia
      if(nbe.gt.0.and.nvbe.gt.0) write(scrfile1) uaib
      close(scrfile1)
      open(scrfile1,file='F12INT1',form='unformatted')
      if(nal.gt.0.and.nval.gt.0) then
         write(scrfile1) vpia(1:nal,:)
         write(scrfile1) caia
      endif
      if(nbe.gt.0.and.nvbe.gt.0) then
         write(scrfile1) vpib(1:nbe,:)
         write(scrfile1) caib
      endif
      close(scrfile1)
      open(scrfile1,file='F12INT2',form='unformatted')
      if(nal.gt.1.and.nval.gt.0) write(scrfile1)
     $(((vaa((a-1)*(a-2)/2+k,ij),k=1,nal),a=nal+1,nb),
     $ij=1,(nal-1)*nal/2)
      if(nal.gt.1.and.nval.gt.1) write(scrfile1)
     $((vaa(kl,ij),kl=1,(nal-1)*nal/2),ij=1,(nal-1)*nal/2)
      if(nbe.gt.1.and.nvbe.gt.0) write(scrfile1)
     $(((vbb((a-1)*(a-2)/2+k,ij),k=1,nbe),a=nbe+1,nb),
     $ij=1,(nbe-1)*nbe/2)
      if(nbe.gt.1.and.nvbe.gt.1) write(scrfile1)
     $((vbb(kl,ij),kl=1,(nbe-1)*nbe/2),ij=1,(nbe-1)*nbe/2)
      if(nal.gt.0.and.nval.gt.0.and.nbe.gt.0) write(scrfile1)
     $((((vab(nal+a,b,i,j),b=1,nbe),a=1,nval),i=1,nal),j=1,nbe)
      if(nal.gt.0.and.nvbe.gt.0.and.nbe.gt.0) write(scrfile1)
     $(((vab(1:nal,nbe+b,i,j),b=1,nvbe),i=1,nal),j=1,nbe) ! before vab(1:nal,nbe+1:nb,:,:)
      if(nal.gt.0.and.nval.gt.0.and.nbe.gt.0.and.nvbe.gt.0)
     $write(scrfile1) (((vab(1:nal,b,i,j),b=1,nbe),i=1,nal),j=1,nbe) ! before vab(1:nal,1:nbe,:,:)
C
      return
      end
C
************************************************************************
      subroutine saveimcs(nb,nal,nval,scrfile1,cabij,vv,uaia,vpia,caia,
     $tscale,tfact,epaa)
************************************************************************
* Save intermediates, closed-shell
************************************************************************
      implicit none
      integer nb,nal,nval,scrfile1,i,j,k,ij,a,b
      real*8 vv(nb,nb,(nal+1)*nal/2),cabij(nval,nval,(nal+1)*nal/2)
      real*8 vpia(nb,nal),uaia(nval,nal),caia(nval,nal),tscale(nal)
      real*8 tfact,epaa(nal*(nal+1)/2)
C
      write(scrfile1) tfact,tscale,epaa
      write(scrfile1)
     $((((2.d0*cabij(a,b,j*(j-1)/2+i)-cabij(b,a,j*(j-1)/2+i),
     $a=1,nval),b=1,nval),i=1,j-1),
     $((0.5d0*cabij(a,b,j*(j-1)/2+j),a=1,nval),b=1,nval),j=1,nal)
      write(scrfile1)
     $((((2.d0*vv(nal+a,nal+b,j*(j-1)/2+i)-vv(nal+b,nal+a,j*(j-1)/2+i),
     $a=1,nval),b=1,nval),i=1,j-1),
     $((0.5d0*vv(nal+a,nal+b,j*(j-1)/2+j),a=1,nval),b=1,nval),j=1,nal)
      write(scrfile1) uaia
      close(scrfile1)
      open(scrfile1,file='F12INT1',form='unformatted')
      write(scrfile1) caia,
     $((((cabij(a,b,j*(j-1)/2+i),a=1,nval),b=1,nval),i=1,j),j=1,nal)
      ij=0
       do j=1,nal
         do i=1,j
          ij=ij+1
          write(scrfile1) transpose(vv(1:nal,1:nal,ij))
         enddo
       enddo
       write(scrfile1) vpia(1:nal,:)
       do j=1,nal
        write(scrfile1) (vv(1:nal,nal+1:nb,j*(j-1)/2+i),i=1,j),
     $          (transpose(vv(nal+1:nb,1:nal,i*(i-1)/2+j)),i=j+1,nal)
       enddo
C
      return
      end
C
