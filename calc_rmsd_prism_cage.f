        program calc_rmsd_prism_cage
        use omp_lib
        implicit real*8(a-h,o-z)
        parameter (nwavetot=50)
        parameter (npermut=720)
        parameter (nmax=6000000)
        parameter (natoms=18)
        parameter (ntmax=15)
        parameter (noxy=6)
        parameter (ndim=54)
        parameter (nisomax=100)
        parameter (nbin=100)
        dimension iopt(nmax),roo_use(ntmax,npermut,nisomax),
     1  pair(ntmax,nmax),psips(ndim,nmax),roo(ntmax,npermut,nisomax),
     1  x(3,natoms,nmax),pair_cage(ntmax,nmax),pair_prism(ntmax,nmax),
     1  weight_tot(nmax),rmsd_cage(nisomax,nmax),psips_reject(ndim,nmax)
     1  ,rmsd_prism(nisomax,nmax),psips_tot(ndim,nmax),weight(nmax),
     1  weight_reject(nmax),num_iso(nmax),min_rmsd_cage(nmax),
     1  min_rmsd_prism(nmax),rval1(nbin**2),box_prism(nbin**2),
     1  box_cage(nbin**2),diff_cage(nmax),diff_prism(nmax),
     1  weight_cage(nmax),weight_prism(nmax),coord_ref(3,natoms,nisomax)
     1  ,diff_other(nmax),tot_box_cage(nbin**2),tot_box_prism(nbin**2),
     1  pair_other(ntmax,nmax),weight_other(nmax),min_rmsd_other(nmax),
     1  rmsd_other(nisomax,nmax),rval2(nbin**2),tot_box_other(nbin**2),
     1  box_other(nbin**2),cage_min(nmax),prism_min(nmax),
     1  other_min(nmax),psips_cage(ndim,nmax),psips_prism(ndim,nmax),
     1  psips_other(ndim,nmax),psips_bad_cage(ndim,nmax),min_rmsd(nmax),
     1  psips_bad_prism(ndim,nmax),rmsd(nisomax,nmax),diff(nmax)
        common/coo/coord(3,natoms,nmax)
        character(len=1024) file_iso,wf,wt,iso_file,iso_psips,isomer,
     1  comp_iso
        open(unit=7,file='wf_final.dat',status='old',form='unformatted')
        open(unit=11,file='cage_prism.dat',status='unknown')
        open(unit=12,file='cage_prism-log.dat',status='unknown')
        open(unit=13,file='isomer_coordinates.dat',status='old')
C		Given a DMC simulation that contains prisms and cages, separate the different 
C		structures to classify each of them as prisms or cages based on the minimum 
C		energy strucutre. Final output is a 2D histogram of the results. 
C		Inputs:
C		Geometries of the walkers in bohr
C		Minimum energy structures of the cage and prism.
C		Outputs:
C		2D histogram file of the (\rho(prism)-\rho(cage)) both in non-log10 and log10 
c       read coordinates of all molecules 
        do i = 1,nbin**2
            tot_box_cage(i) = 0.d0
            tot_box_prism(i) = 0.d0
            tot_box_other(i) = 0.d0
        enddo
        niter1 = 0
        niter2 = 0
        niter3 = 0
        read(7) nwave
        do k = 1,nwave
            read(7) n
            tot_weight = 0.d0
            psips = 0.d0
            do i = 1,n
                read(7) (psips(j,i),j=1,ndim),num_iso(i),weight(i)
                iopt(i) = 1
                tot_weight = tot_weight + weight(i)
            enddo
            iiso = 2
            do i = 1,iiso
                do j = 1,natoms
                    read(13,*) (x(l,j,i),l=1,3)
                enddo
            enddo
            rewind(13)
            do i = 1,iiso
                do j = 1,natoms
                    do l = 1,3
                        x(l,j,i) = x(l,j,i)/0.52917721067d0
                    enddo
                enddo
            enddo
            ntot = 0
            nreject = 0
            wt_reject = 0.d0
            do i = 1,n
                if (iopt(i).eq.1) then
                    ntot = ntot + 1
                    weight_tot(ntot) = weight(i)
                    do j = 1,ndim
                        psips_tot(j,ntot) = psips(j,i)
                    enddo
                else
                    nreject = nreject + 1
                    weight_reject(nreject) = weight(i)
                    wt_reject = wt_reject + weight_reject(nreject)
                    do j = 1,ndim
                        psips_reject(j,nreject) = psips(j,i)
                    enddo
                endif
            enddo
            do i = 1,iiso
                ip = 0
                do j = 1,natoms
                    do l = 1,3
                        ip = ip + 1
                        coord_ref(l,j,i) = x(l,j,i)
                    enddo
                enddo
            enddo
            do i = 1,ntot
                ip = 0
                do j = 1,natoms
                    do l = 1,3
                        ip = ip + 1
                        coord(l,j,i) = psips_tot(ip,i)
                    enddo
                enddo
            enddo
            roo = 0.d0
c       calculate all oo distances for 720 permutations
            do i = 1,iiso
                call calc_ref_pairwise(coord_ref(:,:,i),roo(:,:,i))
                roo_use(:,:,i) = roo(:,:,i)
            enddo
            pair = 0.d0
c       calculate pairwise distances of all 9 structures
            call calc_pairwise(ntot,pair)
c       calculate the rmsd of all 720 permutations for all the pairwise 
            icage = 0
            iprism = 0
            do i = 1,ntot
                do j = 1,iiso
                    call comp_pair(j,pair(:,i),roo_use(:,:,j),ntmax,
     1              rmsd(j,i),i,min_rmsd(i))
                enddo
            enddo
            do i = 1,ntot
                diff(i) = rmsd(2,i)-rmsd(1,i)
            enddo
            do i = 1,ntot
                if (diff(i).gt.0.d0) then
                    icage = icage + 1
                    diff_cage(icage) = diff(i)
                    weight_cage(icage) = weight_tot(i)
                    do j = 1,iiso
                        rmsd_cage(j,icage) = rmsd(j,i)
                    enddo
                else if (diff(i).lt.0.d0) then
                    iprism = iprism + 1
                    diff_prism(iprism) = diff(i)
                    weight_prism(iprism) = weight_tot(i)
                    do j = 1,iiso
                        rmsd_prism(j,iprism) = rmsd(j,i)
                    enddo
                endif
            enddo
            do i = 1,icage
                do j = 1,iiso
                    if (j.eq.1) then
                        cage_min(i) = rmsd_cage(j,i)
                    else 
                        if (rmsd_cage(j,i).lt.cage_min(i)) then
                            cage_min(i) = rmsd_cage(j,i)
                        endif
                    endif
                enddo
            enddo
            do i = 1,iprism
                do j = 1,iiso
                    if (j.eq.1) then
                        prism_min(i) = rmsd_prism(j,i)
                    else 
                        if (rmsd_prism(j,i).lt.prism_min(i)) then
                            prism_min(i) = rmsd_prism(j,i)
                        endif
                    endif
                enddo
            enddo
            if (icage.ne.0) then
                niter1 = niter1 + 1
                call descendant_weight_1(diff_cage(:),cage_min(:),
     1          weight_cage(:),icage,rval1(:),rval2(:),box_cage(:),dx1,
     1          dx2,nbin)
            endif
            if (iprism.ne.0) then
                niter2 = niter2 + 1
                call descendant_weight_1(diff_prism(:),prism_min(:),
     1          weight_prism(:),iprism,rval1(:),rval2(:),box_prism(:),
     1          dx1,dx2,nbin)
            endif
            do i = 1,nbin**2
                tot_box_cage(i) = tot_box_cage(i) + box_cage(i)
                tot_box_prism(i) = tot_box_prism(i) + box_prism(i)
                tot_box_other(i) = tot_box_other(i) + box_other(i)
            enddo
            do i = 1,nbin**2
                box_cage(i) = 0.d0
                box_prism(i) = 0.d0
            enddo
        enddo
        ir = 0
        do i = 1,nbin
            do j = 1,nbin
                ir = ir + 1
                write(11,*) rval1(i),rval2(j),(tot_box_cage(ir)/
     1          dfloat(niter1)),(tot_box_prism(ir)/dfloat(niter2))
                if (((tot_box_cage(ir)/dfloat(niter1)).ne.0.d0).and.
     1          ((tot_box_prism(ir)/dfloat(niter2)).ne.0.d0)) then
                    write(12,*)rval1(i),rval2(j),
     1              log10(tot_box_cage(ir)/dfloat(niter1)),
     1              log10(tot_box_prism(ir)/dfloat(niter2))
                else if (((tot_box_cage(ir)/dfloat(niter1)).eq.0.d0)
     1          .and.((tot_box_prism(ir)/dfloat(niter2)).ne.0.d0)) then
                    write(12,*)rval1(i),rval2(j),-7.d0,
     1              log10(tot_box_prism(ir)/dfloat(niter2))
                else if (((tot_box_cage(ir)/dfloat(niter1)).ne.0.d0)
     1          .and.((tot_box_prism(ir)/dfloat(niter2)).eq.0.d0)) then 
                    write(12,*)rval1(i),rval2(j),
     1              log10(tot_box_cage(ir)/dfloat(niter1)),-7.d0
                else
                    write(12,*)rval1(i),rval2(j),-7.d0,-7.d0
                endif
            enddo
        enddo
        stop
c       loop over all walkers

        end program

        subroutine calc_pairwise(n,pairwise)
        implicit real*8(a-h,o-z)
        parameter (nmax=6000000)
        parameter (ndim=54)
        parameter (np=15)
        parameter (natoms=18)
        parameter (noxy=6)
        common/coo/coord(3,natoms,nmax)
        dimension pairwise(np,nmax),pair(noxy,noxy,nmax),
     1  coord_eq(3,nmax),oxy(3,noxy,nmax)
C		Calculates all of the possible OO distances in the water cluster
C		Inputs:
C		n = number of walkers
C		coord = coordinates of the walker in bohr
C		Output
C		pairwise = all of the OO distances in the water cluster
c       separate psips into arrays regarding water monomers
        do i = 1,n
            do j = 1,3
                coord_eq(j,i) = coord(j,1,i)
            enddo
        enddo
        do i = 1,n
            do j = 1,natoms
                do k = 1,3
                    coord(k,j,i) = coord(k,j,i)-coord_eq(k,i)
                enddo
            enddo
        enddo
c       find which one is the oxygen
        do i = 1,n
            it = 0
            do j = 1,natoms
                if (mod(j,3).eq.1) then
                    it = it + 1
                    do k = 1,3
                        oxy(k,it,i) = coord(k,j,i)
                    enddo
                endif
            enddo
        enddo
c       calculate all the pairwise distances here
        do i = 1,n
            ip = 0
            do j = 1,noxy
                do k = j+1,noxy
                    ip = ip + 1
                    pair(k,j,i) = sqrt((((oxy(1,j,i)-oxy(1,k,i))**2))+
     1              ((oxy(2,j,i)-oxy(2,k,i))**2)+((oxy(3,j,i)-
     1              oxy(3,k,i))**2))
                    pairwise(ip,i) = pair(k,j,i)
                enddo
            enddo
        enddo
        return
        end subroutine

        subroutine comp_pair(n,pair,pair_ref,npair,rmsd_final,it,
     1  min_rmsd)
        implicit real*8(a-h,o-z)
        parameter (nmax=6000000)
        parameter (np=15)
        parameter (npmax=720)
        parameter (noxy=6)
        dimension pair(np,nmax),pair_ref(np,npmax),rmsd_final(nmax),
     1  rmsd(npmax,nmax),rmsd_slice(npmax)
C		Calculate the differences in the OO distances chosen from the 
C		water clusters and comparing to the reference structure.
C		Considers all of the possible 720 permutations of the water
C		hexamer.
C		Inputs:
C		n = number of walkers
C		pair = OO distances of the walker
C		pair_ref = OO distance of the reference structure
C		npair = number of pairwise interactions
C		Outputs:
C		rmsd_final = RMSD comparison between the 2 strucutres
C		it = nth permutation that provides the minimum RMSD.
C		min_rmsd = Minimum value compared with for prism/cage
C		classification.
c       calculate magnitudes of the distances in the pairs
        do i = 1,n
            do j = 1,npmax
                rmsd(j,i) = 0.d0
            enddo
        enddo
        do i = 1,n
            do j = 1,npmax
                do k =1,npair
                    rmsd(j,i)=rmsd(j,i)+((pair(k,i)-pair_ref(k,j))**2)
                enddo
            enddo
        enddo 
        do i=1,n
            do j = 1,npmax
                rmsd(j,i)=(sqrt(2.d0)/sqrt(dfloat(noxy)*
     1          (dfloat(noxy-1))))*sqrt(rmsd(j,i))
            enddo
        enddo
c       find the smallest one for the rmsd
        do i = 1,n
            rmsd_slice(:) = rmsd(:,i)
            rmsd_final(i) = minval(rmsd_slice(:),1)
            if (it.eq.i) then
                min_rmsd = minloc(rmsd_slice(:),1)
            endif
        enddo
        return
        end subroutine

        subroutine calc_ref_pairwise(x,pairwise)
        implicit real*8(a-h,o-z)
        parameter (natoms=18)
        parameter (noxy=6)
        parameter (ntot=15)
        parameter (npmax=720)
        parameter (natomw=3)
        dimension x(3,natoms),pairwise(ntot,npmax),ipermut(noxy,npmax),
     1  oxy(3,natoms,npmax),pair(noxy,noxy,npmax),coord(3,natomw,noxy)
     1  ,coord_iso(3,natoms,noxy,npmax),coord_eq(3,npmax),
     1  coord_permut(3,natoms,npmax)
        open(unit=9,file='water_permut.dat',status='old')
C		Reads in the coordinates of the minimum prism and cage structures and 
C		permutes all of the possible descriptions of the atoms to use
C		for comparison.
C		Inputs:
C		x = coordinates of minimum energy structure
C		water_permut = all 720 permutations to change the ordering of the water hexamer.
C		Outputs:
C		pairwise = OO distances for the minimum energy strucutre
        read(9,*) npermut
        do i = 1,npermut
            read(9,*) (ipermut(j,i),j=1,noxy)
        enddo
        close(9)
        ip = 0
        do i = 1,noxy
            do j = 1,natomw
                ip = ip + 1
                do k = 1,3
                    coord(k,j,i) = x(k,ip)
                enddo
            enddo
        enddo
c       permut these structures 720 times
        do i = 1,npermut
            do j =1,noxy
                do k = 1,natomw
                    do l = 1,3
                        coord_iso(l,k,j,i) = coord(l,k,ipermut(j,i))
                    enddo
                enddo
            enddo
        enddo
        do i = 1,npermut
            it = 0
            do j =1,noxy
                do k = 1,natomw
                    it = it + 1
                    do l = 1,3
                        coord_permut(l,it,i) = coord_iso(l,k,j,i)
                    enddo
                enddo
            enddo
        enddo
        do i = 1,npermut
            do j = 1,3
                coord_eq(j,i) = coord_permut(j,1,i)
            enddo
        enddo
        do i = 1,npermut
            do j = 1,natoms
                do k = 1,3
                    coord_permut(k,j,i) = coord_permut(k,j,i) - 
     1              coord_eq(k,i)
                enddo
            enddo
        enddo
        do i = 1,npermut
            ip = 0
            do j = 1,natoms
                if (mod(j,3).eq.1) then
                    ip = ip + 1
                    do k = 1,3
                        oxy(k,ip,i) = coord_permut(k,j,i)
                    enddo
                endif
            enddo
        enddo
        pair = 0.d0
        do i = 1,npermut
            ip = 0
            do j = 1,noxy
                do k = j+1,noxy
                    ip = ip + 1
                    pair(k,j,i) = sqrt((((oxy(1,j,i)-oxy(1,k,i))**2))+
     1              ((oxy(2,j,i)-oxy(2,k,i))**2)+((oxy(3,j,i)-
     1              oxy(3,k,i))**2))
                    pairwise(ip,i) = pair(k,j,i)
                enddo
            enddo
        enddo
        return
        end subroutine
